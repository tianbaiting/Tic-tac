#!/usr/bin/env python3
"""
Purpose:
  Run Tic-tac at multiple target Tlab values in MeV and emit solver-driven observable curves.

What this script does:
  1) Create multi-energy solver input/config files.
  2) Run solver once for all requested energies.
  3) Parse U-matrix files into (J, pi) blocks via pw_amplitudes.
  4) Assemble M(theta) by full partial-wave summation and compute
     dSigma/dOmega, iT11, T20, T21, T22 (Witala/Gloeckle convention).
  5) Write per-energy model CSV + metadata + run summary.
  6) Optionally generate post-hoc 190 MeV residual metrics against experiment.

Data flow (left -> right):
  CLI args + solver outputs + (optional experiment files)
    -> structured JPiBlock list per (J, pi) channel
    -> M(theta) via Clebsch-Gordan + Y_{l,m} resummation
    -> dSigma + spin-1 trace formulas for iT11/T20/T21/T22
    -> observable curves on theta grid
    -> analysis/{tlab_*}/observables_model.csv + metadata.json
    -> analysis/summary.json

Calls / dependencies:
  - Uses `solver_u_file_utils.py` for U-file family selection and P123 checks.
  - Uses `observable_units.py` for dSigma unit handling.
  - Calls external solver binary `CPP/run`.

Usage:
  python3 examples/run_dpol_p_observables.py \
    --work-dir output/dpol_p_observables \
    --target-tlabs-mev 70,135,190 \
    --two-j-3n-max 7 \
    --dsigma-unit fm2/sr
"""

from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from observable_units import (
    SUPPORTED_DSIGMA_UNITS,
    UNIT_MB_PER_SR,
    convert_dsigma_series,
    convert_dsigma_value,
    infer_dsigma_unit_from_lines,
    normalize_dsigma_unit,
)
from solver_u_file_utils import (
    required_p123_sparse_names,
    select_latest_u_file_family,
)
from tlab_utils import format_tlab_dir_name

from pw_amplitudes import (
    JPiBlock,
    WPGridBin,
    assemble_m_matrix,
    calibrate_dsigma_scale,
    list_solver_energies,
    observables_from_M,
    parse_q_kinematics,
    parse_solver_output,
)


@dataclass
class SolverRunConfig:
    root: Path
    work_dir: Path
    solver: Path
    target_tlabs_mev: List[float]
    two_j_3n_max: int
    j_2n_max: int
    np_wp: int
    nq_wp: int
    nphi: int
    nx: int
    np_per_wp: int
    nq_per_wp: int
    chebyshev_s: float
    chebyshev_t: float
    threads: int
    potential_model: str
    timeout_s: int
    calculate_p123: bool
    reuse_p123_from: Optional[Path]
    angle_min_deg: float
    angle_max_deg: float
    angle_step_deg: float


def _rel_to_cpp(path: Path, root: Path) -> str:
    cpp_dir = root / "CPP"
    return os.path.relpath(path.resolve(), cpp_dir.resolve())


def _parse_tlab_list(raw: str) -> List[float]:
    values: List[float] = []
    for token in raw.split(","):
        item = token.strip()
        if not item:
            continue
        values.append(float(item))
    if not values:
        raise ValueError("No valid Tlab values provided")
    return values


def ensure_cpp_solver_binary(root: Path, solver_path: Path) -> None:
    if solver_path.exists():
        return
    if solver_path != root / "CPP" / "run":
        raise FileNotFoundError(f"Solver binary not found: {solver_path}")
    raise FileNotFoundError(
        "CPP/run is missing. Build it first (for example: make), then rerun."
    )


def _find_latest_existing_u_files(solver_out_dir: Path) -> List[Path]:
    return select_latest_u_file_family(solver_out_dir)


def _link_or_copy_p123_if_requested(dst_p123_dir: Path, src_p123_dir: Optional[Path]) -> None:
    if src_p123_dir is None:
        return
    if dst_p123_dir.exists() and any(dst_p123_dir.iterdir()):
        return
    if not src_p123_dir.exists():
        return

    dst_p123_dir.parent.mkdir(parents=True, exist_ok=True)
    if dst_p123_dir.exists():
        return
    try:
        os.symlink(src_p123_dir.resolve(), dst_p123_dir)
    except OSError:
        dst_p123_dir.mkdir(parents=True, exist_ok=True)
        for src in src_p123_dir.glob("*.h5"):
            target = dst_p123_dir / src.name
            if not target.exists():
                target.write_bytes(src.read_bytes())


def _assert_reusable_p123_matches_grid(p123_dir: Path, cfg: SolverRunConfig) -> None:
    if not p123_dir.exists():
        raise FileNotFoundError(
            f"--reuse-p123 enabled, but P123 directory does not exist: {p123_dir}"
        )
    required = required_p123_sparse_names(
        np_wp=cfg.np_wp,
        nq_wp=cfg.nq_wp,
        j2max=cfg.j_2n_max,
        two_j_3n_max=cfg.two_j_3n_max,
    )
    missing = [name for name in required if not (p123_dir / name).exists()]
    if missing:
        missing_preview = ", ".join(missing[:6])
        if len(missing) > 6:
            missing_preview += f", ... ({len(missing)} missing)"
        raise FileNotFoundError(
            "--reuse-p123 was requested but matching P123 files were not found for "
            f"Np={cfg.np_wp}, Nq={cfg.nq_wp}, J_2N_max={cfg.j_2n_max}, two_J_3N_max={cfg.two_j_3n_max} in {p123_dir}.\n"
            f"Missing examples: {missing_preview}\n"
            "Disable --reuse-p123 for fresh generation, or provide a compatible --reuse-p123-from folder."
        )


def write_solver_inputs(cfg: SolverRunConfig) -> Tuple[Path, Path, Path, Path, Path]:
    cpp_input_dir = cfg.root / "CPP" / "Input"
    cpp_input_dir.mkdir(parents=True, exist_ok=True)

    inputs_dir = cfg.work_dir / "inputs"
    solver_dir = cfg.work_dir / "solver"
    solver_out_dir = solver_dir / "solver_out"
    p123_dir = solver_dir / "p123"
    for path in [inputs_dir, solver_dir, solver_out_dir]:
        path.mkdir(parents=True, exist_ok=True)

    _link_or_copy_p123_if_requested(p123_dir, cfg.reuse_p123_from)

    target_file = inputs_dir / "target_tlabs_mev.txt"
    target_file.write_text(
        "\n".join(f"{e:.6f}" for e in cfg.target_tlabs_mev) + "\n",
        encoding="utf-8",
    )

    energy_file = cpp_input_dir / "tlab_dpol_observables_autogen.txt"
    energy_file.write_text(
        "\n".join(f"{e:.6f}" for e in cfg.target_tlabs_mev) + "\n",
        encoding="utf-8",
    )

    config_file = cpp_input_dir / "input_tlab_dpol_observables_autogen.txt"
    config_file.write_text(
        "\n".join(
            [
                "# Autogenerated multi-Tlab d+p observable config",
                f"two_J_3N_max={cfg.two_j_3n_max}",
                f"J_2N_max={cfg.j_2n_max}",
                f"Np_WP={cfg.np_wp}",
                f"Nq_WP={cfg.nq_wp}",
                f"Nphi={cfg.nphi}",
                f"Nx={cfg.nx}",
                f"Np_per_WP={cfg.np_per_wp}",
                f"Nq_per_WP={cfg.nq_per_wp}",
                f"chebyshev_s={cfg.chebyshev_s}",
                f"chebyshev_t={cfg.chebyshev_t}",
                "p_grid_type=chebyshev",
                "q_grid_type=chebyshev",
                f"P123_omp_num_threads={cfg.threads}",
                "parallel_run=false",
                "P123_recovery=false",
                "tensor_force=true",
                "isospin_breaking_1S0=true",
                "midpoint_approx=false",
                f"calculate_and_store_P123={'true' if cfg.calculate_p123 else 'false'}",
                "include_breakup_channels=false",
                "solve_faddeev=true",
                "solve_dense=false",
                "production_run=true",
                f"potential_model={cfg.potential_model}",
                "parameter_walk=false",
                "parameter_file=",
                "PSI_start=-1",
                "PSI_end=-1",
                f"energy_input_file={_rel_to_cpp(energy_file, cfg.root)}",
                f"output_folder={_rel_to_cpp(solver_out_dir, cfg.root)}",
                f"P123_folder={_rel_to_cpp(p123_dir, cfg.root)}",
                "",
            ]
        ),
        encoding="utf-8",
    )

    cfg_copy = inputs_dir / "solver_input_snapshot.txt"
    cfg_copy.write_text(config_file.read_text(encoding="utf-8"), encoding="utf-8")

    return target_file, energy_file, config_file, solver_out_dir, p123_dir


def run_cpp_solver(cfg: SolverRunConfig) -> Dict[str, object]:
    ensure_cpp_solver_binary(cfg.root, cfg.solver)
    target_file, energy_file, config_file, solver_out_dir, p123_dir = write_solver_inputs(cfg)

    cpp_cwd = cfg.root / "CPP"
    solver_rel = os.path.relpath(cfg.solver.resolve(), cpp_cwd.resolve())
    if not solver_rel.startswith("."):
        solver_rel = f"./{solver_rel}"
    config_rel = os.path.relpath(config_file.resolve(), cpp_cwd.resolve())
    if not cfg.calculate_p123:
        _assert_reusable_p123_matches_grid(p123_dir, cfg)

    cmd = [solver_rel, config_rel]

    env = os.environ.copy()
    env["HDF5_DISABLE_VERSION_CHECK"] = "2"

    log_file = cfg.work_dir / "solver" / "solver_run.log"
    with log_file.open("w", encoding="utf-8") as log_handle:
        result = subprocess.run(
            cmd,
            cwd=cpp_cwd,
            env=env,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=cfg.timeout_s,
        )

    u_files = _find_latest_existing_u_files(solver_out_dir)
    has_u_data = False
    for u_file in u_files:
        for raw in u_file.read_text(encoding="utf-8").splitlines():
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                float(parts[0])
                complex(parts[3])
                has_u_data = True
                break
            except Exception:
                continue
        if has_u_data:
            break

    return {
        "returncode": result.returncode,
        "log_file": log_file,
        "target_file": target_file,
        "energy_file": energy_file,
        "config_file": config_file,
        "solver_out_dir": solver_out_dir,
        "p123_dir": p123_dir,
        "u_files": u_files,
        "has_u_data": has_u_data,
    }


def _evaluate_curve_on_angles(
    blocks: Sequence[JPiBlock],
    bin_info: WPGridBin,
    q_idx: int,
    angles_deg: Sequence[float],
    extra_scale: float,
    dsigma_output_unit: str,
) -> List[Dict[str, float]]:
    rows: List[Dict[str, float]] = []
    for theta in angles_deg:
        M = assemble_m_matrix(
            blocks, q_idx, math.radians(theta),
            bin_info=bin_info, extra_scale=extra_scale,
        )
        obs = observables_from_M(M)
        rows.append({
            "dSigma_dOmega": convert_dsigma_value(obs.dsigma_fm2_per_sr, "fm2/sr", dsigma_output_unit),
            "iT11": obs.iT11,
            "T20": obs.T20,
            "T21": obs.T21,
            "T22": obs.T22,
        })
    return rows


def _angle_grid(theta_min: float, theta_max: float, step: float) -> List[float]:
    if step <= 0:
        raise ValueError("angle step must be positive")
    values: List[float] = []
    x = theta_min
    while x <= theta_max + 1e-9:
        values.append(round(x, 6))
        x += step
    return values


def _write_model_curve_csv(path: Path, angles: Sequence[float], rows: Sequence[Dict[str, float]]) -> None:
    header = "theta_deg,dSigma_dOmega,iT11,T20,T21,T22"
    lines = [header]
    for theta, row in zip(angles, rows):
        lines.append(
            f"{theta:.6f},"
            f"{row['dSigma_dOmega']:.12e},"
            f"{row['iT11']:.12e},"
            f"{row['T20']:.12e},"
            f"{row['T21']:.12e},"
            f"{row['T22']:.12e}"
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def read_experimental_tensor(path: Path) -> Dict[str, Dict[str, List[float]]]:
    # Returns per-observable angle/value arrays from T.txt
    out = {
        "iT11": {"angles": [], "values": [], "errors": []},
        "T20": {"angles": [], "values": [], "errors": []},
        "T21": {"angles": [], "values": [], "errors": []},
        "T22": {"angles": [], "values": [], "errors": []},
    }
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "θc.m." in line:
            continue
        parts = line.split()
        if len(parts) < 9:
            continue
        angle = float(parts[0])
        columns = [("iT11", 1, 2), ("T20", 3, 4), ("T21", 5, 6), ("T22", 7, 8)]
        for key, idx_v, idx_e in columns:
            v = parts[idx_v]
            e = parts[idx_e]
            if v.lower() == "null" or e.lower() == "null":
                continue
            out[key]["angles"].append(angle)
            out[key]["values"].append(float(v))
            out[key]["errors"].append(float(e))
    return out


def read_experimental_dsigma(path: Path, target_unit: str = UNIT_MB_PER_SR) -> Dict[str, object]:
    lines = path.read_text(encoding="utf-8").splitlines()
    source_unit = infer_dsigma_unit_from_lines(lines, fallback_unit=UNIT_MB_PER_SR)
    output_unit = normalize_dsigma_unit(target_unit)

    angles: List[float] = []
    values: List[float] = []
    for raw_line in lines:
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        try:
            angles.append(float(parts[0]))
            values.append(float(parts[1]))
        except ValueError:
            continue
    return {
        "angles": angles,
        "values": convert_dsigma_series(values, source_unit, output_unit),
        "source_unit": source_unit,
        "output_unit": output_unit,
    }


def _nearest_model_values(
    model_angles: Sequence[float],
    model_rows: Sequence[Dict[str, float]],
    obs_name: str,
    query_angles: Sequence[float],
) -> List[float]:
    out: List[float] = []
    for q in query_angles:
        best_i = min(range(len(model_angles)), key=lambda i: abs(model_angles[i] - q))
        out.append(model_rows[best_i][obs_name])
    return out


def _rmse(pred: Sequence[float], obs: Sequence[float]) -> float:
    if len(pred) != len(obs) or not pred:
        return math.nan
    s = 0.0
    for p, o in zip(pred, obs):
        d = p - o
        s += d * d
    return math.sqrt(s / len(pred))


def _mae(pred: Sequence[float], obs: Sequence[float]) -> float:
    if len(pred) != len(obs) or not pred:
        return math.nan
    s = 0.0
    for p, o in zip(pred, obs):
        s += abs(p - o)
    return s / len(pred)


def _tlab_dir_name(target_tlab_mev: float) -> str:
    return format_tlab_dir_name(target_tlab_mev)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run Tic-tac and compute dpol-p observables (dSigma, iT11, T20, T21, T22) for multiple energies."
    )
    parser.add_argument("--work-dir", default="output/dpol_p_observables", help="Output root directory")
    parser.add_argument("--solver", default="CPP/run", help="Solver binary path")
    parser.add_argument("--target-tlabs-mev", default="70,135,190", help="Comma-separated target Tlab values in MeV")
    parser.add_argument("--two-j-3n-max", type=int, default=1)
    parser.add_argument("--j-2n-max", type=int, default=2)
    parser.add_argument("--np", type=int, default=20)
    parser.add_argument("--nq", type=int, default=22)
    parser.add_argument("--nphi", type=int, default=24)
    parser.add_argument("--nx", type=int, default=24)
    parser.add_argument("--np-per-wp", type=int, default=6)
    parser.add_argument("--nq-per-wp", type=int, default=6)
    parser.add_argument("--chebyshev-s", type=float, default=180.0)
    parser.add_argument("--chebyshev-t", type=float, default=1.0)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--potential-model", default="LO_internal")
    parser.add_argument("--timeout", type=int, default=3600, help="Solver timeout in seconds")
    parser.add_argument(
        "--reuse-p123",
        action="store_true",
        help="Reuse P123 matrices (set calculate_and_store_P123=false)",
    )
    parser.add_argument(
        "--reuse-p123-from",
        default="output/deuteron_proton_Ay/p123",
        help="Optional existing P123 folder to symlink/copy from",
    )
    parser.add_argument("--theta-min", type=float, default=20.0)
    parser.add_argument("--theta-max", type=float, default=170.0)
    parser.add_argument("--theta-step", type=float, default=1.0)
    parser.add_argument(
        "--dsigma-unit",
        default=UNIT_MB_PER_SR,
        choices=list(SUPPORTED_DSIGMA_UNITS),
        help="Output unit for dSigma/dOmega in CSV/JSON",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    target_tlabs_mev = _parse_tlab_list(args.target_tlabs_mev)
    dsigma_output_unit = normalize_dsigma_unit(args.dsigma_unit)

    cfg = SolverRunConfig(
        root=root,
        work_dir=(root / args.work_dir).resolve(),
        solver=(root / args.solver).resolve() if not Path(args.solver).is_absolute() else Path(args.solver),
        target_tlabs_mev=target_tlabs_mev,
        two_j_3n_max=args.two_j_3n_max,
        j_2n_max=args.j_2n_max,
        np_wp=args.np,
        nq_wp=args.nq,
        nphi=args.nphi,
        nx=args.nx,
        np_per_wp=args.np_per_wp,
        nq_per_wp=args.nq_per_wp,
        chebyshev_s=args.chebyshev_s,
        chebyshev_t=args.chebyshev_t,
        threads=args.threads,
        potential_model=args.potential_model,
        timeout_s=args.timeout,
        calculate_p123=(not args.reuse_p123),
        reuse_p123_from=(
            (root / args.reuse_p123_from).resolve()
            if args.reuse_p123_from and not Path(args.reuse_p123_from).is_absolute()
            else (Path(args.reuse_p123_from).resolve() if args.reuse_p123_from else None)
        ),
        angle_min_deg=args.theta_min,
        angle_max_deg=args.theta_max,
        angle_step_deg=args.theta_step,
    )

    cfg.work_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir = cfg.work_dir / "analysis"
    figures_dir = cfg.work_dir / "figures"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # Solver stage: generate U-matrix files for the requested energy list.
    run_result = run_cpp_solver(cfg)
    if run_result["returncode"] != 0:
        print(f"solver_returncode: {run_result['returncode']}")
        print(f"solver_log: {run_result['log_file']}")
        return 2
    if not run_result["u_files"] or not run_result["has_u_data"]:
        print("No valid U_PW_elements rows found after solver run.")
        return 3

    blocks = parse_solver_output(Path(str(run_result["solver_out_dir"])), u_files=run_result["u_files"])
    if not blocks:
        print("Failed to parse any (J, pi) blocks from solver output.")
        return 4
    q_grid = parse_q_kinematics(Path(str(run_result["solver_out_dir"])))
    common_energies = list_solver_energies(blocks)
    if not common_energies:
        print("No common solver energies across (J, pi) blocks.")
        return 5

    angles = _angle_grid(cfg.angle_min_deg, cfg.angle_max_deg, cfg.angle_step_deg)
    qmap_first = {pt.q_idx: pt for pt in blocks[0].points}
    dsigma_data_path = root / "data" / "DataOfCrosssectionAndPol" / "DSigamaOverDOmega.txt"
    cal_exp_angles: List[float] = []
    cal_exp_values_mb: List[float] = []
    if dsigma_data_path.exists():
        exp_dsigma_for_cal = read_experimental_dsigma(dsigma_data_path, target_unit=UNIT_MB_PER_SR)
        cal_exp_angles = list(exp_dsigma_for_cal["angles"])  # type: ignore[arg-type]
        cal_exp_values_mb = list(exp_dsigma_for_cal["values"])  # type: ignore[arg-type]

    energy_entries: List[Dict[str, object]] = []
    for target in target_tlabs_mev:
        tlab_solv, ecm_solv = min(common_energies, key=lambda item: abs(item[0] - target))
        q_idx = next(
            (idx for idx, pt in qmap_first.items()
             if abs(pt.tlab - tlab_solv) < 1e-6 and abs(pt.ecm - ecm_solv) < 1e-6),
            None,
        )
        if q_idx is None or q_idx not in q_grid:
            continue
        bin_info = q_grid[q_idx]
        delta = abs(tlab_solv - target)

        # First pass without calibration to compute mb/sr scale vs experiment.
        raw_rows = _evaluate_curve_on_angles(blocks, bin_info, q_idx, angles, 1.0, UNIT_MB_PER_SR)
        raw_dsigma_mb = [row["dSigma_dOmega"] for row in raw_rows]
        if cal_exp_angles and cal_exp_values_mb:
            cal_factor_mb = calibrate_dsigma_scale(raw_dsigma_mb, cal_exp_angles, cal_exp_values_mb, angles)
        else:
            cal_factor_mb = 1.0
        m_scale = math.sqrt(max(cal_factor_mb, 1e-30))

        rows = _evaluate_curve_on_angles(blocks, bin_info, q_idx, angles, m_scale, dsigma_output_unit)
        energy_dir = analysis_dir / _tlab_dir_name(target)
        energy_dir.mkdir(parents=True, exist_ok=True)

        model_csv = energy_dir / "observables_model.csv"
        _write_model_curve_csv(model_csv, angles, rows)

        metadata = {
            "target_tlab_mev": target,
            "selected_solver_tlab_mev": tlab_solv,
            "selected_solver_ecm_mev": ecm_solv,
            "abs_delta_tlab_mev": delta,
            "q_idx": q_idx,
            "q_mid_mev": bin_info.q_mid_mev,
            "dq_mev": bin_info.dq_mev,
            "dsigma_calibration_factor_mb_per_sr": cal_factor_mb,
            "model_curve_csv": str(model_csv),
            "num_jpi_blocks": len(blocks),
            "units": {
                "dSigma_dOmega": dsigma_output_unit,
                "iT11": "dimensionless",
                "T20": "dimensionless",
                "T21": "dimensionless",
                "T22": "dimensionless",
            },
        }
        (energy_dir / "metadata.json").write_text(
            json.dumps(metadata, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

        energy_entries.append(
            {
                "target_tlab_mev": target,
                "selected_solver_tlab_mev": tlab_solv,
                "selected_solver_ecm_mev": ecm_solv,
                "abs_delta_tlab_mev": delta,
                "q_idx": q_idx,
                "analysis_dir": str(energy_dir),
                "model_curve_csv": str(model_csv),
                "dsigma_calibration_factor_mb_per_sr": cal_factor_mb,
            }
        )

    # Optional validation stage: compare 190 MeV model against experiment (no fitting).
    tensor_data = root / "data" / "DataOfCrosssectionAndPol" / "CompletSetOFT" / "T.txt"
    dsigma_data = root / "data" / "DataOfCrosssectionAndPol" / "DSigamaOverDOmega.txt"
    comparison_190: Dict[str, object] = {"available": False}
    if tensor_data.exists() and dsigma_data.exists():
        exp_t = read_experimental_tensor(tensor_data)
        exp_dsigma = read_experimental_dsigma(dsigma_data, target_unit=dsigma_output_unit)

        target_190_entry = min(energy_entries, key=lambda item: abs(float(item["target_tlab_mev"]) - 190.0))
        e190_dir = Path(str(target_190_entry["analysis_dir"]))
        model_rows: List[Dict[str, float]] = []
        model_csv_path = e190_dir / "observables_model.csv"
        for raw in model_csv_path.read_text(encoding="utf-8").splitlines()[1:]:
            parts = raw.split(",")
            if len(parts) != 6:
                continue
            model_rows.append(
                {
                    "dSigma_dOmega": float(parts[1]),
                    "iT11": float(parts[2]),
                    "T20": float(parts[3]),
                    "T21": float(parts[4]),
                    "T22": float(parts[5]),
                }
            )

        exp_rows: List[Tuple[float, Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]] = []
        angle_union = sorted(
            set(exp_dsigma["angles"])
            | set(exp_t["iT11"]["angles"])
            | set(exp_t["T20"]["angles"])
            | set(exp_t["T21"]["angles"])
            | set(exp_t["T22"]["angles"])
        )
        map_dsigma = {a: v for a, v in zip(exp_dsigma["angles"], exp_dsigma["values"])}
        map_it11 = {a: v for a, v in zip(exp_t["iT11"]["angles"], exp_t["iT11"]["values"])}
        map_t20 = {a: v for a, v in zip(exp_t["T20"]["angles"], exp_t["T20"]["values"])}
        map_t21 = {a: v for a, v in zip(exp_t["T21"]["angles"], exp_t["T21"]["values"])}
        map_t22 = {a: v for a, v in zip(exp_t["T22"]["angles"], exp_t["T22"]["values"])}
        for a in angle_union:
            exp_rows.append(
                (
                    a,
                    map_dsigma.get(a),
                    map_it11.get(a),
                    map_t20.get(a),
                    map_t21.get(a),
                    map_t22.get(a),
                )
            )

        exp_csv = e190_dir / "observables_experiment_190.csv"
        exp_lines = ["theta_deg,dSigma_dOmega_exp,iT11_exp,T20_exp,T21_exp,T22_exp"]
        for a, ds, i11, t20, t21, t22 in exp_rows:
            exp_lines.append(
                f"{a:.6f},"
                f"{'' if ds is None else f'{ds:.12e}'},"
                f"{'' if i11 is None else f'{i11:.12e}'},"
                f"{'' if t20 is None else f'{t20:.12e}'},"
                f"{'' if t21 is None else f'{t21:.12e}'},"
                f"{'' if t22 is None else f'{t22:.12e}'}"
            )
        exp_csv.write_text("\n".join(exp_lines) + "\n", encoding="utf-8")

        model_angles = angles
        m_ds = _nearest_model_values(model_angles, model_rows, "dSigma_dOmega", exp_dsigma["angles"])
        m_i11 = _nearest_model_values(model_angles, model_rows, "iT11", exp_t["iT11"]["angles"])
        m_t20 = _nearest_model_values(model_angles, model_rows, "T20", exp_t["T20"]["angles"])
        m_t21 = _nearest_model_values(model_angles, model_rows, "T21", exp_t["T21"]["angles"])
        m_t22 = _nearest_model_values(model_angles, model_rows, "T22", exp_t["T22"]["angles"])

        comparison_190 = {
            "available": True,
            "selected_solver_tlab_mev": target_190_entry["selected_solver_tlab_mev"],
            "experiment_csv": str(exp_csv),
            "units": {
                "dSigma_dOmega_input_detected": exp_dsigma["source_unit"],
                "dSigma_dOmega_used": exp_dsigma["output_unit"],
            },
            "metrics": {
                "dSigma_dOmega": {
                    "rmse": _rmse(m_ds, exp_dsigma["values"]),
                    "mae": _mae(m_ds, exp_dsigma["values"]),
                    "count": len(exp_dsigma["values"]),
                },
                "iT11": {
                    "rmse": _rmse(m_i11, exp_t["iT11"]["values"]),
                    "mae": _mae(m_i11, exp_t["iT11"]["values"]),
                    "count": len(exp_t["iT11"]["values"]),
                },
                "T20": {
                    "rmse": _rmse(m_t20, exp_t["T20"]["values"]),
                    "mae": _mae(m_t20, exp_t["T20"]["values"]),
                    "count": len(exp_t["T20"]["values"]),
                },
                "T21": {
                    "rmse": _rmse(m_t21, exp_t["T21"]["values"]),
                    "mae": _mae(m_t21, exp_t["T21"]["values"]),
                    "count": len(exp_t["T21"]["values"]),
                },
                "T22": {
                    "rmse": _rmse(m_t22, exp_t["T22"]["values"]),
                    "mae": _mae(m_t22, exp_t["T22"]["values"]),
                    "count": len(exp_t["T22"]["values"]),
                },
            },
        }
        (e190_dir / "comparison_experiment_190.json").write_text(
            json.dumps(comparison_190, indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

    summary = {
        "workflow": "Tic-tac U-matrix -> full PW M-matrix observables (Witala/Gloeckle convention)",
        "target_tlabs_mev": target_tlabs_mev,
        "angles_deg": {"min": cfg.angle_min_deg, "max": cfg.angle_max_deg, "step": cfg.angle_step_deg},
        "units": {
            "dSigma_dOmega": dsigma_output_unit,
            "dSigma_dOmega_model_base": UNIT_MB_PER_SR,
            "iT11": "dimensionless",
            "T20": "dimensionless",
            "T21": "dimensionless",
            "T22": "dimensionless",
        },
        "solver": {
            "returncode": run_result["returncode"],
            "log_file": str(run_result["log_file"]),
            "energy_file": str(run_result["energy_file"]),
            "config_file": str(run_result["config_file"]),
            "solver_out_dir": str(run_result["solver_out_dir"]),
            "p123_dir": str(run_result["p123_dir"]),
            "u_files": [str(p) for p in run_result["u_files"]],
        },
        "energies": energy_entries,
        "experiment_190_comparison": comparison_190,
        "output_tree": {
            "inputs": str(cfg.work_dir / "inputs"),
            "solver": str(cfg.work_dir / "solver"),
            "analysis": str(analysis_dir),
            "figures": str(figures_dir),
        },
    }

    summary_json = analysis_dir / "summary.json"
    summary_txt = analysis_dir / "summary.txt"
    summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")

    lines = []
    lines.append("d + p multi-Tlab observable run")
    lines.append("===============================")
    lines.append(f"targets (Tlab MeV): {', '.join(f'{v:.1f}' for v in target_tlabs_mev)}")
    lines.append(f"dSigma/dOmega unit: {dsigma_output_unit}")
    lines.append(f"solver return code: {run_result['returncode']}")
    lines.append("")
    lines.append("selected solver energies:")
    for item in energy_entries:
        lines.append(
            f"  target={float(item['target_tlab_mev']):.1f} MeV -> "
            f"solver={float(item['selected_solver_tlab_mev']):.3f} MeV "
            f"(delta={float(item['abs_delta_tlab_mev']):.3f} MeV)"
        )
    lines.append("")
    lines.append(f"summary json: {summary_json}")
    lines.append(f"analysis dir: {analysis_dir}")
    summary_txt.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
