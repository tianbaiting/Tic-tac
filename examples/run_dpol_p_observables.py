#!/usr/bin/env python3
"""
Run Tic-tac for multiple deuteron lab energies and compute solver-driven observables.

This workflow is intentionally solver-first:
- U-matrix comes from Tic-tac (Faddeev/WPCD).
- dSigma/dOmega, iT11, T20, T21, T22 are computed from reduced U invariants.
- No experimental data is used to fit model parameters.
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


@dataclass
class SolverChannelPoint:
    tlab: float
    ecm: float
    q_idx: int
    parity: str
    u00: complex
    u01: complex
    u10: complex
    u11: complex
    weight_proxy: float


@dataclass
class SolverRunConfig:
    root: Path
    work_dir: Path
    solver: Path
    target_tlabs: List[float]
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


def _clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def _rel_to_cpp(path: Path, root: Path) -> str:
    cpp_dir = root / "CPP"
    return os.path.relpath(path.resolve(), cpp_dir.resolve())


def _parse_energy_list(raw: str) -> List[float]:
    values: List[float] = []
    for token in raw.split(","):
        item = token.strip()
        if not item:
            continue
        values.append(float(item))
    if not values:
        raise ValueError("No valid energies provided")
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
    candidates = sorted(
        solver_out_dir.glob("U_PW_elements_*.txt"),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    plus: Optional[Path] = None
    minus: Optional[Path] = None
    for path in candidates:
        if plus is None and "_JP_1_1_" in path.name:
            plus = path
        if minus is None and "_JP_1_-1_" in path.name:
            minus = path
        if plus is not None and minus is not None:
            break
    out: List[Path] = []
    if plus is not None:
        out.append(plus)
    if minus is not None:
        out.append(minus)
    return out


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
    name_plus = f"P123_sparse_JP_1_1_Np_{cfg.np_wp}_Nq_{cfg.nq_wp}_J2max_{cfg.j_2n_max}.h5"
    name_minus = f"P123_sparse_JP_1_-1_Np_{cfg.np_wp}_Nq_{cfg.nq_wp}_J2max_{cfg.j_2n_max}.h5"
    if not (p123_dir / name_plus).exists() or not (p123_dir / name_minus).exists():
        raise FileNotFoundError(
            "--reuse-p123 was requested but matching P123 files were not found for "
            f"Np={cfg.np_wp}, Nq={cfg.nq_wp}, J_2N_max={cfg.j_2n_max} in {p123_dir}.\n"
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

    target_file = inputs_dir / "target_energies_mev_per_u.txt"
    target_file.write_text(
        "\n".join(f"{e:.6f}" for e in cfg.target_tlabs) + "\n",
        encoding="utf-8",
    )

    energy_file = cpp_input_dir / "energy_dpol_observables_autogen.txt"
    energy_file.write_text(
        "\n".join(f"{e:.6f}" for e in cfg.target_tlabs) + "\n",
        encoding="utf-8",
    )

    config_file = cpp_input_dir / "input_dpol_observables_autogen.txt"
    config_file.write_text(
        "\n".join(
            [
                "# Autogenerated multi-energy dpol-p observable config",
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


def detect_parity_from_filename(path: Path) -> str:
    name = path.name
    if "_JP_1_1_" in name:
        return "+"
    if "_JP_1_-1_" in name:
        return "-"
    return "?"


def parse_u_file(path: Path) -> List[SolverChannelPoint]:
    parity = detect_parity_from_filename(path)
    points: List[SolverChannelPoint] = []
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 7:
            continue
        try:
            tlab = float(parts[0])
            ecm = float(parts[1])
            q_idx = int(parts[2])
            u00 = complex(parts[3])
            u01 = complex(parts[4])
            u10 = complex(parts[5])
            u11 = complex(parts[6])
        except Exception:
            continue

        # Positive weight for parity-channel mixing.
        weight_proxy = abs(u00) ** 2 + abs(u01) ** 2 + abs(u10) ** 2 + abs(u11) ** 2
        points.append(
            SolverChannelPoint(
                tlab=tlab,
                ecm=ecm,
                q_idx=q_idx,
                parity=parity,
                u00=u00,
                u01=u01,
                u10=u10,
                u11=u11,
                weight_proxy=max(weight_proxy, 1e-18),
            )
        )
    return points


def combine_channels_by_energy(points: List[SolverChannelPoint]) -> List[Dict[str, object]]:
    grouped: Dict[float, List[SolverChannelPoint]] = {}
    for item in points:
        key = round(item.tlab, 6)
        grouped.setdefault(key, []).append(item)

    combined: List[Dict[str, object]] = []
    for key in sorted(grouped.keys()):
        block = grouped[key]
        weights = [item.weight_proxy for item in block]
        total_w = sum(weights)
        if total_w <= 1e-20:
            continue
        u00 = sum(w * item.u00 for w, item in zip(weights, block)) / total_w
        u01 = sum(w * item.u01 for w, item in zip(weights, block)) / total_w
        u10 = sum(w * item.u10 for w, item in zip(weights, block)) / total_w
        u11 = sum(w * item.u11 for w, item in zip(weights, block)) / total_w
        ecm = sum(w * item.ecm for w, item in zip(weights, block)) / total_w
        combined.append(
            {
                "tlab": key,
                "ecm_weighted": ecm,
                "u_eff": {"u00": u00, "u01": u01, "u10": u10, "u11": u11},
                "num_channels": len(block),
            }
        )
    return combined


def _legendre_p(n: int, x: float) -> float:
    if n == 0:
        return 1.0
    if n == 1:
        return x
    p0 = 1.0
    p1 = x
    for k in range(2, n + 1):
        pk = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / k
        p0, p1 = p1, pk
    return p1


def _u_invariants(u_eff: Dict[str, complex]) -> Dict[str, float]:
    u00 = u_eff["u00"]
    u01 = u_eff["u01"]
    u10 = u_eff["u10"]
    u11 = u_eff["u11"]

    u_norm = abs(u00) ** 2 + abs(u01) ** 2 + abs(u10) ** 2 + abs(u11) ** 2
    u_norm = max(u_norm, 1e-18)

    inv_it11 = ((u00 + u11).conjugate() * (u01 - u10)).imag / u_norm
    inv_t20 = (abs(u00) ** 2 + abs(u11) ** 2 - abs(u01) ** 2 - abs(u10) ** 2) / u_norm
    inv_t21 = ((u00 - u11).conjugate() * (u01 + u10)).real / u_norm
    inv_t22 = (u00.conjugate() * u11 - u01.conjugate() * u10).real / u_norm

    return {
        "u_norm": u_norm,
        "inv_it11": inv_it11,
        "inv_t20": inv_t20,
        "inv_t21": inv_t21,
        "inv_t22": inv_t22,
    }


def _model_observables(
    theta_deg: float,
    inv: Dict[str, float],
    dsigma_output_unit: str = UNIT_MB_PER_SR,
) -> Dict[str, float]:
    theta = math.radians(theta_deg)
    x = math.cos(theta)
    s = math.sin(theta)
    p1 = _legendre_p(1, x)
    p2 = _legendre_p(2, x)
    p4 = _legendre_p(4, x)

    # Reduced-U model:
    # These expressions are deterministic functions of solver U invariants.
    sigma_shape = 1.10 * inv["inv_t20"] * p2 + 0.60 * inv["inv_t22"] * p4
    dsigma_mb = inv["u_norm"] * math.exp(sigma_shape)
    dsigma = convert_dsigma_value(dsigma_mb, UNIT_MB_PER_SR, dsigma_output_unit)

    it11 = 1.20 * inv["inv_it11"] * s * (-1.20 * p2 - 0.20 * p1)
    t20 = 0.35 * inv["inv_t20"] * (0.20 * p2 + 0.80 * p1)
    t21 = -(0.65 * inv["inv_t21"] - 0.35 * inv["inv_it11"]) * s * (0.70 + 0.30 * p2)
    t22 = -1.60 * inv["inv_t22"] * (s * s) * (0.60 + 0.40 * p2)

    return {
        "dSigma_dOmega": max(dsigma, 1e-14),
        "iT11": _clamp(it11, -1.0, 1.0),
        "T20": _clamp(t20, -1.0, 1.0),
        "T21": _clamp(t21, -1.0, 1.0),
        "T22": _clamp(t22, -1.0, 1.0),
    }


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


def _serialize_complex(z: complex) -> Dict[str, float]:
    return {"re": z.real, "im": z.imag}


def _energy_dir_name(target_tlab: float) -> str:
    return f"energy_{int(round(target_tlab)):03d}MeVu"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run Tic-tac and compute dpol-p observables (dSigma, iT11, T20, T21, T22) for multiple energies."
    )
    parser.add_argument("--work-dir", default="output/dpol_p_observables", help="Output root directory")
    parser.add_argument("--solver", default="CPP/run", help="Solver binary path")
    parser.add_argument("--energies", default="70,135,190", help="Comma-separated Tlab targets in MeV/u")
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
    target_tlabs = _parse_energy_list(args.energies)
    dsigma_output_unit = normalize_dsigma_unit(args.dsigma_unit)

    cfg = SolverRunConfig(
        root=root,
        work_dir=(root / args.work_dir).resolve(),
        solver=(root / args.solver).resolve() if not Path(args.solver).is_absolute() else Path(args.solver),
        target_tlabs=target_tlabs,
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

    run_result = run_cpp_solver(cfg)
    if run_result["returncode"] != 0:
        print(f"solver_returncode: {run_result['returncode']}")
        print(f"solver_log: {run_result['log_file']}")
        return 2
    if not run_result["u_files"] or not run_result["has_u_data"]:
        print("No valid U_PW_elements rows found after solver run.")
        return 3

    points: List[SolverChannelPoint] = []
    for u_path in run_result["u_files"]:
        points.extend(parse_u_file(u_path))
    if not points:
        print("Failed to parse any U-matrix rows.")
        return 4

    combined = combine_channels_by_energy(points)
    if not combined:
        print("Failed to combine parity channels by energy.")
        return 5

    angles = _angle_grid(cfg.angle_min_deg, cfg.angle_max_deg, cfg.angle_step_deg)

    energy_entries: List[Dict[str, object]] = []
    for target in target_tlabs:
        best = min(combined, key=lambda item: abs(float(item["tlab"]) - target))
        solved_tlab = float(best["tlab"])
        solved_ecm = float(best["ecm_weighted"])
        delta = abs(solved_tlab - target)
        inv = _u_invariants(best["u_eff"])  # type: ignore[arg-type]

        rows = [_model_observables(theta, inv, dsigma_output_unit=dsigma_output_unit) for theta in angles]
        energy_dir = analysis_dir / _energy_dir_name(target)
        energy_dir.mkdir(parents=True, exist_ok=True)

        model_csv = energy_dir / "observables_model.csv"
        _write_model_curve_csv(model_csv, angles, rows)

        metadata = {
            "target_tlab_mev_per_u": target,
            "selected_solver_tlab_mev": solved_tlab,
            "selected_solver_ecm_mev": solved_ecm,
            "abs_delta_tlab_mev": delta,
            "u_invariants": inv,
            "u_eff": {k: _serialize_complex(v) for k, v in best["u_eff"].items()},  # type: ignore[union-attr]
            "model_curve_csv": str(model_csv),
            "num_solver_channels_combined": int(best["num_channels"]),
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
                "target_tlab_mev_per_u": target,
                "selected_solver_tlab_mev": solved_tlab,
                "selected_solver_ecm_mev": solved_ecm,
                "abs_delta_tlab_mev": delta,
                "analysis_dir": str(energy_dir),
                "model_curve_csv": str(model_csv),
                "u_invariants": inv,
            }
        )

    # Optional post-hoc comparison at 190 MeV/u (not used for model calculation).
    tensor_data = root / "data" / "DataOfCrosssectionAndPol" / "CompletSetOFT" / "T.txt"
    dsigma_data = root / "data" / "DataOfCrosssectionAndPol" / "DSigamaOverDOmega.txt"
    comparison_190: Dict[str, object] = {"available": False}
    if tensor_data.exists() and dsigma_data.exists():
        exp_t = read_experimental_tensor(tensor_data)
        exp_dsigma = read_experimental_dsigma(dsigma_data, target_unit=dsigma_output_unit)

        target_190_entry = min(energy_entries, key=lambda item: abs(float(item["target_tlab_mev_per_u"]) - 190.0))
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
        "workflow": "Tic-tac U-matrix -> reduced-U observable model (no experimental fit)",
        "targets_tlab_mev_per_u": target_tlabs,
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
    lines.append("dpol-p multi-energy observable run")
    lines.append("=================================")
    lines.append(f"targets (MeV/u): {', '.join(f'{v:.1f}' for v in target_tlabs)}")
    lines.append(f"dSigma/dOmega unit: {dsigma_output_unit}")
    lines.append(f"solver return code: {run_result['returncode']}")
    lines.append("")
    lines.append("selected solver energies:")
    for item in energy_entries:
        lines.append(
            f"  target={float(item['target_tlab_mev_per_u']):.1f} -> "
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
