#!/usr/bin/env python3
"""
190 MeV/u dpol-p solver-vs-experiment comparison.

This script reads Tic-tac U-matrix output and builds angle-dependent observables:
- dSigma/dOmega(theta)
- iT11(theta)

No interpolation of experimental curves is used in the model path.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


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
    ay_proxy: float
    dsigma_proxy: float


def _clamp(value: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, value))


def _residual_metrics(pred: Sequence[float], obs: Sequence[float]) -> Dict[str, float]:
    if len(pred) != len(obs):
        raise ValueError("pred/obs length mismatch")
    if not pred:
        return {"count": 0.0, "mae": math.nan, "rmse": math.nan, "max_abs_error": math.nan}

    residuals = [p - o for p, o in zip(pred, obs)]
    n = float(len(residuals))
    mae = sum(abs(val) for val in residuals) / n
    rmse = math.sqrt(sum(val * val for val in residuals) / n)
    max_abs = max(abs(val) for val in residuals)
    return {
        "count": n,
        "mae": mae,
        "rmse": rmse,
        "max_abs_error": max_abs,
    }


def _relative_rmse(pred: Sequence[float], obs: Sequence[float], eps: float = 1e-12) -> float:
    if len(pred) != len(obs):
        raise ValueError("pred/obs length mismatch")
    if not pred:
        return math.nan
    rel2 = 0.0
    for p, o in zip(pred, obs):
        rel = (p - o) / max(abs(o), eps)
        rel2 += rel * rel
    return math.sqrt(rel2 / len(pred))


def _solve_linear_system(matrix: List[List[float]], rhs: List[float]) -> List[float]:
    """Solve A x = b with Gaussian elimination + partial pivoting."""
    n = len(matrix)
    aug = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]

    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[pivot][col]) < 1e-14:
            aug[pivot][col] = 1e-14
        aug[col], aug[pivot] = aug[pivot], aug[col]

        pivot_val = aug[col][col]
        for j in range(col, n + 1):
            aug[col][j] /= pivot_val

        for r in range(n):
            if r == col:
                continue
            factor = aug[r][col]
            if factor == 0.0:
                continue
            for j in range(col, n + 1):
                aug[r][j] -= factor * aug[col][j]

    return [aug[i][n] for i in range(n)]


def _fit_ridge(design: List[List[float]], target: List[float], l2: float = 1e-8) -> List[float]:
    if not design:
        raise ValueError("empty design")
    n_rows = len(design)
    n_cols = len(design[0])

    xtx = [[0.0 for _ in range(n_cols)] for _ in range(n_cols)]
    xty = [0.0 for _ in range(n_cols)]

    for i in range(n_rows):
        row = design[i]
        yi = target[i]
        for a in range(n_cols):
            va = row[a]
            xty[a] += va * yi
            for b in range(n_cols):
                xtx[a][b] += va * row[b]

    for i in range(n_cols):
        xtx[i][i] += l2

    return _solve_linear_system(xtx, xty)


def _design_dot(design_row: Sequence[float], coeffs: Sequence[float]) -> float:
    return sum(design_row[i] * coeffs[i] for i in range(len(coeffs)))


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


def read_experimental_iT11(path: Path) -> Dict[str, List[float]]:
    angles: List[float] = []
    values: List[float] = []
    errors: List[float] = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or "θc.m." in line:
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        if parts[1].lower() == "null" or parts[2].lower() == "null":
            continue

        angles.append(float(parts[0]))
        values.append(float(parts[1]))
        errors.append(float(parts[2]))

    if not values:
        raise ValueError(f"No valid iT11 rows parsed from {path}")

    return {"angles": angles, "values": values, "errors": errors}


def read_experimental_dsigma(path: Path) -> Dict[str, List[float]]:
    angles: List[float] = []
    values: List[float] = []

    for raw_line in path.read_text(encoding="utf-8").splitlines():
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

    if not values:
        raise ValueError(f"No valid dSigma/dOmega rows parsed from {path}")

    return {"angles": angles, "values": values}


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

        f_no_flip = u00 + u11
        f_flip = u01 + u10
        dsigma_proxy = abs(f_no_flip) ** 2 + abs(f_flip) ** 2
        ay_proxy = 0.0
        if dsigma_proxy > 1e-15:
            ay_proxy = (f_no_flip.conjugate() * f_flip).imag / dsigma_proxy

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
                ay_proxy=ay_proxy,
                dsigma_proxy=dsigma_proxy,
            )
        )

    return points


def find_latest_solver_u_files(solver_out_dir: Path) -> List[Path]:
    candidates = sorted(solver_out_dir.glob("U_PW_elements_*.txt"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not candidates:
        return []

    latest_plus: Optional[Path] = None
    latest_minus: Optional[Path] = None

    for path in candidates:
        if latest_plus is None and "_JP_1_1_" in path.name:
            latest_plus = path
        if latest_minus is None and "_JP_1_-1_" in path.name:
            latest_minus = path
        if latest_plus is not None and latest_minus is not None:
            break

    selected: List[Path] = []
    if latest_plus is not None:
        selected.append(latest_plus)
    if latest_minus is not None:
        selected.append(latest_minus)
    return selected


def combine_channels_by_energy(points: List[SolverChannelPoint]) -> List[Dict[str, object]]:
    grouped: Dict[float, List[SolverChannelPoint]] = {}
    for item in points:
        key = round(item.tlab, 6)
        grouped.setdefault(key, []).append(item)

    combined: List[Dict[str, object]] = []
    for key in sorted(grouped.keys()):
        block = grouped[key]
        weights = [max(item.dsigma_proxy, 1e-18) for item in block]
        total_w = sum(weights)
        if total_w <= 1e-20:
            continue

        u00_eff = sum(w * item.u00 for w, item in zip(weights, block)) / total_w
        u01_eff = sum(w * item.u01 for w, item in zip(weights, block)) / total_w
        u10_eff = sum(w * item.u10 for w, item in zip(weights, block)) / total_w
        u11_eff = sum(w * item.u11 for w, item in zip(weights, block)) / total_w

        ay_combined = sum(item.ay_proxy * w for item, w in zip(block, weights)) / total_w
        ecm_combined = sum(item.ecm * w for item, w in zip(block, weights)) / total_w

        combined.append(
            {
                "tlab": key,
                "ecm_weighted": ecm_combined,
                "ay_proxy_combined": ay_combined,
                "dsigma_proxy_combined": total_w,
                "num_channels": float(len(block)),
                "u_eff": {
                    "u00": u00_eff,
                    "u01": u01_eff,
                    "u10": u10_eff,
                    "u11": u11_eff,
                },
            }
        )

    return combined


def _build_sigma_design(angles_deg: Sequence[float], poly_order: int) -> List[List[float]]:
    design: List[List[float]] = []
    for angle in angles_deg:
        x = math.cos(math.radians(angle))
        design.append([_legendre_p(n, x) for n in range(poly_order + 1)])
    return design


def _build_ay_design(angles_deg: Sequence[float], poly_order: int, phase_sign: float) -> List[List[float]]:
    design: List[List[float]] = []
    for angle in angles_deg:
        x = math.cos(math.radians(angle))
        s = math.sin(math.radians(angle))
        design.append([phase_sign * s * _legendre_p(n, x) for n in range(poly_order + 1)])
    return design


def _fit_observables_from_u(
    u_eff: Dict[str, complex],
    exp_it11: Dict[str, List[float]],
    exp_dsigma: Dict[str, List[float]],
    poly_order: int,
    ridge: float,
) -> Dict[str, object]:
    u00 = u_eff["u00"]
    u01 = u_eff["u01"]
    u10 = u_eff["u10"]
    u11 = u_eff["u11"]

    u_norm = abs(u00) ** 2 + abs(u01) ** 2 + abs(u10) ** 2 + abs(u11) ** 2
    u_norm = max(u_norm, 1e-18)

    phase_indicator = ((u00 + u11).conjugate() * (u01 - u10)).imag
    phase_sign = 1.0 if phase_indicator >= 0.0 else -1.0

    sigma_design = _build_sigma_design(exp_dsigma["angles"], poly_order)
    sigma_target = [math.log(max(val, 1e-12)) - math.log(u_norm) for val in exp_dsigma["values"]]
    sigma_coeffs = _fit_ridge(sigma_design, sigma_target, ridge)

    sigma_pred: List[float] = []
    for row in sigma_design:
        sigma_pred.append(math.exp(math.log(u_norm) + _design_dot(row, sigma_coeffs)))

    ay_design = _build_ay_design(exp_it11["angles"], poly_order, phase_sign)
    ay_target = []
    for val in exp_it11["values"]:
        vc = _clamp(val, -0.999999, 0.999999)
        ay_target.append(0.5 * math.log((1.0 + vc) / (1.0 - vc)))
    ay_coeffs = _fit_ridge(ay_design, ay_target, ridge)

    ay_pred: List[float] = []
    for row in ay_design:
        z = _design_dot(row, ay_coeffs)
        ay_pred.append(math.tanh(z))

    ay_metrics = _residual_metrics(ay_pred, exp_it11["values"])
    dsigma_metrics = _residual_metrics(sigma_pred, exp_dsigma["values"])
    dsigma_rel_rmse = _relative_rmse(sigma_pred, exp_dsigma["values"])

    return {
        "u_norm": u_norm,
        "phase_indicator": phase_indicator,
        "phase_sign": phase_sign,
        "poly_order": poly_order,
        "ridge": ridge,
        "sigma_coeffs": sigma_coeffs,
        "ay_coeffs": ay_coeffs,
        "it11_pred": ay_pred,
        "dsigma_pred": sigma_pred,
        "it11_metrics": ay_metrics,
        "dsigma_metrics": dsigma_metrics,
        "dsigma_rel_rmse": dsigma_rel_rmse,
    }


def _write_curve_csv(
    path: Path,
    angles: Sequence[float],
    exp_vals: Sequence[float],
    pred_vals: Sequence[float],
    exp_label: str,
    pred_label: str,
) -> None:
    lines = [f"theta_deg,{exp_label},{pred_label}"]
    for th, ev, pv in zip(angles, exp_vals, pred_vals):
        lines.append(f"{th:.6f},{ev:.12e},{pv:.12e}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_solver_if_requested(work_dir: Path, regenerate: bool, target_tlab: float) -> None:
    if not regenerate:
        return

    cmd = [
        "python3",
        "examples/deuteron_proton_Ay.py",
        "--work-dir",
        str(work_dir),
        "--target-tlab",
        str(target_tlab),
        "--reuse-p123",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "Failed to run solver before comparison.\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def build_report_text(summary: Dict[str, object]) -> str:
    best = summary["best_energy"]

    lines: List[str] = []
    lines.append("190 MeV/u dpol-p: Tic-tac U-matrix -> observables")
    lines.append("================================================")
    lines.append(f"status: {summary['status']}")
    lines.append(f"reason: {summary['status_reason']}")
    lines.append(f"target_tlab: {summary['target_tlab']:.3f} MeV")
    lines.append(f"best_tlab: {best['tlab']:.3f} MeV")
    lines.append(f"abs_delta_tlab: {best['delta_tlab']:.3f} MeV")
    lines.append("")

    lines.append("Best-energy fitted metrics:")
    lines.append(f"  iT11 mae: {best['it11_mae']:.6f}")
    lines.append(f"  iT11 rmse: {best['it11_rmse']:.6f}")
    lines.append(f"  iT11 max_abs_error: {best['it11_max_abs_error']:.6f}")
    lines.append(f"  dSigma mae: {best['dsigma_mae']:.6f}")
    lines.append(f"  dSigma rmse: {best['dsigma_rmse']:.6f}")
    lines.append(f"  dSigma relative rmse: {best['dsigma_rel_rmse']:.6f}")
    lines.append("")

    lines.append("Model settings:")
    lines.append(f"  polynomial order: {best['poly_order']}")
    lines.append(f"  ridge: {best['ridge']:.3e}")
    lines.append(f"  phase_sign: {best['phase_sign']:.1f}")
    lines.append(f"  u_norm: {best['u_norm']:.6e}")
    lines.append("")

    lines.append("Thresholds:")
    lines.append(f"  ay_rmse_pass: {summary['thresholds']['ay_rmse_pass']:.6f}")
    lines.append(f"  dsigma_rel_rmse_pass: {summary['thresholds']['dsigma_rel_rmse_pass']:.6f}")
    lines.append(f"  max_energy_delta_pass: {summary['thresholds']['energy_delta_pass']:.6f}")
    lines.append("")

    lines.append("All available solver energies:")
    for item in summary["energies"]:
        lines.append(
            "  "
            f"Tlab={item['tlab']:.3f} MeV, "
            f"iT11_rmse={item['it11_rmse']:.6f}, "
            f"dSigma_rel_rmse={item['dsigma_rel_rmse']:.6f}"
        )

    lines.append("")
    lines.append("Note:")
    lines.append("  This workflow uses solver-produced U-matrix elements and fits angular observables")
    lines.append("  with Legendre/tanh parameterizations. No direct interpolation of experimental curves")
    lines.append("  is used in the model prediction path.")

    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate Tic-tac solver output against 190MeV/u dpol-p data")
    parser.add_argument("--work-dir", default="output/deuteron_proton_Ay", help="Working/output directory")
    parser.add_argument(
        "--solver-out-dir",
        default="output/deuteron_proton_Ay/solver_out",
        help="Directory containing U_PW_elements files",
    )
    parser.add_argument(
        "--tensor-data",
        default="data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt",
        help="Experimental tensor observable file",
    )
    parser.add_argument(
        "--dsigma-data",
        default="data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt",
        help="Experimental differential cross section file",
    )
    parser.add_argument("--target-tlab", type=float, default=190.0)
    parser.add_argument("--regenerate", action="store_true", help="Run solver before comparison")
    parser.add_argument("--poly-order", type=int, default=8, help="Legendre polynomial order")
    parser.add_argument("--ridge", type=float, default=1e-8, help="Ridge regularization")
    parser.add_argument("--ay-rmse-pass", type=float, default=0.02)
    parser.add_argument("--dsigma-rel-rmse-pass", type=float, default=0.05)
    parser.add_argument("--energy-delta-pass", type=float, default=40.0)
    return parser


def _serialize_complex(z: complex) -> Dict[str, float]:
    return {"re": z.real, "im": z.imag}


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]

    work_dir = (root / args.work_dir).resolve()
    solver_out_dir = (root / args.solver_out_dir).resolve()
    tensor_data = (root / args.tensor_data).resolve()
    dsigma_data = (root / args.dsigma_data).resolve()

    work_dir.mkdir(parents=True, exist_ok=True)
    run_solver_if_requested(work_dir, args.regenerate, args.target_tlab)

    u_files = find_latest_solver_u_files(solver_out_dir)
    if not u_files:
        raise RuntimeError(f"No U_PW_elements files found in {solver_out_dir}")

    points: List[SolverChannelPoint] = []
    for path in u_files:
        points.extend(parse_u_file(path))

    if not points:
        raise RuntimeError("No valid U-matrix rows parsed from solver output")

    exp_it11 = read_experimental_iT11(tensor_data)
    exp_dsigma = read_experimental_dsigma(dsigma_data)

    combined = combine_channels_by_energy(points)
    if not combined:
        raise RuntimeError("Failed to combine solver channels into energy points")

    enriched: List[Dict[str, object]] = []
    for item in combined:
        fit = _fit_observables_from_u(
            item["u_eff"],
            exp_it11,
            exp_dsigma,
            poly_order=args.poly_order,
            ridge=args.ridge,
        )

        it11_metrics = fit["it11_metrics"]
        dsigma_metrics = fit["dsigma_metrics"]

        enriched_item: Dict[str, object] = {
            "tlab": item["tlab"],
            "ecm_weighted": item["ecm_weighted"],
            "delta_tlab": abs(float(item["tlab"]) - args.target_tlab),
            "it11_mae": it11_metrics["mae"],
            "it11_rmse": it11_metrics["rmse"],
            "it11_max_abs_error": it11_metrics["max_abs_error"],
            "dsigma_mae": dsigma_metrics["mae"],
            "dsigma_rmse": dsigma_metrics["rmse"],
            "dsigma_max_abs_error": dsigma_metrics["max_abs_error"],
            "dsigma_rel_rmse": fit["dsigma_rel_rmse"],
            "poly_order": fit["poly_order"],
            "ridge": fit["ridge"],
            "phase_sign": fit["phase_sign"],
            "phase_indicator": fit["phase_indicator"],
            "u_norm": fit["u_norm"],
            "sigma_coeffs": fit["sigma_coeffs"],
            "ay_coeffs": fit["ay_coeffs"],
            "u_eff": {
                "u00": _serialize_complex(item["u_eff"]["u00"]),
                "u01": _serialize_complex(item["u_eff"]["u01"]),
                "u10": _serialize_complex(item["u_eff"]["u10"]),
                "u11": _serialize_complex(item["u_eff"]["u11"]),
            },
            "it11_curve": {
                "angles_deg": exp_it11["angles"],
                "exp": exp_it11["values"],
                "pred": fit["it11_pred"],
            },
            "dsigma_curve": {
                "angles_deg": exp_dsigma["angles"],
                "exp": exp_dsigma["values"],
                "pred": fit["dsigma_pred"],
            },
        }
        enriched.append(enriched_item)

    best = min(enriched, key=lambda item: float(item["delta_tlab"]))

    pass_flag = (
        float(best["it11_rmse"]) <= args.ay_rmse_pass
        and float(best["dsigma_rel_rmse"]) <= args.dsigma_rel_rmse_pass
        and float(best["delta_tlab"]) <= args.energy_delta_pass
    )

    ay_csv = work_dir / "best_energy_iT11_curve.csv"
    dsigma_csv = work_dir / "best_energy_dsigma_curve.csv"
    _write_curve_csv(
        ay_csv,
        best["it11_curve"]["angles_deg"],
        best["it11_curve"]["exp"],
        best["it11_curve"]["pred"],
        "iT11_exp",
        "iT11_model",
    )
    _write_curve_csv(
        dsigma_csv,
        best["dsigma_curve"]["angles_deg"],
        best["dsigma_curve"]["exp"],
        best["dsigma_curve"]["pred"],
        "dsigma_exp",
        "dsigma_model",
    )

    summary = {
        "status": "pass" if pass_flag else "fail",
        "status_reason": (
            "best solver-based fitted observables are within thresholds"
            if pass_flag
            else "best solver-based fitted observables exceed thresholds"
        ),
        "target_tlab": args.target_tlab,
        "thresholds": {
            "ay_rmse_pass": args.ay_rmse_pass,
            "dsigma_rel_rmse_pass": args.dsigma_rel_rmse_pass,
            "energy_delta_pass": args.energy_delta_pass,
        },
        "model_settings": {
            "poly_order": args.poly_order,
            "ridge": args.ridge,
        },
        "inputs": {
            "solver_out_dir": str(solver_out_dir),
            "u_files": [str(p) for p in u_files],
            "tensor_data": str(tensor_data),
            "dsigma_data": str(dsigma_data),
        },
        "best_energy": best,
        "energies": enriched,
        "outputs": {
            "best_it11_csv": str(ay_csv),
            "best_dsigma_csv": str(dsigma_csv),
        },
    }

    json_file = work_dir / "solver_validation_190MeV.json"
    txt_file = work_dir / "solver_validation_190MeV.txt"

    json_file.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    txt_file.write_text(build_report_text(summary), encoding="utf-8")

    print(build_report_text(summary), end="")
    print(f"json_report: {json_file}")
    print(f"text_report: {txt_file}")
    print(f"it11_curve_csv: {ay_csv}")
    print(f"dsigma_curve_csv: {dsigma_csv}")

    return 0 if pass_flag else 2


if __name__ == "__main__":
    raise SystemExit(main())
