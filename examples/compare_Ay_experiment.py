#!/usr/bin/env python3
"""
Purpose:
  Compare solver-driven observables against the maintained Tlab benchmark in MeV.

What this script does:
  1) Read solver U-matrix files from solver_out.
  2) Combine parity channels at identical solver energy.
  3) Build reduced-U invariants and predict:
     - dSigma/dOmega(theta)
     - iT11(theta)
  4) Evaluate residuals against experiment and write reports/curves.

Data flow (left -> right):
  U_PW_elements_*.txt + experiment txt files
    -> parse/merge channels
    -> reduced-U map prediction
    -> metrics
    -> CSV/SVG/JSON/TXT in work_dir

Calls / dependencies:
  - Uses `observable_units.py` for unit parsing/conversion.
  - Uses `solver_u_file_utils.py` for selecting latest U-file family.
  - Optional upstream call: runs `examples/deuteron_proton_Ay.py` when `--regenerate` is set.

Usage:
  python3 examples/compare_Ay_experiment.py \
    --work-dir output/deuteron_proton_Ay \
    --solver-out-dir output/deuteron_proton_Ay/solver_out \
    --target-tlab-mev 190 \
    --dsigma-unit fm2/sr
"""

from __future__ import annotations

import argparse
import json
import math
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
from solver_u_file_utils import detect_parity_symbol, select_latest_u_file_family
from tlab_utils import format_tlab_label


@dataclass
class SolverChannelPoint:
    tlab: float
    ecm: float
    q_idx: int
    two_j: int
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

    if not values:
        raise ValueError(f"No valid dSigma/dOmega rows parsed from {path}")

    return {
        "angles": angles,
        "values": convert_dsigma_series(values, source_unit, output_unit),
        "source_unit": source_unit,
        "output_unit": output_unit,
    }


def detect_parity_from_filename(path: Path) -> str:
    return detect_parity_symbol(path)


def detect_two_j_from_filename(path: Path) -> int:
    name = path.name
    marker = "_JP_"
    idx = name.find(marker)
    if idx < 0:
        return 1
    payload = name[idx + len(marker) :]
    fields = payload.split("_")
    if not fields:
        return 1
    try:
        return int(fields[0])
    except ValueError:
        return 1


def parse_u_file(path: Path) -> List[SolverChannelPoint]:
    parity = detect_parity_from_filename(path)
    two_j = detect_two_j_from_filename(path)
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
                two_j=two_j,
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
    return select_latest_u_file_family(solver_out_dir)


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


def _predict_it11(theta_deg: float, inv_it11: float) -> float:
    theta = math.radians(theta_deg)
    x = math.cos(theta)
    s = math.sin(theta)
    p1 = _legendre_p(1, x)
    p2 = _legendre_p(2, x)
    value = 1.20 * inv_it11 * s * (-1.20 * p2 - 0.20 * p1)
    return _clamp(value, -1.0, 1.0)


def _predict_dsigma(theta_deg: float, u_norm: float, inv_t20: float, inv_t22: float) -> float:
    theta = math.radians(theta_deg)
    x = math.cos(theta)
    p2 = _legendre_p(2, x)
    p4 = _legendre_p(4, x)
    sigma_shape = 1.10 * inv_t20 * p2 + 0.60 * inv_t22 * p4
    return max(u_norm * math.exp(sigma_shape), 1e-14)


def _predict_observables_from_u(
    u_eff: Dict[str, complex],
    exp_it11: Dict[str, List[float]],
    exp_dsigma: Dict[str, object],
    dsigma_output_unit: str = UNIT_MB_PER_SR,
) -> Dict[str, object]:
    u00 = u_eff["u00"]
    u01 = u_eff["u01"]
    u10 = u_eff["u10"]
    u11 = u_eff["u11"]

    u_norm = abs(u00) ** 2 + abs(u01) ** 2 + abs(u10) ** 2 + abs(u11) ** 2
    u_norm = max(u_norm, 1e-18)

    phase_indicator = ((u00 + u11).conjugate() * (u01 - u10)).imag
    phase_sign = 1.0 if phase_indicator >= 0.0 else -1.0

    inv_it11 = phase_indicator / u_norm
    inv_t20 = (abs(u00) ** 2 + abs(u11) ** 2 - abs(u01) ** 2 - abs(u10) ** 2) / u_norm
    inv_t21 = ((u00 - u11).conjugate() * (u01 + u10)).real / u_norm
    inv_t22 = (u00.conjugate() * u11 - u01.conjugate() * u10).real / u_norm

    ay_pred = [_predict_it11(angle, inv_it11) for angle in exp_it11["angles"]]
    sigma_mb = [_predict_dsigma(angle, u_norm, inv_t20, inv_t22) for angle in exp_dsigma["angles"]]
    sigma_pred = [
        convert_dsigma_value(value, UNIT_MB_PER_SR, dsigma_output_unit) for value in sigma_mb
    ]

    ay_metrics = _residual_metrics(ay_pred, exp_it11["values"])
    dsigma_metrics = _residual_metrics(sigma_pred, exp_dsigma["values"])
    dsigma_rel_rmse = _relative_rmse(sigma_pred, exp_dsigma["values"])

    return {
        "u_norm": u_norm,
        "phase_indicator": phase_indicator,
        "phase_sign": phase_sign,
        "inv_it11": inv_it11,
        "inv_t20": inv_t20,
        "inv_t21": inv_t21,
        "inv_t22": inv_t22,
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


def _unit_to_axis_suffix(unit: str) -> str:
    normalized = normalize_dsigma_unit(unit)
    if normalized == UNIT_MB_PER_SR:
        return "mb/sr"
    return "fm^2/sr"


def _unit_to_csv_suffix(unit: str) -> str:
    normalized = normalize_dsigma_unit(unit)
    if normalized == UNIT_MB_PER_SR:
        return "mb_per_sr"
    return "fm2_per_sr"


def _write_comparison_svg(
    path: Path,
    *,
    title: str,
    x_label: str,
    y_label: str,
    angles: Sequence[float],
    exp_vals: Sequence[float],
    pred_vals: Sequence[float],
    use_log_y: bool = False,
) -> None:
    if len(angles) != len(exp_vals) or len(angles) != len(pred_vals):
        raise ValueError("angles/exp/pred length mismatch")
    if not angles:
        raise ValueError("empty series")

    width = 980
    height = 640
    left = 90
    right = 30
    top = 60
    bottom = 90
    plot_w = width - left - right
    plot_h = height - top - bottom

    x_min = min(angles)
    x_max = max(angles)
    if x_max <= x_min:
        x_max = x_min + 1.0

    y_exp = list(exp_vals)
    y_pred = list(pred_vals)
    if use_log_y:
        y_exp = [math.log10(max(val, 1e-12)) for val in y_exp]
        y_pred = [math.log10(max(val, 1e-12)) for val in y_pred]

    y_min = min(min(y_exp), min(y_pred))
    y_max = max(max(y_exp), max(y_pred))
    if y_max <= y_min:
        y_max = y_min + 1.0
    y_pad = 0.08 * (y_max - y_min)
    y_min -= y_pad
    y_max += y_pad

    def x_to_px(x: float) -> float:
        return left + (x - x_min) * plot_w / (x_max - x_min)

    def y_to_px(y: float) -> float:
        return top + (y_max - y) * plot_h / (y_max - y_min)

    def _make_polyline(xs: Sequence[float], ys: Sequence[float]) -> str:
        pts = [f"{x_to_px(x):.3f},{y_to_px(y):.3f}" for x, y in zip(xs, ys)]
        return " ".join(pts)

    grid_lines: List[str] = []
    tick_labels: List[str] = []

    num_y_ticks = 6
    for i in range(num_y_ticks + 1):
        frac = i / num_y_ticks
        y_val = y_min + frac * (y_max - y_min)
        y_px = y_to_px(y_val)
        grid_lines.append(
            f"<line x1='{left}' y1='{y_px:.3f}' x2='{left + plot_w}' y2='{y_px:.3f}' "
            "stroke='#e5e7eb' stroke-width='1'/>"
        )
        label = f"{(10 ** y_val):.3g}" if use_log_y else f"{y_val:.3f}"
        tick_labels.append(
            f"<text x='{left - 10}' y='{y_px + 4:.3f}' text-anchor='end' "
            "font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>"
            f"{label}</text>"
        )

    num_x_ticks = 8
    for i in range(num_x_ticks + 1):
        frac = i / num_x_ticks
        x_val = x_min + frac * (x_max - x_min)
        x_px = x_to_px(x_val)
        grid_lines.append(
            f"<line x1='{x_px:.3f}' y1='{top}' x2='{x_px:.3f}' y2='{top + plot_h}' "
            "stroke='#f3f4f6' stroke-width='1'/>"
        )
        tick_labels.append(
            f"<text x='{x_px:.3f}' y='{top + plot_h + 24}' text-anchor='middle' "
            "font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>"
            f"{x_val:.1f}</text>"
        )

    model_poly = _make_polyline(angles, y_pred)
    exp_points = []
    for x, y in zip(angles, y_exp):
        exp_points.append(
            f"<circle cx='{x_to_px(x):.3f}' cy='{y_to_px(y):.3f}' r='3.2' "
            "fill='#1f2937' stroke='white' stroke-width='0.8'/>"
        )

    y_axis_text = y_label + (" (log10 scale)" if use_log_y else "")
    svg = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect x='0' y='0' width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2:.1f}' y='32' text-anchor='middle' font-family='DejaVu Sans, Arial, sans-serif' "
        "font-size='22' font-weight='700' fill='#111827'>"
        f"{title}</text>",
        *grid_lines,
        f"<rect x='{left}' y='{top}' width='{plot_w}' height='{plot_h}' fill='none' stroke='#111827' stroke-width='1.2'/>",
        f"<polyline fill='none' stroke='#dc2626' stroke-width='2.4' points='{model_poly}'/>",
        *exp_points,
        *tick_labels,
        f"<text x='{left + plot_w/2:.1f}' y='{height - 28}' text-anchor='middle' font-family='DejaVu Sans, Arial, sans-serif' "
        "font-size='16' fill='#111827'>"
        f"{x_label}</text>",
        f"<text x='24' y='{top + plot_h/2:.1f}' text-anchor='middle' transform='rotate(-90 24 {top + plot_h/2:.1f})' "
        "font-family='DejaVu Sans, Arial, sans-serif' font-size='16' fill='#111827'>"
        f"{y_axis_text}</text>",
        f"<line x1='{left + 18}' y1='{top + 18}' x2='{left + 64}' y2='{top + 18}' stroke='#dc2626' stroke-width='2.4'/>",
        f"<text x='{left + 70}' y='{top + 22}' font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>model</text>",
        f"<circle cx='{left + 38}' cy='{top + 40}' r='3.2' fill='#1f2937' stroke='white' stroke-width='0.8'/>",
        f"<text x='{left + 70}' y='{top + 44}' font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>experiment</text>",
        "</svg>",
    ]
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def run_solver_if_requested(work_dir: Path, regenerate: bool, target_tlab_mev: float) -> None:
    if not regenerate:
        return

    cmd = [
        "python3",
        "examples/deuteron_proton_Ay.py",
        "--work-dir",
        str(work_dir),
        "--target-tlab-mev",
        str(target_tlab_mev),
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
    dsigma_unit = str(summary["units"]["dsigma_output"])  # type: ignore[index]
    dsigma_unit_text = _unit_to_axis_suffix(dsigma_unit)

    lines: List[str] = []
    lines.append("d + p validation from Tic-tac U-matrix")
    lines.append("======================================")
    lines.append(f"status: {summary['status']}")
    lines.append(f"reason: {summary['status_reason']}")
    lines.append(f"target_tlab_mev: {summary['target_tlab_mev']:.3f}")
    lines.append(f"best_tlab_mev: {best['tlab_mev']:.3f}")
    lines.append(f"abs_delta_tlab_mev: {best['abs_delta_tlab_mev']:.3f}")
    lines.append("")

    lines.append("Best-energy metrics:")
    lines.append(f"  iT11 mae: {best['it11_mae']:.6f}")
    lines.append(f"  iT11 rmse: {best['it11_rmse']:.6f}")
    lines.append(f"  iT11 max_abs_error: {best['it11_max_abs_error']:.6f}")
    lines.append(f"  dSigma mae ({dsigma_unit_text}): {best['dsigma_mae']:.6f}")
    lines.append(f"  dSigma rmse ({dsigma_unit_text}): {best['dsigma_rmse']:.6f}")
    lines.append(f"  dSigma relative rmse: {best['dsigma_rel_rmse']:.6f}")
    lines.append("")

    lines.append("Model settings:")
    lines.append("  observable model: reduced-U deterministic map")
    lines.append("  dSigma: u_norm * exp(1.10*inv_t20*P2 + 0.60*inv_t22*P4)")
    lines.append("  iT11: clamp(1.20*inv_it11*sin(theta)*(-1.20*P2 - 0.20*P1), -1, 1)")
    lines.append(f"  dSigma base unit (model internal): {_unit_to_axis_suffix(UNIT_MB_PER_SR)}")
    lines.append(f"  dSigma reported unit: {dsigma_unit_text}")
    lines.append(f"  phase_sign: {best['phase_sign']:.1f}")
    lines.append(f"  u_norm: {best['u_norm']:.6e}")
    lines.append(f"  inv_it11: {best['inv_it11']:.6e}")
    lines.append(f"  inv_t20: {best['inv_t20']:.6e}")
    lines.append(f"  inv_t22: {best['inv_t22']:.6e}")
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
            f"Tlab={item['tlab_mev']:.3f} MeV, "
            f"iT11_rmse={item['it11_rmse']:.6f}, "
            f"dSigma_rel_rmse={item['dsigma_rel_rmse']:.6f}"
        )

    lines.append("")
    lines.append("Note:")
    lines.append("  This workflow predicts observables from solver-produced U-matrix elements only.")
    lines.append("  Experimental data is used only for residual evaluation.")
    lines.append("  No interpolation or fitting to experimental curves is performed in prediction.")
    lines.append("  External references often call this benchmark '190 MeV/u'; this repository uses solver input Tlab [MeV].")

    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate Tic-tac solver output against a target Tlab benchmark in MeV")
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
    parser.add_argument(
        "--dsigma-unit",
        default=UNIT_MB_PER_SR,
        choices=list(SUPPORTED_DSIGMA_UNITS),
        help="Output unit for dSigma/dOmega in reports, CSV, and plots",
    )
    parser.add_argument("--target-tlab-mev", type=float, default=190.0)
    parser.add_argument("--regenerate", action="store_true", help="Run solver before comparison")
    parser.add_argument("--ay-rmse-pass", type=float, default=0.02)
    parser.add_argument("--dsigma-rel-rmse-pass", type=float, default=0.05)
    parser.add_argument("--energy-delta-pass", type=float, default=40.0)
    return parser


def _serialize_complex(z: complex) -> Dict[str, float]:
    return {"re": z.real, "im": z.imag}


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    dsigma_output_unit = normalize_dsigma_unit(args.dsigma_unit)

    # Resolve all IO roots first so the rest of the script only handles absolute paths.
    work_dir = (root / args.work_dir).resolve()
    solver_out_dir = (root / args.solver_out_dir).resolve()
    tensor_data = (root / args.tensor_data).resolve()
    dsigma_data = (root / args.dsigma_data).resolve()

    work_dir.mkdir(parents=True, exist_ok=True)
    run_solver_if_requested(work_dir, args.regenerate, args.target_tlab_mev)

    # Solver branch: U files -> parsed channel points.
    u_files = find_latest_solver_u_files(solver_out_dir)
    if not u_files:
        raise RuntimeError(f"No U_PW_elements files found in {solver_out_dir}")

    points: List[SolverChannelPoint] = []
    for path in u_files:
        points.extend(parse_u_file(path))

    if not points:
        raise RuntimeError("No valid U-matrix rows parsed from solver output")

    # Experiment branch: raw measurement tables -> typed arrays.
    exp_it11 = read_experimental_iT11(tensor_data)
    exp_dsigma = read_experimental_dsigma(dsigma_data, target_unit=dsigma_output_unit)

    # Join branch: combine solver channels by energy and evaluate every available energy point.
    combined = combine_channels_by_energy(points)
    if not combined:
        raise RuntimeError("Failed to combine solver channels into energy points")

    enriched: List[Dict[str, object]] = []
    for item in combined:
        fit = _predict_observables_from_u(
            item["u_eff"],
            exp_it11,
            exp_dsigma,
            dsigma_output_unit=dsigma_output_unit,
        )

        it11_metrics = fit["it11_metrics"]
        dsigma_metrics = fit["dsigma_metrics"]

        enriched_item: Dict[str, object] = {
            "tlab_mev": item["tlab"],
            "ecm_weighted_mev": item["ecm_weighted"],
            "abs_delta_tlab_mev": abs(float(item["tlab"]) - args.target_tlab_mev),
            "it11_mae": it11_metrics["mae"],
            "it11_rmse": it11_metrics["rmse"],
            "it11_max_abs_error": it11_metrics["max_abs_error"],
            "dsigma_mae": dsigma_metrics["mae"],
            "dsigma_rmse": dsigma_metrics["rmse"],
            "dsigma_max_abs_error": dsigma_metrics["max_abs_error"],
            "dsigma_rel_rmse": fit["dsigma_rel_rmse"],
            "phase_sign": fit["phase_sign"],
            "phase_indicator": fit["phase_indicator"],
            "u_norm": fit["u_norm"],
            "inv_it11": fit["inv_it11"],
            "inv_t20": fit["inv_t20"],
            "inv_t21": fit["inv_t21"],
            "inv_t22": fit["inv_t22"],
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

    # Select closest available solver energy to the requested target and emit artifacts.
    best = min(enriched, key=lambda item: float(item["abs_delta_tlab_mev"]))

    pass_flag = (
        float(best["it11_rmse"]) <= args.ay_rmse_pass
        and float(best["dsigma_rel_rmse"]) <= args.dsigma_rel_rmse_pass
        and float(best["abs_delta_tlab_mev"]) <= args.energy_delta_pass
    )

    ay_csv = work_dir / "best_energy_iT11_curve.csv"
    dsigma_csv = work_dir / "best_energy_dsigma_curve.csv"
    ay_svg = work_dir / "best_energy_iT11_comparison.svg"
    dsigma_svg = work_dir / "best_energy_dsigma_comparison.svg"
    _write_curve_csv(
        ay_csv,
        best["it11_curve"]["angles_deg"],
        best["it11_curve"]["exp"],
        best["it11_curve"]["pred"],
        "iT11_exp",
        "iT11_model",
    )
    dsigma_csv_suffix = _unit_to_csv_suffix(dsigma_output_unit)
    _write_curve_csv(
        dsigma_csv,
        best["dsigma_curve"]["angles_deg"],
        best["dsigma_curve"]["exp"],
        best["dsigma_curve"]["pred"],
        f"dsigma_exp_{dsigma_csv_suffix}",
        f"dsigma_model_{dsigma_csv_suffix}",
    )
    _write_comparison_svg(
        ay_svg,
        title=f"iT11 Comparison (best Tlab={best['tlab_mev']:.3f} MeV)",
        x_label="theta_cm (deg)",
        y_label="iT11",
        angles=best["it11_curve"]["angles_deg"],
        exp_vals=best["it11_curve"]["exp"],
        pred_vals=best["it11_curve"]["pred"],
        use_log_y=False,
    )
    _write_comparison_svg(
        dsigma_svg,
        title=f"dSigma/dOmega Comparison (best Tlab={best['tlab_mev']:.3f} MeV)",
        x_label="theta_cm (deg)",
        y_label=f"dSigma/dOmega [{_unit_to_axis_suffix(dsigma_output_unit)}]",
        angles=best["dsigma_curve"]["angles_deg"],
        exp_vals=best["dsigma_curve"]["exp"],
        pred_vals=best["dsigma_curve"]["pred"],
        use_log_y=True,
    )

    summary = {
        "status": "pass" if pass_flag else "fail",
        "status_reason": (
            "best solver-based predicted observables are within thresholds"
            if pass_flag
            else "best solver-based predicted observables exceed thresholds"
        ),
        "target_tlab_mev": args.target_tlab_mev,
        "thresholds": {
            "ay_rmse_pass": args.ay_rmse_pass,
            "dsigma_rel_rmse_pass": args.dsigma_rel_rmse_pass,
            "energy_delta_pass": args.energy_delta_pass,
        },
        "model_settings": {
            "observable_model": "reduced-U deterministic map",
            "dsigma_formula": "u_norm * exp(1.10*inv_t20*P2 + 0.60*inv_t22*P4)",
            "it11_formula": "clamp(1.20*inv_it11*sin(theta)*(-1.20*P2 - 0.20*P1), -1, 1)",
        },
        "units": {
            "dsigma_input_detected": exp_dsigma["source_unit"],
            "dsigma_output": exp_dsigma["output_unit"],
            "dsigma_model_base": UNIT_MB_PER_SR,
            "it11": "dimensionless",
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
            "best_it11_svg": str(ay_svg),
            "best_dsigma_svg": str(dsigma_svg),
        },
    }

    tlab_label = format_tlab_label(args.target_tlab_mev)
    json_file = work_dir / f"solver_validation_tlab_{tlab_label}.json"
    txt_file = work_dir / f"solver_validation_tlab_{tlab_label}.txt"

    json_file.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    txt_file.write_text(build_report_text(summary), encoding="utf-8")

    print(build_report_text(summary), end="")
    print(f"json_report: {json_file}")
    print(f"text_report: {txt_file}")
    print(f"it11_curve_csv: {ay_csv}")
    print(f"dsigma_curve_csv: {dsigma_csv}")
    print(f"it11_curve_svg: {ay_svg}")
    print(f"dsigma_curve_svg: {dsigma_svg}")

    return 0 if pass_flag else 2


if __name__ == "__main__":
    raise SystemExit(main())
