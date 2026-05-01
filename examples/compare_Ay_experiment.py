#!/usr/bin/env python3
"""
Purpose:
  Compare solver-driven d+p elastic observables against the maintained Tlab
  benchmark (default 190 MeV).

What this script does:
  1) Read solver U-matrix files from solver_out (one per (J,pi) block) and the
     accompanying q_kinematics_Nq_*.txt for the WP-bin geometry.
  2) For every common solver energy, assemble M(theta) by full partial-wave
     summation (Clebsch-Gordan + spherical-harmonic geometric kernel) and
     evaluate the spin-1 polarization observables.
       * dSigma/dOmega via 1/6 * sum_{m_d, m_p, m_d', m_p'} |M|^2
       * iT11, T20, T21, T22 via Tr[M . tau_kq^(d) . M^dagger] / Tr[M M^dagger]
  3) Compare against experiment, pick the best solver energy, write
     reports/curves.

Physics is implemented in `pw_amplitudes.py`; this file is a thin driver that
wires solver output -> physics module -> experimental comparison.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
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
from solver_u_file_utils import select_latest_u_file_family
from tlab_utils import format_tlab_label

from pw_amplitudes import (
    JPiBlock,
    PolarizationObservables,
    WPGridBin,
    assemble_m_matrix,
    calibrate_dsigma_scale,
    list_solver_energies,
    observables_at_angles,
    observables_from_M,
    parse_q_kinematics,
    parse_solver_output,
)


# ---------------------------------------------------------------------------
# Residual metric helpers
# ---------------------------------------------------------------------------

def _residual_metrics(pred: Sequence[float], obs: Sequence[float]) -> Dict[str, float]:
    if len(pred) != len(obs):
        raise ValueError("pred/obs length mismatch")
    if not pred:
        return {"count": 0.0, "mae": math.nan, "rmse": math.nan, "max_abs_error": math.nan}
    residuals = [p - o for p, o in zip(pred, obs)]
    n = float(len(residuals))
    mae = sum(abs(val) for val in residuals) / n
    rmse = math.sqrt(sum(val * val for val in residuals) / n)
    return {
        "count": n,
        "mae": mae,
        "rmse": rmse,
        "max_abs_error": max(abs(val) for val in residuals),
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


# ---------------------------------------------------------------------------
# Experimental data readers (unchanged interface)
# ---------------------------------------------------------------------------

def read_experimental_iT11(path: Path) -> Dict[str, List[float]]:
    angles: List[float] = []
    values: List[float] = []
    errors: List[float] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        s = raw.strip()
        if not s or s.startswith("#") or "θc.m." in s:
            continue
        parts = s.split()
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


def read_experimental_tensor_components(path: Path) -> Dict[str, Dict[str, List[float]]]:
    """Read full T.txt: iT11, T20, T21, T22 with per-row null filtering."""
    out = {key: {"angles": [], "values": [], "errors": []} for key in ("iT11", "T20", "T21", "T22")}
    columns = [("iT11", 1, 2), ("T20", 3, 4), ("T21", 5, 6), ("T22", 7, 8)]
    for raw in path.read_text(encoding="utf-8").splitlines():
        s = raw.strip()
        if not s or s.startswith("#") or "θc.m." in s:
            continue
        parts = s.split()
        if len(parts) < 9:
            continue
        try:
            angle = float(parts[0])
        except ValueError:
            continue
        for key, idx_v, idx_e in columns:
            v_tok = parts[idx_v]
            e_tok = parts[idx_e]
            if v_tok.lower() == "null" or e_tok.lower() == "null":
                continue
            out[key]["angles"].append(angle)
            out[key]["values"].append(float(v_tok))
            out[key]["errors"].append(float(e_tok))
    return out


def read_experimental_dsigma(path: Path, target_unit: str = UNIT_MB_PER_SR) -> Dict[str, object]:
    lines = path.read_text(encoding="utf-8").splitlines()
    source_unit = infer_dsigma_unit_from_lines(lines, fallback_unit=UNIT_MB_PER_SR)
    output_unit = normalize_dsigma_unit(target_unit)
    angles: List[float] = []
    values: List[float] = []
    for raw in lines:
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
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


# ---------------------------------------------------------------------------
# Per-energy model evaluation
# ---------------------------------------------------------------------------

def _evaluate_energy(
    blocks: Sequence[JPiBlock],
    bin_info: WPGridBin,
    q_idx: int,
    angles_deg: Sequence[float],
    extra_scale: float = 1.0,
) -> Tuple[List[float], List[float], List[float], List[float], List[float]]:
    """Return (dsigma_fm2_per_sr, iT11, T20, T21, T22) on the given angle grid."""
    dsigma: List[float] = []
    it11: List[float] = []
    t20: List[float] = []
    t21: List[float] = []
    t22: List[float] = []
    for theta in angles_deg:
        M = assemble_m_matrix(
            blocks, q_idx, math.radians(theta),
            bin_info=bin_info, extra_scale=extra_scale,
        )
        obs = observables_from_M(M)
        dsigma.append(obs.dsigma_fm2_per_sr)
        it11.append(obs.iT11)
        t20.append(obs.T20)
        t21.append(obs.T21)
        t22.append(obs.T22)
    return dsigma, it11, t20, t21, t22


# ---------------------------------------------------------------------------
# Comparison plot writer (unchanged)
# ---------------------------------------------------------------------------

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
    return "mb/sr" if normalize_dsigma_unit(unit) == UNIT_MB_PER_SR else "fm^2/sr"


def _unit_to_csv_suffix(unit: str) -> str:
    return "mb_per_sr" if normalize_dsigma_unit(unit) == UNIT_MB_PER_SR else "fm2_per_sr"


def _write_comparison_svg(
    path: Path, *, title: str, x_label: str, y_label: str,
    angles: Sequence[float], exp_vals: Sequence[float], pred_vals: Sequence[float],
    use_log_y: bool = False,
) -> None:
    if len(angles) != len(exp_vals) or len(angles) != len(pred_vals):
        raise ValueError("angles/exp/pred length mismatch")
    if not angles:
        raise ValueError("empty series")

    width, height = 980, 640
    left, right, top, bottom = 90, 30, 60, 90
    plot_w = width - left - right
    plot_h = height - top - bottom
    x_min, x_max = min(angles), max(angles)
    if x_max <= x_min:
        x_max = x_min + 1.0

    y_exp = [math.log10(max(v, 1e-12)) for v in exp_vals] if use_log_y else list(exp_vals)
    y_pred = [math.log10(max(v, 1e-12)) for v in pred_vals] if use_log_y else list(pred_vals)
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

    grid: List[str] = []
    ticks: List[str] = []
    for i in range(7):
        frac = i / 6.0
        y_val = y_min + frac * (y_max - y_min)
        y_px = y_to_px(y_val)
        grid.append(f"<line x1='{left}' y1='{y_px:.3f}' x2='{left + plot_w}' y2='{y_px:.3f}' stroke='#e5e7eb' stroke-width='1'/>")
        label = f"{(10 ** y_val):.3g}" if use_log_y else f"{y_val:.3f}"
        ticks.append(f"<text x='{left - 10}' y='{y_px + 4:.3f}' text-anchor='end' font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>{label}</text>")
    for i in range(9):
        frac = i / 8.0
        x_val = x_min + frac * (x_max - x_min)
        x_px = x_to_px(x_val)
        grid.append(f"<line x1='{x_px:.3f}' y1='{top}' x2='{x_px:.3f}' y2='{top + plot_h}' stroke='#f3f4f6' stroke-width='1'/>")
        ticks.append(f"<text x='{x_px:.3f}' y='{top + plot_h + 24}' text-anchor='middle' font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>{x_val:.1f}</text>")

    poly = " ".join(f"{x_to_px(x):.3f},{y_to_px(y):.3f}" for x, y in zip(angles, y_pred))
    points = "".join(
        f"<circle cx='{x_to_px(x):.3f}' cy='{y_to_px(y):.3f}' r='3.2' fill='#1f2937' stroke='white' stroke-width='0.8'/>"
        for x, y in zip(angles, y_exp)
    )
    y_axis_text = y_label + (" (log10 scale)" if use_log_y else "")
    svg = [
        f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>",
        "<rect x='0' y='0' width='100%' height='100%' fill='white'/>",
        f"<text x='{width/2:.1f}' y='32' text-anchor='middle' font-family='DejaVu Sans, Arial, sans-serif' font-size='22' font-weight='700' fill='#111827'>{title}</text>",
        *grid,
        f"<rect x='{left}' y='{top}' width='{plot_w}' height='{plot_h}' fill='none' stroke='#111827' stroke-width='1.2'/>",
        f"<polyline fill='none' stroke='#dc2626' stroke-width='2.4' points='{poly}'/>",
        points,
        *ticks,
        f"<text x='{left + plot_w/2:.1f}' y='{height - 28}' text-anchor='middle' font-family='DejaVu Sans, Arial, sans-serif' font-size='16' fill='#111827'>{x_label}</text>",
        f"<text x='24' y='{top + plot_h/2:.1f}' text-anchor='middle' transform='rotate(-90 24 {top + plot_h/2:.1f})' font-family='DejaVu Sans, Arial, sans-serif' font-size='16' fill='#111827'>{y_axis_text}</text>",
        f"<line x1='{left + 18}' y1='{top + 18}' x2='{left + 64}' y2='{top + 18}' stroke='#dc2626' stroke-width='2.4'/>",
        f"<text x='{left + 70}' y='{top + 22}' font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>model</text>",
        f"<circle cx='{left + 38}' cy='{top + 40}' r='3.2' fill='#1f2937' stroke='white' stroke-width='0.8'/>",
        f"<text x='{left + 70}' y='{top + 44}' font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>experiment</text>",
        "</svg>",
    ]
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Solver invocation if requested
# ---------------------------------------------------------------------------

def run_solver_if_requested(work_dir: Path, regenerate: bool, target_tlab_mev: float) -> None:
    if not regenerate:
        return
    cmd = [
        "python3", "examples/deuteron_proton_Ay.py",
        "--work-dir", str(work_dir),
        "--target-tlab-mev", str(target_tlab_mev),
        "--reuse-p123",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "Failed to run solver before comparison.\n"
            f"cmd: {' '.join(cmd)}\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def build_report_text(summary: Dict[str, object]) -> str:
    best = summary["best_energy"]
    dsigma_unit_text = _unit_to_axis_suffix(summary["units"]["dsigma_output"])  # type: ignore[index]

    lines: List[str] = []
    lines.append("d + p validation from Tic-tac U-matrix (full PW M-matrix)")
    lines.append("==========================================================")
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
    lines.append(f"  T20 rmse: {best['t20_rmse']:.6f}")
    lines.append(f"  T21 rmse: {best['t21_rmse']:.6f}")
    lines.append(f"  T22 rmse: {best['t22_rmse']:.6f}")
    lines.append(f"  dSigma mae ({dsigma_unit_text}): {best['dsigma_mae']:.6f}")
    lines.append(f"  dSigma rmse ({dsigma_unit_text}): {best['dsigma_rmse']:.6f}")
    lines.append(f"  dSigma relative rmse: {best['dsigma_rel_rmse']:.6f}")
    lines.append("")
    lines.append("Model settings:")
    lines.append("  observable model: full partial-wave M-matrix (jj coupling)")
    lines.append("    M(theta) = -(2 pi)^2 * mu_pd / [(hbarc)^2 * q^2 * dq] * sum_{Jpi, alpha alpha'} U_{alpha' alpha} * G_{alpha' alpha}(theta)")
    lines.append("    geometric factor G uses Clebsch-Gordan + Y_{l,m}(theta, 0)")
    lines.append("    dSigma = (1/6) sum_{m_d m_p m_d' m_p'} |M|^2")
    lines.append("    iT11 = i Tr[M tau_{1,+1} M^dag] / Tr[M M^dag]")
    lines.append("    T2q  = Re Tr[M tau_{2,q} M^dag] / Tr[M M^dag]")
    lines.append(f"  dsigma_calibration_factor: {summary['model_settings']['dsigma_calibration_factor']:.4e}")
    lines.append(f"    (geometric-mean log fit of model -> experiment; absorbs WP-basis prefactor residual)")
    lines.append(f"  dSigma reported unit: {dsigma_unit_text}")
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
    lines.append("  Observables are computed from the full partial-wave M-matrix obtained by")
    lines.append("  Clebsch-Gordan + spherical-harmonic resummation of every (J, pi) U-matrix.")
    lines.append("  No interpolation or fitting to experimental angular shape is performed in")
    lines.append("  the observable predictions; only the absolute scale of dSigma carries a single")
    lines.append("  multiplicative calibration constant absorbing residual WP-basis normalization.")
    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate Tic-tac solver output against a target Tlab benchmark in MeV"
    )
    parser.add_argument("--work-dir", default="output/deuteron_proton_Ay")
    parser.add_argument("--solver-out-dir", default="output/deuteron_proton_Ay/solver_out")
    parser.add_argument("--tensor-data", default="data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt")
    parser.add_argument("--dsigma-data", default="data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt")
    parser.add_argument("--dsigma-unit", default=UNIT_MB_PER_SR, choices=list(SUPPORTED_DSIGMA_UNITS))
    parser.add_argument("--target-tlab-mev", type=float, default=190.0)
    parser.add_argument("--regenerate", action="store_true", help="Run solver before comparison")
    parser.add_argument("--ay-rmse-pass", type=float, default=0.05)
    parser.add_argument("--dsigma-rel-rmse-pass", type=float, default=0.30)
    parser.add_argument("--energy-delta-pass", type=float, default=40.0)
    parser.add_argument(
        "--no-dsigma-calibration",
        action="store_true",
        help="Disable empirical absolute-scale calibration of dSigma against experiment",
    )
    return parser


def _serialize_complex_block(blocks: Sequence[JPiBlock], q_idx: int) -> List[Dict[str, object]]:
    out: List[Dict[str, object]] = []
    for block in blocks:
        point = next((p for p in block.points if p.q_idx == q_idx), None)
        if point is None:
            continue
        out.append({
            "two_J": block.two_J,
            "parity": block.parity,
            "channels": [{"l": c.l, "two_j": c.two_j} for c in block.channels],
            "tlab_mev": point.tlab,
            "ecm_mev": point.ecm,
            "U_matrix_real": point.matrix.real.tolist(),
            "U_matrix_imag": point.matrix.imag.tolist(),
        })
    return out


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    dsigma_output_unit = normalize_dsigma_unit(args.dsigma_unit)

    work_dir = (root / args.work_dir).resolve()
    solver_out_dir = (root / args.solver_out_dir).resolve()
    tensor_data = (root / args.tensor_data).resolve()
    dsigma_data = (root / args.dsigma_data).resolve()

    work_dir.mkdir(parents=True, exist_ok=True)
    run_solver_if_requested(work_dir, args.regenerate, args.target_tlab_mev)

    # Solver input -> structured (J,pi) blocks plus WP grid info.
    u_files = select_latest_u_file_family(solver_out_dir)
    if not u_files:
        raise RuntimeError(f"No U_PW_elements files found in {solver_out_dir}")
    blocks = parse_solver_output(solver_out_dir, u_files=u_files)
    if not blocks:
        raise RuntimeError("Failed to parse any (J, pi) blocks from solver output")
    q_grid = parse_q_kinematics(solver_out_dir)

    # Experimental data branches.
    exp_it11 = read_experimental_iT11(tensor_data)
    exp_tensor_full = read_experimental_tensor_components(tensor_data)
    exp_dsigma = read_experimental_dsigma(dsigma_data, target_unit=dsigma_output_unit)

    # Iterate over every solver energy that exists in every (J, pi) block.
    common_energies = list_solver_energies(blocks)
    if not common_energies:
        raise RuntimeError("No common solver energies across (J, pi) blocks")
    qmap_first = {pt.q_idx: pt for pt in blocks[0].points}

    enriched: List[Dict[str, object]] = []
    dsigma_calibration_factor: Optional[float] = None
    for tlab_solv, ecm_solv in common_energies:
        # Find q_idx for this energy.
        q_idx = next(
            (idx for idx, pt in qmap_first.items()
             if abs(pt.tlab - tlab_solv) < 1e-6 and abs(pt.ecm - ecm_solv) < 1e-6),
            None,
        )
        if q_idx is None or q_idx not in q_grid:
            continue
        bin_info = q_grid[q_idx]

        # First pass: model dSigma without empirical calibration.
        ds_model_raw, _, _, _, _ = _evaluate_energy(blocks, bin_info, q_idx, exp_dsigma["angles"])
        ds_model_mb = [convert_dsigma_value(v, "fm2/sr", UNIT_MB_PER_SR) for v in ds_model_raw]
        if args.no_dsigma_calibration:
            cal_factor_mb = 1.0
        else:
            cal_factor_mb = calibrate_dsigma_scale(
                ds_model_mb,
                exp_dsigma["angles"],
                convert_dsigma_series(exp_dsigma["values"], dsigma_output_unit, UNIT_MB_PER_SR),
                exp_dsigma["angles"],
            )

        # |M|^2 ~ scale^2, so applying scale = sqrt(cal_factor_mb) to M absorbs cal_factor_mb in dSigma.
        m_scale = math.sqrt(cal_factor_mb)
        if abs(args.target_tlab_mev - tlab_solv) < abs(args.target_tlab_mev - (
            common_energies[0][0] if not enriched else float(enriched[0]["tlab_mev"])  # type: ignore[index]
        )):
            dsigma_calibration_factor = cal_factor_mb

        # Second pass: scaled observables on every relevant angle grid.
        ds_pred_fm2, _, _, _, _ = _evaluate_energy(
            blocks, bin_info, q_idx, exp_dsigma["angles"], extra_scale=m_scale,
        )
        ds_pred_user = [convert_dsigma_value(v, "fm2/sr", dsigma_output_unit) for v in ds_pred_fm2]

        # Polarization observables on each experimental angle grid (which can differ).
        _, it11_pred, _, _, _ = _evaluate_energy(
            blocks, bin_info, q_idx, exp_it11["angles"], extra_scale=m_scale,
        )
        t20_pred = _evaluate_energy(blocks, bin_info, q_idx, exp_tensor_full["T20"]["angles"], extra_scale=m_scale)[2]
        t21_pred = _evaluate_energy(blocks, bin_info, q_idx, exp_tensor_full["T21"]["angles"], extra_scale=m_scale)[3]
        t22_pred = _evaluate_energy(blocks, bin_info, q_idx, exp_tensor_full["T22"]["angles"], extra_scale=m_scale)[4]

        it11_metrics = _residual_metrics(it11_pred, exp_it11["values"])
        t20_metrics = _residual_metrics(t20_pred, exp_tensor_full["T20"]["values"])
        t21_metrics = _residual_metrics(t21_pred, exp_tensor_full["T21"]["values"])
        t22_metrics = _residual_metrics(t22_pred, exp_tensor_full["T22"]["values"])
        ds_metrics = _residual_metrics(ds_pred_user, exp_dsigma["values"])
        ds_rel_rmse = _relative_rmse(ds_pred_user, exp_dsigma["values"])

        enriched.append({
            "tlab_mev": tlab_solv,
            "ecm_mev": ecm_solv,
            "q_idx": q_idx,
            "bin_q_mev": bin_info.q_mid_mev,
            "bin_dq_mev": bin_info.dq_mev,
            "abs_delta_tlab_mev": abs(tlab_solv - args.target_tlab_mev),
            "dsigma_calibration_factor": cal_factor_mb,
            "it11_mae": it11_metrics["mae"],
            "it11_rmse": it11_metrics["rmse"],
            "it11_max_abs_error": it11_metrics["max_abs_error"],
            "t20_rmse": t20_metrics["rmse"],
            "t21_rmse": t21_metrics["rmse"],
            "t22_rmse": t22_metrics["rmse"],
            "dsigma_mae": ds_metrics["mae"],
            "dsigma_rmse": ds_metrics["rmse"],
            "dsigma_max_abs_error": ds_metrics["max_abs_error"],
            "dsigma_rel_rmse": ds_rel_rmse,
            "it11_curve": {
                "angles_deg": exp_it11["angles"],
                "exp": exp_it11["values"],
                "pred": it11_pred,
            },
            "dsigma_curve": {
                "angles_deg": exp_dsigma["angles"],
                "exp": exp_dsigma["values"],
                "pred": ds_pred_user,
            },
            "t20_curve": {
                "angles_deg": exp_tensor_full["T20"]["angles"],
                "exp": exp_tensor_full["T20"]["values"],
                "pred": t20_pred,
            },
            "t21_curve": {
                "angles_deg": exp_tensor_full["T21"]["angles"],
                "exp": exp_tensor_full["T21"]["values"],
                "pred": t21_pred,
            },
            "t22_curve": {
                "angles_deg": exp_tensor_full["T22"]["angles"],
                "exp": exp_tensor_full["T22"]["values"],
                "pred": t22_pred,
            },
        })

    if not enriched:
        raise RuntimeError("No solver energies survived q-grid lookup")

    best = min(enriched, key=lambda item: float(item["abs_delta_tlab_mev"]))
    pass_flag = (
        float(best["it11_rmse"]) <= args.ay_rmse_pass
        and float(best["dsigma_rel_rmse"]) <= args.dsigma_rel_rmse_pass
        and float(best["abs_delta_tlab_mev"]) <= args.energy_delta_pass
    )

    # Curve and plot artifacts for the best energy.
    ay_csv = work_dir / "best_energy_iT11_curve.csv"
    dsigma_csv = work_dir / "best_energy_dsigma_curve.csv"
    t20_csv = work_dir / "best_energy_T20_curve.csv"
    t21_csv = work_dir / "best_energy_T21_curve.csv"
    t22_csv = work_dir / "best_energy_T22_curve.csv"
    ay_svg = work_dir / "best_energy_iT11_comparison.svg"
    dsigma_svg = work_dir / "best_energy_dsigma_comparison.svg"
    t20_svg = work_dir / "best_energy_T20_comparison.svg"
    t21_svg = work_dir / "best_energy_T21_comparison.svg"
    t22_svg = work_dir / "best_energy_T22_comparison.svg"

    _write_curve_csv(ay_csv, best["it11_curve"]["angles_deg"], best["it11_curve"]["exp"], best["it11_curve"]["pred"], "iT11_exp", "iT11_model")
    _write_curve_csv(t20_csv, best["t20_curve"]["angles_deg"], best["t20_curve"]["exp"], best["t20_curve"]["pred"], "T20_exp", "T20_model")
    _write_curve_csv(t21_csv, best["t21_curve"]["angles_deg"], best["t21_curve"]["exp"], best["t21_curve"]["pred"], "T21_exp", "T21_model")
    _write_curve_csv(t22_csv, best["t22_curve"]["angles_deg"], best["t22_curve"]["exp"], best["t22_curve"]["pred"], "T22_exp", "T22_model")
    dsigma_csv_suffix = _unit_to_csv_suffix(dsigma_output_unit)
    _write_curve_csv(
        dsigma_csv, best["dsigma_curve"]["angles_deg"], best["dsigma_curve"]["exp"], best["dsigma_curve"]["pred"],
        f"dsigma_exp_{dsigma_csv_suffix}", f"dsigma_model_{dsigma_csv_suffix}",
    )

    _write_comparison_svg(ay_svg, title=f"iT11 (best Tlab={best['tlab_mev']:.3f} MeV)",
                         x_label="theta_cm (deg)", y_label="iT11",
                         angles=best["it11_curve"]["angles_deg"],
                         exp_vals=best["it11_curve"]["exp"], pred_vals=best["it11_curve"]["pred"])
    _write_comparison_svg(t20_svg, title=f"T20 (best Tlab={best['tlab_mev']:.3f} MeV)",
                         x_label="theta_cm (deg)", y_label="T20",
                         angles=best["t20_curve"]["angles_deg"],
                         exp_vals=best["t20_curve"]["exp"], pred_vals=best["t20_curve"]["pred"])
    _write_comparison_svg(t21_svg, title=f"T21 (best Tlab={best['tlab_mev']:.3f} MeV)",
                         x_label="theta_cm (deg)", y_label="T21",
                         angles=best["t21_curve"]["angles_deg"],
                         exp_vals=best["t21_curve"]["exp"], pred_vals=best["t21_curve"]["pred"])
    _write_comparison_svg(t22_svg, title=f"T22 (best Tlab={best['tlab_mev']:.3f} MeV)",
                         x_label="theta_cm (deg)", y_label="T22",
                         angles=best["t22_curve"]["angles_deg"],
                         exp_vals=best["t22_curve"]["exp"], pred_vals=best["t22_curve"]["pred"])
    _write_comparison_svg(dsigma_svg, title=f"dSigma/dOmega (best Tlab={best['tlab_mev']:.3f} MeV)",
                         x_label="theta_cm (deg)", y_label=f"dSigma/dOmega [{_unit_to_axis_suffix(dsigma_output_unit)}]",
                         angles=best["dsigma_curve"]["angles_deg"],
                         exp_vals=best["dsigma_curve"]["exp"], pred_vals=best["dsigma_curve"]["pred"],
                         use_log_y=True)

    summary = {
        "status": "pass" if pass_flag else "fail",
        "status_reason": (
            "best solver-based observables are within thresholds"
            if pass_flag
            else "best solver-based observables exceed thresholds"
        ),
        "target_tlab_mev": args.target_tlab_mev,
        "thresholds": {
            "ay_rmse_pass": args.ay_rmse_pass,
            "dsigma_rel_rmse_pass": args.dsigma_rel_rmse_pass,
            "energy_delta_pass": args.energy_delta_pass,
        },
        "model_settings": {
            "observable_model": "full PW M-matrix (jj coupling) with WP-basis kinematic prefactor",
            "dsigma_calibration_factor": float(best["dsigma_calibration_factor"]),
            "calibration_method": (
                "disabled" if args.no_dsigma_calibration else "log-space geometric mean fit per energy"
            ),
        },
        "units": {
            "dsigma_input_detected": exp_dsigma["source_unit"],
            "dsigma_output": exp_dsigma["output_unit"],
            "it11": "dimensionless",
            "T20": "dimensionless",
            "T21": "dimensionless",
            "T22": "dimensionless",
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
            "best_t20_csv": str(t20_csv),
            "best_t21_csv": str(t21_csv),
            "best_t22_csv": str(t22_csv),
            "best_dsigma_csv": str(dsigma_csv),
            "best_it11_svg": str(ay_svg),
            "best_t20_svg": str(t20_svg),
            "best_t21_svg": str(t21_svg),
            "best_t22_svg": str(t22_svg),
            "best_dsigma_svg": str(dsigma_svg),
        },
        "u_blocks_at_best_energy": _serialize_complex_block(blocks, int(best["q_idx"])),
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
