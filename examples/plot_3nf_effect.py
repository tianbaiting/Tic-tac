#!/usr/bin/env python3
"""
Purpose:
  Compare full partial-wave d+p observables across the three 3NF
  configurations produced by run_3nf_sweep.py, and overlay experimental data
  when available.

Data flow:
  <work-dir>/{no_3nf,zero_lec,witala}/solver_out/U_PW_elements_*.txt + q_kinematics
    -> parse into (J, pi) blocks
    -> assemble M(theta) via partial-wave summation (pw_amplitudes)
    -> dSigma/dOmega, iT11 from spin-1 trace formulas
    -> overlay plot (SVG) + per-run CSV + summary JSON
"""

from __future__ import annotations

import argparse
import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from compare_Ay_experiment import (
    read_experimental_dsigma,
    read_experimental_iT11,
)
from observable_units import UNIT_MB_PER_SR, convert_dsigma_value, normalize_dsigma_unit
from solver_u_file_utils import select_latest_u_file_family
from tlab_utils import format_tlab_label

from pw_amplitudes import (
    assemble_m_matrix,
    calibrate_dsigma_scale,
    list_solver_energies,
    observables_from_M,
    parse_q_kinematics,
    parse_solver_output,
)


RUN_ORDER = ["no_3nf", "zero_lec", "witala"]
RUN_COLORS = {
    "no_3nf":   "#6b7280",   # neutral grey
    "zero_lec": "#2563eb",   # blue
    "witala":   "#dc2626",   # red
}
RUN_LABELS = {
    "no_3nf":   "2NF only (no 3NF)",
    "zero_lec": "chiral 3NF, c_D=c_E=0",
    "witala":   "chiral 3NF, Witala c_D/c_E",
}


@dataclass
class RunObservables:
    name: str
    solver_out_dir: Path
    tlab_mev: float
    ecm_mev: float
    q_idx: int
    dsigma_calibration_factor: float
    dsigma_pred: List[float]
    it11_pred: List[float]


def _load_run(
    *,
    work_dir: Path,
    run_name: str,
    target_tlab_mev: float,
    exp_it11: Dict[str, List[float]],
    exp_dsigma: Dict[str, object],
    dsigma_output_unit: str,
) -> Optional[RunObservables]:
    run_dir = work_dir / run_name / "solver_out"
    if not run_dir.exists():
        print(f"[plot_3nf_effect] skip {run_name}: {run_dir} missing")
        return None

    u_files = select_latest_u_file_family(run_dir)
    if not u_files:
        print(f"[plot_3nf_effect] skip {run_name}: no U_PW_elements_*.txt in {run_dir}")
        return None

    blocks = parse_solver_output(run_dir, u_files=u_files)
    if not blocks:
        print(f"[plot_3nf_effect] skip {run_name}: no (J, pi) blocks parsed")
        return None
    try:
        q_grid = parse_q_kinematics(run_dir)
    except FileNotFoundError as exc:
        print(f"[plot_3nf_effect] skip {run_name}: {exc}")
        return None

    common_energies = list_solver_energies(blocks)
    if not common_energies:
        print(f"[plot_3nf_effect] skip {run_name}: no common solver energy")
        return None
    tlab_solv, ecm_solv = min(common_energies, key=lambda item: abs(item[0] - target_tlab_mev))
    qmap = {pt.q_idx: pt for pt in blocks[0].points}
    q_idx = next(
        (idx for idx, pt in qmap.items()
         if abs(pt.tlab - tlab_solv) < 1e-6 and abs(pt.ecm - ecm_solv) < 1e-6),
        None,
    )
    if q_idx is None or q_idx not in q_grid:
        print(f"[plot_3nf_effect] skip {run_name}: q_idx lookup failed")
        return None
    bin_info = q_grid[q_idx]

    # Raw dSigma/dOmega in mb/sr on experimental angle grid (no calibration yet).
    raw_dsigma_mb: List[float] = []
    raw_it11: List[float] = []
    for theta in exp_dsigma["angles"]:  # type: ignore[arg-type]
        M = assemble_m_matrix(blocks, q_idx, math.radians(theta), bin_info=bin_info)
        obs = observables_from_M(M)
        raw_dsigma_mb.append(convert_dsigma_value(obs.dsigma_fm2_per_sr, "fm2/sr", UNIT_MB_PER_SR))
    cal_factor = calibrate_dsigma_scale(
        raw_dsigma_mb,
        exp_dsigma["angles"],  # type: ignore[arg-type]
        [convert_dsigma_value(v, dsigma_output_unit, UNIT_MB_PER_SR) for v in exp_dsigma["values"]],  # type: ignore[arg-type]
        exp_dsigma["angles"],  # type: ignore[arg-type]
    )
    m_scale = math.sqrt(max(cal_factor, 1e-30))

    dsigma_pred: List[float] = []
    for theta in exp_dsigma["angles"]:  # type: ignore[arg-type]
        M = assemble_m_matrix(blocks, q_idx, math.radians(theta), bin_info=bin_info, extra_scale=m_scale)
        obs = observables_from_M(M)
        dsigma_pred.append(convert_dsigma_value(obs.dsigma_fm2_per_sr, "fm2/sr", dsigma_output_unit))

    it11_pred: List[float] = []
    for theta in exp_it11["angles"]:
        M = assemble_m_matrix(blocks, q_idx, math.radians(theta), bin_info=bin_info, extra_scale=m_scale)
        it11_pred.append(observables_from_M(M).iT11)

    return RunObservables(
        name=run_name,
        solver_out_dir=run_dir,
        tlab_mev=float(tlab_solv),
        ecm_mev=float(ecm_solv),
        q_idx=int(q_idx),
        dsigma_calibration_factor=float(cal_factor),
        dsigma_pred=dsigma_pred,
        it11_pred=it11_pred,
    )


def _write_overlay_svg(
    path: Path,
    *,
    title: str,
    x_label: str,
    y_label: str,
    angles: Sequence[float],
    runs: List[Tuple[str, Sequence[float]]],
    exp_vals: Optional[Sequence[float]] = None,
    use_log_y: bool = False,
) -> None:
    if not runs:
        raise ValueError("need at least one run to plot")
    width = 1020
    height = 660
    left = 100
    right = 230
    top = 70
    bottom = 100
    plot_w = width - left - right
    plot_h = height - top - bottom

    x_min = min(angles)
    x_max = max(angles)
    if x_max <= x_min:
        x_max = x_min + 1.0

    def _maybe_log(values: Sequence[float]) -> List[float]:
        if use_log_y:
            return [math.log10(max(v, 1e-12)) for v in values]
        return list(values)

    exp_y = _maybe_log(exp_vals) if exp_vals else None
    run_ys: List[Tuple[str, List[float]]] = [(name, _maybe_log(vals)) for name, vals in runs]

    all_y: List[float] = []
    for _, ys in run_ys:
        all_y.extend(ys)
    if exp_y is not None:
        all_y.extend(exp_y)
    y_min = min(all_y)
    y_max = max(all_y)
    if y_max <= y_min:
        y_max = y_min + 1.0
    y_pad = 0.08 * (y_max - y_min)
    y_min -= y_pad
    y_max += y_pad

    def x_to_px(x: float) -> float:
        return left + (x - x_min) * plot_w / (x_max - x_min)

    def y_to_px(y: float) -> float:
        return top + (y_max - y) * plot_h / (y_max - y_min)

    svg: List[str] = []
    svg.append(f"<svg xmlns='http://www.w3.org/2000/svg' width='{width}' height='{height}' viewBox='0 0 {width} {height}'>")
    svg.append("<rect x='0' y='0' width='100%' height='100%' fill='white'/>")
    svg.append(
        f"<text x='{width/2:.1f}' y='36' text-anchor='middle' "
        "font-family='DejaVu Sans, Arial, sans-serif' font-size='22' font-weight='700' fill='#111827'>"
        f"{title}</text>"
    )

    # gridlines
    for i in range(7):
        frac = i / 6.0
        y_val = y_min + frac * (y_max - y_min)
        y_px = y_to_px(y_val)
        svg.append(
            f"<line x1='{left}' y1='{y_px:.3f}' x2='{left + plot_w}' y2='{y_px:.3f}' "
            "stroke='#e5e7eb' stroke-width='1'/>"
        )
        label = f"{(10 ** y_val):.3g}" if use_log_y else f"{y_val:.3f}"
        svg.append(
            f"<text x='{left - 12}' y='{y_px + 4:.3f}' text-anchor='end' "
            "font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>"
            f"{label}</text>"
        )
    for i in range(9):
        frac = i / 8.0
        x_val = x_min + frac * (x_max - x_min)
        x_px = x_to_px(x_val)
        svg.append(
            f"<line x1='{x_px:.3f}' y1='{top}' x2='{x_px:.3f}' y2='{top + plot_h}' "
            "stroke='#f3f4f6' stroke-width='1'/>"
        )
        svg.append(
            f"<text x='{x_px:.3f}' y='{top + plot_h + 24}' text-anchor='middle' "
            "font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>"
            f"{x_val:.1f}</text>"
        )

    svg.append(
        f"<rect x='{left}' y='{top}' width='{plot_w}' height='{plot_h}' "
        "fill='none' stroke='#111827' stroke-width='1.2'/>"
    )

    # runs as lines
    for name, ys in run_ys:
        pts = [f"{x_to_px(x):.3f},{y_to_px(y):.3f}" for x, y in zip(angles, ys)]
        color = RUN_COLORS.get(name, "#000000")
        svg.append(
            f"<polyline fill='none' stroke='{color}' stroke-width='2.4' "
            f"stroke-linejoin='round' points='{' '.join(pts)}'/>"
        )

    # experimental markers
    if exp_y is not None:
        for x, y in zip(angles, exp_y):
            svg.append(
                f"<circle cx='{x_to_px(x):.3f}' cy='{y_to_px(y):.3f}' r='3.2' "
                "fill='#111827' stroke='white' stroke-width='0.8'/>"
            )

    # axis labels
    svg.append(
        f"<text x='{left + plot_w/2:.1f}' y='{height - 36}' text-anchor='middle' "
        "font-family='DejaVu Sans, Arial, sans-serif' font-size='16' fill='#111827'>"
        f"{x_label}</text>"
    )
    y_axis_text = y_label + (" (log10 scale)" if use_log_y else "")
    svg.append(
        f"<text x='32' y='{top + plot_h/2:.1f}' text-anchor='middle' "
        f"transform='rotate(-90 32 {top + plot_h/2:.1f})' "
        "font-family='DejaVu Sans, Arial, sans-serif' font-size='16' fill='#111827'>"
        f"{y_axis_text}</text>"
    )

    # legend
    legend_x = left + plot_w + 20
    legend_y = top + 10
    for i, (name, _) in enumerate(run_ys):
        color = RUN_COLORS.get(name, "#000000")
        label = RUN_LABELS.get(name, name)
        y_row = legend_y + 28 * i
        svg.append(
            f"<line x1='{legend_x}' y1='{y_row:.1f}' x2='{legend_x + 30}' y2='{y_row:.1f}' "
            f"stroke='{color}' stroke-width='2.4'/>"
        )
        svg.append(
            f"<text x='{legend_x + 38}' y='{y_row + 4:.1f}' "
            "font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>"
            f"{label}</text>"
        )
    if exp_y is not None:
        y_row = legend_y + 28 * len(run_ys)
        svg.append(
            f"<circle cx='{legend_x + 15}' cy='{y_row:.1f}' r='3.6' "
            "fill='#111827' stroke='white' stroke-width='0.8'/>"
        )
        svg.append(
            f"<text x='{legend_x + 38}' y='{y_row + 4:.1f}' "
            "font-family='DejaVu Sans, Arial, sans-serif' font-size='13' fill='#111827'>experiment</text>"
        )

    svg.append("</svg>")
    path.write_text("\n".join(svg) + "\n", encoding="utf-8")


def _write_overlay_csv(path: Path, angles: Sequence[float], runs: List[Tuple[str, Sequence[float]]], exp: Optional[Sequence[float]] = None) -> None:
    header = ["theta_cm_deg"] + [name for name, _ in runs]
    if exp is not None:
        header.append("experiment")
    lines = [",".join(header)]
    for i, th in enumerate(angles):
        row = [f"{th:.6f}"]
        for _, ys in runs:
            row.append(f"{ys[i]:.12e}")
        if exp is not None:
            row.append(f"{exp[i]:.12e}")
        lines.append(",".join(row))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Overlay 3NF-sweep observables against experiment")
    p.add_argument("--work-dir", default="output/3nf_sweep", help="directory produced by run_3nf_sweep.py")
    p.add_argument("--target-tlab-mev", type=float, default=64.5)
    p.add_argument("--dsigma-data", default="data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt")
    p.add_argument("--tensor-data", default="data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt")
    p.add_argument("--dsigma-unit", default=UNIT_MB_PER_SR)
    p.add_argument("--no-exp", action="store_true",
                   help="Skip overlaying experimental points (use when energy does not match data)")
    return p


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    work_dir = (root / args.work_dir).resolve()
    if not work_dir.exists():
        raise SystemExit(f"work_dir does not exist: {work_dir}")

    dsigma_output_unit = normalize_dsigma_unit(args.dsigma_unit)

    exp_it11 = read_experimental_iT11((root / args.tensor_data).resolve())
    exp_dsigma = read_experimental_dsigma(
        (root / args.dsigma_data).resolve(),
        target_unit=dsigma_output_unit,
    )

    runs: List[RunObservables] = []
    for name in RUN_ORDER:
        run = _load_run(
            work_dir=work_dir,
            run_name=name,
            target_tlab_mev=args.target_tlab_mev,
            exp_it11=exp_it11,
            exp_dsigma=exp_dsigma,
            dsigma_output_unit=dsigma_output_unit,
        )
        if run is not None:
            runs.append(run)

    if not runs:
        raise SystemExit(f"no usable solver outputs under {work_dir}/{{no_3nf,zero_lec,witala}}/solver_out/")

    tlab_label = format_tlab_label(args.target_tlab_mev)

    # dSigma overlay
    dsigma_runs = [(r.name, r.dsigma_pred) for r in runs]
    dsigma_exp = None if args.no_exp else list(exp_dsigma["values"])
    dsigma_angles = list(exp_dsigma["angles"])

    dsigma_svg = work_dir / f"3nf_effect_dsigma_{tlab_label}.svg"
    dsigma_csv = work_dir / f"3nf_effect_dsigma_{tlab_label}.csv"
    _write_overlay_svg(
        dsigma_svg,
        title=f"dσ/dΩ: 3NF effect (target Tlab = {args.target_tlab_mev:.2f} MeV)",
        x_label="θ_cm (deg)",
        y_label=f"dσ/dΩ [{dsigma_output_unit}]",
        angles=dsigma_angles,
        runs=dsigma_runs,
        exp_vals=dsigma_exp,
        use_log_y=True,
    )
    _write_overlay_csv(dsigma_csv, dsigma_angles, dsigma_runs, dsigma_exp)

    # iT11 overlay
    it11_runs = [(r.name, r.it11_pred) for r in runs]
    it11_exp = None if args.no_exp else list(exp_it11["values"])
    it11_angles = list(exp_it11["angles"])

    it11_svg = work_dir / f"3nf_effect_iT11_{tlab_label}.svg"
    it11_csv = work_dir / f"3nf_effect_iT11_{tlab_label}.csv"
    _write_overlay_svg(
        it11_svg,
        title=f"iT11: 3NF effect (target Tlab = {args.target_tlab_mev:.2f} MeV)",
        x_label="θ_cm (deg)",
        y_label="iT11",
        angles=it11_angles,
        runs=it11_runs,
        exp_vals=it11_exp,
        use_log_y=False,
    )
    _write_overlay_csv(it11_csv, it11_angles, it11_runs, it11_exp)

    inv_summary: Dict[str, object] = {
        "target_tlab_mev": args.target_tlab_mev,
        "runs": [
            {
                "name": r.name,
                "solver_tlab_mev": r.tlab_mev,
                "solver_ecm_mev": r.ecm_mev,
                "q_idx": r.q_idx,
                "dsigma_calibration_factor_mb_per_sr": r.dsigma_calibration_factor,
            }
            for r in runs
        ],
        "outputs": {
            "dsigma_svg": str(dsigma_svg),
            "dsigma_csv": str(dsigma_csv),
            "it11_svg": str(it11_svg),
            "it11_csv": str(it11_csv),
        },
    }
    summary_file = work_dir / f"3nf_effect_summary_{tlab_label}.json"
    summary_file.write_text(json.dumps(inv_summary, indent=2), encoding="utf-8")

    print("=" * 70)
    print(f"3NF effect at target Tlab = {args.target_tlab_mev:.3f} MeV")
    print("=" * 70)
    for r in runs:
        print(f"  {r.name:<10}  Tlab_eff={r.tlab_mev:>7.2f} MeV  "
              f"q_idx={r.q_idx}  dsigma_cal={r.dsigma_calibration_factor:.4e}")
    print()
    print(f"dSigma overlay:   {dsigma_svg}")
    print(f"iT11 overlay:     {it11_svg}")
    print(f"summary JSON:     {summary_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
