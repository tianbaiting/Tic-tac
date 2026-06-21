#!/usr/bin/env python3
"""Quantify how accurately Tic-tac reproduces the Miller Gate 2 nd Ay(n)
benchmark, with vs without 3NF.

For each solver output directory given on the command line, compute Ay(n) at
Tlab = 10/35/67 MeV on the experimental theta grid, then report MAE / RMSE /
rel-RMSE / achieved-fraction against experiment.

Reference data:
    data/miller_paper3_fig6/Ay_n_expt_Elab{10,35,67}MeV.csv
    data/miller_paper3_fig6/Ay_n_N3LO_Elab{10,35,67}MeV.csv

Example:
    python3 examples/compare_3nf_accuracy.py \\
        --label "2NF J3N<=9/2 (v4)"  --dir CPP/Output/miller_gate2_v4_J3N9_J2N3 \\
        --label "2NF J3N<=13/2 (v5)" --dir CPP/Output/miller_gate2_v5_J3N13_J2N3 \\
        --label "3NF chiral N2LO J3N<=9/2 (labenpg)" --dir CPP/Output/labenpg_3NF_J3N9 \\
        --json-output output/miller_gate2_3nf_accuracy.json
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "examples"))

import pw_amplitudes as pw  # noqa: E402

DATA = REPO / "data" / "miller_paper3_fig6"

ENERGIES: List[Tuple[int, str]] = [
    (10, "Elab10MeV"),
    (35, "Elab35MeV"),
    (67, "Elab67MeV"),
]


def load_reference(tag: str) -> Tuple[List[Tuple[float, float, float]],
                                      List[Tuple[float, float]]]:
    """Return (experimental [(theta, Ay, err)], N3LO [(theta, Ay)])."""
    expt: List[Tuple[float, float, float]] = []
    with (DATA / f"Ay_n_expt_{tag}.csv").open() as f:
        for row in csv.reader(f):
            if not row or row[0].startswith("#") or row[0] == "theta_cm_deg":
                continue
            theta = float(row[0])
            value = float(row[1])
            err = float(row[2]) if len(row) > 2 else 0.0
            expt.append((theta, value, err))
    n3lo: List[Tuple[float, float]] = []
    with (DATA / f"Ay_n_N3LO_{tag}.csv").open() as f:
        for row in csv.reader(f):
            if not row or row[0].startswith("#") or row[0] == "theta_cm_deg":
                continue
            n3lo.append((float(row[0]), float(row[1])))
    return expt, n3lo


def interp_linear(grid_x: List[float], grid_y: List[float],
                  x_q: float) -> float:
    """Linear interpolation; returns 0.0 outside [grid_x[0], grid_x[-1]]."""
    if x_q <= grid_x[0]:
        return grid_y[0] if x_q == grid_x[0] else 0.0
    if x_q >= grid_x[-1]:
        return grid_y[-1] if x_q == grid_x[-1] else 0.0
    for i in range(len(grid_x) - 1):
        if grid_x[i] <= x_q <= grid_x[i + 1]:
            x0, x1 = grid_x[i], grid_x[i + 1]
            y0, y1 = grid_y[i], grid_y[i + 1]
            if x1 == x0:
                return y0
            return y0 + (y1 - y0) * (x_q - x0) / (x1 - x0)
    return 0.0


def evaluate_run(run_dir: Path) -> Dict:
    """Compute Ay(n) curves + accuracy metrics for one solver output dir."""
    blocks = pw.parse_solver_output(run_dir)
    q_kin = pw.parse_q_kinematics(run_dir)
    points = sorted({(p.q_idx, p.tlab) for b in blocks for p in b.points})

    theta_grid = [float(t) for t in range(0, 181, 1)]
    result: Dict = {"run_dir": str(run_dir), "energies": {}}

    for target_mev, tag in ENERGIES:
        expt, n3lo = load_reference(tag)
        q_idx, tlab = min(points, key=lambda qt: abs(qt[1] - target_mev))
        obs = pw.observables_at_angles(
            blocks, q_idx, theta_grid, bin_info=q_kin[q_idx])
        ay_curve = [o.Ay_n for o in obs]

        # Match model to each experimental theta
        preds_at_expt: List[float] = []
        exp_values: List[float] = []
        exp_errs: List[float] = []
        exp_thetas: List[float] = []
        for theta, value, err in expt:
            preds_at_expt.append(
                interp_linear(theta_grid, ay_curve, theta))
            exp_values.append(value)
            exp_errs.append(err)
            exp_thetas.append(theta)

        # Metrics
        diffs = [p - e for p, e in zip(preds_at_expt, exp_values)]
        mae = float(np.mean(np.abs(diffs))) if diffs else float("nan")
        rmse = float(math.sqrt(np.mean(np.square(diffs)))) if diffs else float("nan")
        max_abs = float(max(np.abs(diffs))) if diffs else float("nan")
        exp_rms = float(math.sqrt(np.mean(np.square(exp_values)))) if exp_values else 1.0
        rel_rmse = rmse / exp_rms if exp_rms > 0 else float("nan")
        # Mean absolute experimental |Ay| for "fraction of signal captured"
        exp_abs_mean = float(np.mean(np.abs(exp_values))) if exp_values else 1.0
        # "Achieved fraction" = 1 - MAE / exp_abs_mean (clipped to [0, inf))
        achieved = max(0.0, 1.0 - mae / exp_abs_mean) if exp_abs_mean > 0 else 0.0

        # Extremum of predicted curve (angle and value)
        ext_idx = int(np.argmax(np.abs(ay_curve)))
        ext_theta = theta_grid[ext_idx]
        ext_val = ay_curve[ext_idx]

        # N3LO extremum for context
        if n3lo:
            n3_th = [t for t, _ in n3lo]
            n3_v = [v for _, v in n3lo]
            n3_ext_idx = int(np.argmax(np.abs(n3_v)))
            n3_ext = (n3_th[n3_ext_idx], n3_v[n3_ext_idx])
        else:
            n3_ext = None

        result["energies"][target_mev] = {
            "solver_tlab_mev": float(tlab),
            "q_idx": int(q_idx),
            "theta_grid": theta_grid,
            "ay_curve": ay_curve,
            "exp_thetas": exp_thetas,
            "exp_values": exp_values,
            "exp_errs": exp_errs,
            "preds_at_expt": preds_at_expt,
            "metrics": {
                "mae": mae,
                "rmse": rmse,
                "max_abs_error": max_abs,
                "rel_rmse": rel_rmse,
                "achieved_fraction": achieved,
                "exp_abs_signal_mean": exp_abs_mean,
            },
            "pred_extremum": {"theta_deg": ext_theta, "value": ext_val},
            "n3lo_extremum": (
                {"theta_deg": n3_ext[0], "value": n3_ext[1]}
                if n3_ext else None
            ),
        }
    return result


def print_table(results: List[Tuple[str, Dict]]) -> None:
    """Pretty-print accuracy comparison table."""
    print()
    print("=" * 100)
    print("Miller Gate 2  nd  Ay(n)  accuracy  (Tlab = 10 / 35 / 67 MeV)")
    print("Reference: Miller PRC 107, 014002 (2023) Fig. 6  +  EXFOR experimental markers")
    print("=" * 100)
    for label, res in results:
        print(f"\n--- {label} ---")
        print(f"    dir: {res['run_dir']}")
        hdr = (f"    {'Tlab':>5} {'solverTlab':>11} {'MAE':>8} {'RMSE':>8} "
               f"{'relRMSE':>9} {'maxErr':>8} {'achFrac':>8}  "
               f"{'predExt':>14} {'expt|signal|':>12}")
        print(hdr)
        for target_mev, _ in ENERGIES:
            e = res["energies"][target_mev]
            m = e["metrics"]
            ext = e["pred_extremum"]
            print(
                f"    {target_mev:>5} {e['solver_tlab_mev']:>11.2f} "
                f"{m['mae']:>8.4f} {m['rmse']:>8.4f} "
                f"{m['rel_rmse']:>9.3f} {m['max_abs_error']:>8.4f} "
                f"{m['achieved_fraction']:>8.3f}  "
                f"({ext['theta_deg']:>4.0f}, {ext['value']:+.3f})   "
                f"{m['exp_abs_signal_mean']:>10.3f}"
            )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", action="append", default=[],
                    help="Display label for the following --dir (repeatable).")
    ap.add_argument("--dir", action="append", default=[],
                    help="Solver output directory (repeatable).")
    ap.add_argument("--json-output", default=None,
                    help="Write full metrics to this JSON path.")
    args = ap.parse_args()

    if len(args.label) != len(args.dir):
        ap.error("--label and --dir must be given the same number of times "
                 f"(got {len(args.label)} labels, {len(args.dir)} dirs).")

    results: List[Tuple[str, Dict]] = []
    for label, d in zip(args.label, args.dir):
        run_dir = Path(d)
        if not run_dir.is_absolute():
            run_dir = REPO / run_dir
        if not run_dir.exists():
            print(f"WARN: {run_dir} does not exist, skipping.", file=sys.stderr)
            continue
        print(f"Evaluating: {label} -> {run_dir}", file=sys.stderr)
        results.append((label, evaluate_run(run_dir)))

    if not results:
        print("No runs to evaluate.", file=sys.stderr)
        return 1

    print_table(results)

    if args.json_output:
        out_path = Path(args.json_output)
        if not out_path.is_absolute():
            out_path = REPO / out_path
        out_path.parent.mkdir(parents=True, exist_ok=True)
        # Strip the heavy per-theta curve arrays for the JSON dump; keep metrics.
        slim = []
        for label, res in results:
            slim_entry = {"label": label, "run_dir": res["run_dir"],
                          "energies": {}}
            for tlab, e in res["energies"].items():
                slim_entry["energies"][str(tlab)] = {
                    "solver_tlab_mev": e["solver_tlab_mev"],
                    "q_idx": e["q_idx"],
                    "metrics": e["metrics"],
                    "pred_extremum": e["pred_extremum"],
                    "n3lo_extremum": e["n3lo_extremum"],
                    "exp_thetas": e["exp_thetas"],
                    "exp_values": e["exp_values"],
                    "exp_errs": e["exp_errs"],
                    "preds_at_expt": e["preds_at_expt"],
                }
            slim.append(slim_entry)
        out_path.write_text(json.dumps(slim, indent=2))
        print(f"\nJSON metrics written to {out_path}", file=sys.stderr)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
