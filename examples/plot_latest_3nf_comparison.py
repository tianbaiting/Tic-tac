#!/usr/bin/env python3
"""Plot the latest computed Ay(n) data: experiment vs 2NF vs 2NF+3NF.

Reads the pre-computed comparison CSV produced by the solver pipeline
(together with its JSON summary for the energy bracket / RMSE) and overlays
the digitized experimental points.  Pure plotting -- no solver calls.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402


def load_csv_curve(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta, ay_2nf, ay_3nf = [], [], []
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle):
            if not row or row[0].startswith("#") or row[0] == "theta_cm_deg":
                continue
            theta.append(float(row[0]))
            ay_2nf.append(float(row[1]))
            ay_3nf.append(float(row[2]))
    return np.asarray(theta), np.asarray(ay_2nf), np.asarray(ay_3nf)


def load_experiment(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta, ay, err = [], [], []
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle):
            if not row or row[0].startswith("#") or row[0] == "theta_cm_deg":
                continue
            theta.append(float(row[0]))
            ay.append(float(row[1]))
            err.append(float(row[2]))
    return np.asarray(theta), np.asarray(ay), np.asarray(err)


def rmse(theta_grid: np.ndarray, curve: np.ndarray,
         exp_theta: np.ndarray, exp_ay: np.ndarray) -> float:
    pred = np.interp(exp_theta, theta_grid, curve)
    return float(math.sqrt(np.mean(np.square(pred - exp_ay))))


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    p = argparse.ArgumentParser()
    p.add_argument("--csv", type=Path,
                   default=repo / "output/ay_3nf_current/Ay_2nf_vs_3nf_10MeV.csv")
    p.add_argument("--json", type=Path,
                   default=repo / "output/ay_3nf_current/Ay_2nf_vs_3nf_10MeV.json")
    p.add_argument("--experiment", type=Path,
                   default=repo / "data/miller_paper3_fig6/Ay_n_expt_Elab10MeV.csv")
    p.add_argument("--output", type=Path,
                   default=repo / "output/ay_3nf_current/Ay_3nf_comparison_latest")
    args = p.parse_args()

    theta, ay_2nf, ay_3nf = load_csv_curve(args.csv)
    exp_theta, exp_ay, exp_err = load_experiment(args.experiment)

    summary = {}
    if args.json.exists():
        summary = json.loads(args.json.read_text(encoding="utf-8"))

    rmse_2nf = rmse(theta, ay_2nf, exp_theta, exp_ay)
    rmse_3nf = rmse(theta, ay_3nf, exp_theta, exp_ay)
    delta = ay_3nf - ay_2nf

    fig, (ax, dax) = plt.subplots(
        2, 1, figsize=(10.0, 8.0), sharex=True,
        gridspec_kw={"height_ratios": [3.0, 1.1], "hspace": 0.10},
    )

    ax.errorbar(
        exp_theta, exp_ay, yerr=exp_err, fmt="o", color="black",
        markersize=5.5, capsize=3, linewidth=1.1,
        label="Experiment (Bunker 1968, digitized)", zorder=4,
    )
    ax.plot(theta, ay_2nf, color="#2563eb", linewidth=2.4,
            label=rf"2NF only  (RMSE={rmse_2nf:.3f})")
    ax.plot(theta, ay_3nf, color="#dc2626", linewidth=2.4, linestyle="--",
            label=rf"2NF + approximate 3NF  (RMSE={rmse_3nf:.3f})")
    ax.axhline(0.0, color="#9ca3af", linewidth=0.8)
    ax.set_ylabel(r"$A_y(n)$")
    ax.set_xlim(0.0, 180.0)
    ax.grid(alpha=0.25)
    ax.legend(loc="upper left", framealpha=0.95)

    dax.plot(theta, delta, color="#7c3aed", linewidth=2.1)
    dax.fill_between(theta, 0.0, delta, color="#7c3aed", alpha=0.18)
    dax.axhline(0.0, color="black", linewidth=0.8)
    dax.set_xlabel(r"$\theta_{\mathrm{c.m.}}$ (deg)")
    dax.set_ylabel(r"$\Delta A_y$")
    dax.set_xlim(0.0, 180.0)
    dax.grid(alpha=0.25)

    interp = summary.get("energy_interpolation", {})
    lo = interp.get("lower_tlab_mev")
    hi = interp.get("upper_tlab_mev")
    bracket = ""
    if lo is not None and hi is not None:
        bracket = rf"  ($T_{{lab}}$: {lo:.2f}$\leftrightarrow${hi:.2f} MeV)"
    tlab = summary.get("target_tlab_mev", 10.0)
    fig.suptitle(
        rf"$nd$ $A_y(n)$ at $T_{{lab}}={tlab:.1f}$ MeV: "
        rf"experiment vs. 2NF / 2NF+3NF" + bracket,
        fontsize=14,
    )
    fig.text(
        0.5, 0.012,
        "Experimental points digitized from Miller & Ekström, PRC 107, 014002 (2023), Fig. 6.",
        ha="center", fontsize=8.5, color="#4b5563",
    )
    fig.subplots_adjust(top=0.90, bottom=0.11)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output.with_suffix(".png"), dpi=190, bbox_inches="tight")
    fig.savefig(args.output.with_suffix(".svg"), bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {args.output.with_suffix('.png')}")
    print(f"RMSE 2NF={rmse_2nf:.4f}  2NF+3NF={rmse_3nf:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
