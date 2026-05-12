#!/usr/bin/env python3
"""Plot Miller Gate 2 PW convergence + Ay_n comparison.

Two outputs:
  - miller_gate2_angular_panels.{svg,png}: three side-by-side panels for
    Tlab = 10/35/67 MeV showing v5 (J_3N<=13/2) vs Miller N3LO digitized
    curve vs experimental markers.
  - miller_gate2_pw_convergence.{svg,png}: Ay_n at theta=90 deg as a function
    of two_J_3N_max truncation, one line per energy.

Reads v1..v5 from CPP/Output/miller_gate2_*; Miller reference data from
data/miller_paper3_fig6/Ay_n_{expt,N3LO}_Elab{10,35,67}MeV.csv.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "examples"))

import pw_amplitudes as pw  # noqa: E402

DATA = REPO / "data" / "miller_paper3_fig6"
OUTPUT_DIR = REPO / "output" / "miller_gate2_plots"

# Solver runs in order of increasing truncation.
RUNS: List[Tuple[str, int, Path]] = [
    ("v1: J_3N<=1/2",  1,  REPO / "CPP" / "Output" / "miller_gate2_np20_2nf"),
    ("v3: J_3N<=3/2",  3,  REPO / "CPP" / "Output" / "miller_gate2_v3_J3N3_J2N3"),
    ("v2: J_3N<=5/2",  5,  REPO / "CPP" / "Output" / "miller_gate2_v2_J3N5_J2N3"),
    ("v4: J_3N<=9/2",  9,  REPO / "CPP" / "Output" / "miller_gate2_v4_J3N9_J2N3"),
    ("v5: J_3N<=13/2", 13, REPO / "CPP" / "Output" / "miller_gate2_v5_J3N13_J2N3"),
]

ENERGIES: List[Tuple[int, str]] = [
    (10, "Elab10MeV"),
    (35, "Elab35MeV"),
    (67, "Elab67MeV"),
]


def load_run(run_dir: Path):
    blocks = pw.parse_solver_output(run_dir)
    q_kin = pw.parse_q_kinematics(run_dir)
    points = sorted({(p.q_idx, p.tlab) for b in blocks for p in b.points})
    return blocks, q_kin, points


def load_reference(tag: str) -> Tuple[List[Tuple[float, float, float]], List[Tuple[float, float]]]:
    """Return (experimental [(theta, Ay, err)], N3LO [(theta, Ay)]) for this energy tag."""
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


def compute_angular_curve(blocks, q_kin, points, target_mev: int,
                          theta_grid: List[float]) -> Tuple[float, np.ndarray]:
    q_idx, tlab = min(points, key=lambda qt: abs(qt[1] - target_mev))
    obs = pw.observables_at_angles(blocks, q_idx, theta_grid, bin_info=q_kin[q_idx])
    return tlab, np.array([o.Ay_n for o in obs])


def plot_angular_panels(output_path: Path) -> None:
    """Three side-by-side panels: Ay(n) vs theta_cm at 10/35/67 MeV."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=False)
    theta_grid = list(np.arange(0, 181, 1))

    # Load v5 (final) once; we'll overlay v1, v2, v4 for context.
    runs_to_plot = [(label, jmax, *load_run(rd)) for label, jmax, rd in RUNS]
    colors_legacy = ["#cfd8dc", "#90a4ae", "#607d8b", "#37474f"]  # gradient grays for v1..v4
    color_final = "#1565c0"

    for ax, (target_mev, tag) in zip(axes, ENERGIES):
        expt, n3lo = load_reference(tag)

        # Plot v1, v3, v2, v4 as faint convergence overlay (last entry plotted bold).
        for run_idx, (label, jmax, blocks, q_kin, points) in enumerate(runs_to_plot[:-1]):
            tlab, ay = compute_angular_curve(blocks, q_kin, points, target_mev, theta_grid)
            ax.plot(theta_grid, ay, color=colors_legacy[run_idx],
                    linewidth=0.9, alpha=0.7, label=label)

        # v5 — bold
        label, jmax, blocks, q_kin, points = runs_to_plot[-1]
        tlab_v5, ay_v5 = compute_angular_curve(blocks, q_kin, points, target_mev, theta_grid)
        ax.plot(theta_grid, ay_v5, color=color_final, linewidth=2.2,
                label=f"{label} (this work)")

        # Miller N3LO digitized
        if n3lo:
            n3_th = [t for t, _ in n3lo]
            n3_v  = [v for _, v in n3lo]
            ax.plot(n3_th, n3_v, color="#c62828", linewidth=1.6, linestyle="--",
                    label="Miller N3LO 2NF (Fig. 6)")

        # Experimental markers
        if expt:
            xs = [t for t, _, _ in expt]
            ys = [v for _, v, _ in expt]
            es = [e for _, _, e in expt]
            ax.errorbar(xs, ys, yerr=es, fmt="o", color="black", markersize=5,
                        capsize=3, label="Expt (Miller Fig. 6)")

        ax.axhline(0.0, color="gray", linewidth=0.5)
        ax.set_xlabel(r"$\theta_{c.m.}$ (deg)")
        ax.set_xlim(0, 180)
        ax.set_xticks(range(0, 181, 30))
        ax.grid(alpha=0.3)
        ax.set_title(rf"$T_{{lab}}\approx {target_mev}$ MeV   (solver $T_{{lab}}={tlab_v5:.2f}$ MeV)")

    axes[0].set_ylabel(r"$A_y(n)$ (Madison)")
    axes[0].legend(loc="lower left", fontsize=8, framealpha=0.95)

    fig.suptitle(r"Miller Gate 2 — $nd$ elastic $A_y(n)$ at $N_p=20$, $J_{2N}\leq 3$, Nijmegen-I 2NF",
                 y=1.02, fontsize=12)
    fig.tight_layout()
    fig.savefig(output_path.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)


def plot_pw_convergence(output_path: Path) -> None:
    """Ay_n at theta=90 deg vs two_J_3N_max truncation."""
    fig, ax = plt.subplots(figsize=(7.5, 5))

    runs_to_plot = [(label, jmax, *load_run(rd)) for label, jmax, rd in RUNS]

    energy_colors = {10: "#1b5e20", 35: "#e65100", 67: "#1565c0"}
    expt90 = {}
    n3lo90 = {}

    for target_mev, tag in ENERGIES:
        expt, n3lo = load_reference(tag)
        e90 = next((v for t, v, _ in expt if abs(t - 90) < 1), None)
        n90 = next((v for t, v in n3lo if abs(t - 90) < 3), None)
        expt90[target_mev] = e90
        n3lo90[target_mev] = n90

        # Trajectory across truncations
        traj_x = []
        traj_y = []
        for label, jmax, blocks, q_kin, points in runs_to_plot:
            _, ay = compute_angular_curve(blocks, q_kin, points, target_mev, [90.0])
            traj_x.append(jmax / 2.0)  # convert two_J_3N_max -> J_3N_max
            traj_y.append(float(ay[0]))
        c = energy_colors[target_mev]
        ax.plot(traj_x, traj_y, marker="o", color=c, linewidth=2,
                markersize=7, label=rf"$T_{{lab}}={target_mev}$ MeV (this work)")
        if e90 is not None:
            ax.axhline(e90, color=c, linewidth=1.0, linestyle=":", alpha=0.6)
            ax.text(traj_x[-1] + 0.15, e90, f"expt {target_mev}: {e90:+.3f}",
                    fontsize=8, color=c, va="center")
        if n90 is not None:
            ax.axhline(n90, color=c, linewidth=1.0, linestyle="--", alpha=0.5)

    ax.set_xlabel(r"$J_{3N,\max}$ truncation")
    ax.set_ylabel(r"$A_y(n)$ at $\theta_{c.m.}=90^\circ$")
    ax.set_title("Partial-wave convergence of $A_y(n)$ (Np=20, 2NF Nijmegen-I)")
    ax.set_xticks([0.5, 1.5, 2.5, 4.5, 6.5])
    ax.set_xticklabels(["1/2", "3/2", "5/2", "9/2", "13/2"])
    ax.axhline(0.0, color="gray", linewidth=0.5)
    ax.grid(alpha=0.3)
    ax.legend(loc="best", fontsize=9)

    # Caption: dashed lines = Miller N3LO at theta=90, dotted = expt
    info = ("dashed: Miller N3LO at $\\theta=90^\\circ$   "
            "dotted: experimental Ay(n)@90")
    ax.text(0.5, -0.13, info, transform=ax.transAxes, ha="center", fontsize=8,
            color="gray")

    fig.tight_layout()
    fig.savefig(output_path.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(output_path.with_suffix(".png"), bbox_inches="tight", dpi=150)
    plt.close(fig)


def main() -> int:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_angular_panels(OUTPUT_DIR / "miller_gate2_angular_panels")
    plot_pw_convergence(OUTPUT_DIR / "miller_gate2_pw_convergence")
    print(f"Wrote plots to {OUTPUT_DIR}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
