#!/usr/bin/env python3
"""Plot a diagnostic Ay(n) comparison for paired 2NF and 2NF+3NF runs.

The script interpolates observables between two solver energies bracketing the
requested laboratory energy.  It is intended for transparent diagnostics, not
for hiding an unconverged momentum or partial-wave grid.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "examples"))

import pw_amplitudes as pw  # noqa: E402


@dataclass
class CurveResult:
    curve: np.ndarray
    lower_tlab: float
    upper_tlab: float
    interpolation_weight: float
    convergence_counts: Dict[str, int]


def load_experiment(path: Path) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    theta: List[float] = []
    ay: List[float] = []
    err: List[float] = []
    with path.open(encoding="utf-8") as handle:
        for row in csv.reader(handle):
            if not row or row[0].startswith("#") or row[0] == "theta_cm_deg":
                continue
            theta.append(float(row[0]))
            ay.append(float(row[1]))
            err.append(float(row[2]))
    if not theta:
        raise ValueError(f"no experimental data in {path}")
    return np.asarray(theta), np.asarray(ay), np.asarray(err)


def convergence_counts(run_dir: Path) -> Dict[str, int]:
    counts = {"converged": 0, "truncated": 0, "unknown": 0}
    for path in sorted(run_dir.glob("U_PW_convergence_*.txt")):
        for raw in path.read_text(encoding="utf-8").splitlines():
            fields = raw.split()
            if not fields or fields[0].startswith("#") or len(fields) < 5:
                continue
            status = int(fields[3])
            if status == 1:
                counts["converged"] += 1
            elif status == 2:
                counts["truncated"] += 1
            else:
                counts["unknown"] += 1
    return counts


def energy_to_q_idx(blocks: Sequence[pw.JPiBlock], tlab: float, ecm: float) -> int:
    candidates = blocks[0].points
    match = min(candidates, key=lambda point: abs(point.tlab - tlab) + abs(point.ecm - ecm))
    if abs(match.tlab - tlab) > 1e-6 or abs(match.ecm - ecm) > 1e-6:
        raise ValueError(f"could not map solver energy Tlab={tlab} Ecm={ecm} to q index")
    return int(match.q_idx)


def interpolate_curve(run_dir: Path, target_tlab: float,
                      angles: Sequence[float]) -> CurveResult:
    blocks = pw.parse_solver_output(run_dir)
    q_kinematics = pw.parse_q_kinematics(run_dir)
    energies = sorted(pw.list_solver_energies(blocks))
    lower = [entry for entry in energies if entry[0] <= target_tlab]
    upper = [entry for entry in energies if entry[0] >= target_tlab]
    if not lower or not upper:
        available = ", ".join(f"{entry[0]:.3f}" for entry in energies)
        raise ValueError(
            f"target Tlab={target_tlab} is not bracketed; available energies: {available}"
        )
    lo = lower[-1]
    hi = upper[0]
    q_lo = energy_to_q_idx(blocks, *lo)
    q_hi = energy_to_q_idx(blocks, *hi)

    obs_lo = pw.observables_at_angles(
        blocks, q_lo, angles, bin_info=q_kinematics[q_lo]
    )
    curve_lo = np.asarray([obs.Ay_n for obs in obs_lo])
    if q_lo == q_hi:
        weight = 0.0
        curve = curve_lo
    else:
        obs_hi = pw.observables_at_angles(
            blocks, q_hi, angles, bin_info=q_kinematics[q_hi]
        )
        curve_hi = np.asarray([obs.Ay_n for obs in obs_hi])
        weight = (target_tlab - lo[0]) / (hi[0] - lo[0])
        curve = (1.0 - weight) * curve_lo + weight * curve_hi

    return CurveResult(
        curve=curve,
        lower_tlab=float(lo[0]),
        upper_tlab=float(hi[0]),
        interpolation_weight=float(weight),
        convergence_counts=convergence_counts(run_dir),
    )


def rmse_at_experiment(theta_grid: np.ndarray, curve: np.ndarray,
                       exp_theta: np.ndarray, exp_ay: np.ndarray) -> float:
    prediction = np.interp(exp_theta, theta_grid, curve)
    return float(math.sqrt(np.mean(np.square(prediction - exp_ay))))


def write_csv(path: Path, theta: np.ndarray, ay_2nf: np.ndarray,
              ay_3nf: np.ndarray) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["theta_cm_deg", "Ay_2NF", "Ay_2NF_plus_approx_3NF", "delta_Ay"])
        for values in zip(theta, ay_2nf, ay_3nf):
            angle, no_3nf, with_3nf = values
            writer.writerow([
                f"{angle:.6f}", f"{no_3nf:.12e}", f"{with_3nf:.12e}",
                f"{with_3nf - no_3nf:.12e}",
            ])


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--no-3nf-dir", type=Path, required=True)
    parser.add_argument("--with-3nf-dir", type=Path, required=True)
    parser.add_argument("--experiment", type=Path, required=True)
    parser.add_argument("--target-tlab-mev", type=float, default=10.0)
    parser.add_argument("--output", type=Path, required=True,
                        help="output path without extension")
    parser.add_argument("--binding-energy-mev", type=float, default=float("nan"))
    args = parser.parse_args()

    theta = np.arange(0.0, 181.0, 1.0)
    exp_theta, exp_ay, exp_err = load_experiment(args.experiment)
    no_3nf = interpolate_curve(args.no_3nf_dir, args.target_tlab_mev, theta)
    with_3nf = interpolate_curve(args.with_3nf_dir, args.target_tlab_mev, theta)
    if (abs(no_3nf.lower_tlab - with_3nf.lower_tlab) > 1e-9
            or abs(no_3nf.upper_tlab - with_3nf.upper_tlab) > 1e-9):
        raise ValueError("2NF and 3NF runs do not use the same energy bracket")

    rmse_2nf = rmse_at_experiment(theta, no_3nf.curve, exp_theta, exp_ay)
    rmse_3nf = rmse_at_experiment(theta, with_3nf.curve, exp_theta, exp_ay)
    delta = with_3nf.curve - no_3nf.curve

    fig, (axis, delta_axis) = plt.subplots(
        2, 1, figsize=(10.2, 8.2), sharex=True,
        gridspec_kw={"height_ratios": [3.1, 1.15], "hspace": 0.08},
    )
    axis.errorbar(
        exp_theta, exp_ay, yerr=exp_err, fmt="o", color="black",
        markersize=5.5, capsize=3, linewidth=1.1,
        label="Experiment (digitized EXFOR markers)", zorder=4,
    )
    axis.plot(
        theta, no_3nf.curve, color="#4b5563", linewidth=2.4,
        label=rf"2NF only  (RMSE={rmse_2nf:.3f})",
    )
    axis.plot(
        theta, with_3nf.curve, color="#dc2626", linewidth=2.4,
        label=rf"2NF + approximate 3NF  (RMSE={rmse_3nf:.3f})",
    )
    axis.axhline(0.0, color="#9ca3af", linewidth=0.8)
    axis.set_ylabel(r"$A_y(n)$")
    axis.set_xlim(0.0, 180.0)
    axis.grid(alpha=0.25)
    axis.legend(loc="best", framealpha=0.95)

    delta_axis.plot(theta, delta, color="#7c3aed", linewidth=2.1)
    delta_axis.fill_between(theta, 0.0, delta, color="#7c3aed", alpha=0.18)
    delta_axis.axhline(0.0, color="black", linewidth=0.8)
    delta_axis.set_xlabel(r"$\theta_{c.m.}$ (deg)")
    delta_axis.set_ylabel(r"$\Delta A_y$")
    delta_axis.grid(alpha=0.25)

    bracket = (
        rf"linear energy interpolation: $T_{{lab}}={no_3nf.lower_tlab:.2f}$"
        rf" $\leftrightarrow$ {no_3nf.upper_tlab:.2f} MeV"
    )
    binding = (
        f", E_d={args.binding_energy_mev:.3f} MeV"
        if math.isfinite(args.binding_energy_mev) else ""
    )
    warning = (
        "DIAGNOSTIC, NOT A CONVERGED PREDICTION\n"
        rf"$N_p=N_q=10$, $J_{{3N}}\leq 1/2$, W1 cell N=2{binding}; "
        "most Padé elements are max-order truncated."
    )
    fig.suptitle(
        rf"Effect of the current approximate 3NF on $nd$ $A_y(n)$ at "
        rf"$T_{{lab}}={args.target_tlab_mev:.1f}$ MeV" + "\n" + bracket,
        fontsize=14,
    )
    axis.text(
        0.02, 0.03, warning, transform=axis.transAxes, fontsize=9,
        va="bottom", ha="left",
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "#fff7ed",
              "edgecolor": "#f97316", "alpha": 0.95},
    )
    fig.text(
        0.5, 0.012,
        "Experimental points are draft values digitized from Miller et al., PRC 107, 014002 (2023), Fig. 6; stated digitization precision ≈ ±0.01.",
        ha="center", fontsize=8.5, color="#4b5563",
    )
    fig.subplots_adjust(top=0.86, bottom=0.11)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    png_path = args.output.with_suffix(".png")
    svg_path = args.output.with_suffix(".svg")
    csv_path = args.output.with_suffix(".csv")
    json_path = args.output.with_suffix(".json")
    fig.savefig(png_path, dpi=190, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    write_csv(csv_path, theta, no_3nf.curve, with_3nf.curve)

    summary = {
        "target_tlab_mev": args.target_tlab_mev,
        "energy_interpolation": {
            "lower_tlab_mev": no_3nf.lower_tlab,
            "upper_tlab_mev": no_3nf.upper_tlab,
            "upper_weight": no_3nf.interpolation_weight,
        },
        "rmse_vs_digitized_experiment": {"2nf": rmse_2nf, "2nf_plus_3nf": rmse_3nf},
        "delta_Ay": {
            "max_abs": float(np.max(np.abs(delta))),
            "theta_at_max_abs_deg": float(theta[int(np.argmax(np.abs(delta)))]),
            "signed_value_at_max_abs": float(delta[int(np.argmax(np.abs(delta)))]),
        },
        "pade_counts": {
            "2nf": no_3nf.convergence_counts,
            "2nf_plus_3nf": with_3nf.convergence_counts,
        },
        "diagnostic_only": True,
        "binding_energy_mev": args.binding_energy_mev,
        "experiment": str(args.experiment),
        "outputs": {"png": str(png_path), "svg": str(svg_path), "csv": str(csv_path)},
    }
    json_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
