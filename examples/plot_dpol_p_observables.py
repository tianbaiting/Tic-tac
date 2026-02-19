#!/usr/bin/env python3
"""Plot solver-driven dpol-p observables for 70/135/190 MeV/u style workflows."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


OBS_KEYS = ["dSigma_dOmega", "iT11", "T20", "T21", "T22"]
OBS_LABELS = {
    "dSigma_dOmega": r"$d\sigma/d\Omega$",
    "iT11": "iT11",
    "T20": "T20",
    "T21": "T21",
    "T22": "T22",
}


def read_model_csv(path: Path) -> Dict[str, List[float]]:
    out = {k: [] for k in ["theta_deg"] + OBS_KEYS}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            out["theta_deg"].append(float(row["theta_deg"]))
            for k in OBS_KEYS:
                out[k].append(float(row[k]))
    return out


def read_experiment_csv(path: Path) -> Dict[str, List[float]]:
    out = {k: [] for k in ["theta_deg"] + [f"{x}_exp" for x in OBS_KEYS]}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            out["theta_deg"].append(float(row["theta_deg"]))
            for k in OBS_KEYS:
                v = row.get(f"{k}_exp", "")
                if v is None or v.strip() == "":
                    out[f"{k}_exp"].append(math.nan)
                else:
                    out[f"{k}_exp"].append(float(v))
    return out


def _valid_pairs(xs: List[float], ys: List[float]) -> Tuple[List[float], List[float]]:
    px: List[float] = []
    py: List[float] = []
    for x, y in zip(xs, ys):
        if math.isfinite(y):
            px.append(x)
            py.append(y)
    return px, py


def _finite_values(values: List[float]) -> List[float]:
    return [v for v in values if math.isfinite(v)]


def _auto_linear_ylim(values: List[float]) -> Tuple[float, float]:
    finite = _finite_values(values)
    if not finite:
        return (-1.0, 1.0)

    y_min = min(finite)
    y_max = max(finite)
    if y_max <= y_min:
        center = y_min
        half_span = max(0.05, 0.1 * max(1.0, abs(center)))
        return (center - half_span, center + half_span)

    span = y_max - y_min
    pad = 0.10 * span
    return (y_min - pad, y_max + pad)


def plot_single_energy(
    *,
    model: Dict[str, List[float]],
    exp: Optional[Dict[str, List[float]]],
    target_tlab: float,
    solver_tlab: float,
    out_file: Path,
    dsigma_unit_text: str,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), dpi=140)
    ax_list = [axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1]]
    axes[1, 2].axis("off")

    fig.suptitle(
        f"d + p observables (target {target_tlab:.1f} MeV/u, solver {solver_tlab:.3f} MeV)",
        fontsize=14,
        fontweight="bold",
    )

    theta = model["theta_deg"]
    for ax, key in zip(ax_list, OBS_KEYS):
        ax.plot(theta, model[key], color="#dc2626", lw=2.0, label="Faddeev simulation")
        if exp is not None:
            ex, ey = _valid_pairs(exp["theta_deg"], exp[f"{key}_exp"])
            if ex:
                ax.scatter(ex, ey, s=20, color="#111827", alpha=0.9, label="Experiment")
        if key == "dSigma_dOmega":
            ax.set_title(f"{OBS_LABELS[key]} [{dsigma_unit_text}]")
        else:
            ax.set_title(OBS_LABELS[key])
        ax.set_xlabel(r"$\theta_{\mathrm{c.m.}}$ (deg)")
        ax.grid(True, alpha=0.25)
        if key == "dSigma_dOmega":
            ax.set_yscale("log")
        else:
            y_for_limits = _finite_values(model[key])
            if exp is not None:
                y_for_limits.extend(_finite_values(exp[f"{key}_exp"]))
            ax.set_ylim(*_auto_linear_ylim(y_for_limits))
        ax.legend(loc="best", fontsize=9)

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file)
    plt.close(fig)


def plot_overview(
    *,
    all_model: List[Tuple[float, float, Dict[str, List[float]]]],
    out_file: Path,
    dsigma_unit_text: str,
) -> None:
    fig, axes = plt.subplots(2, 3, figsize=(14, 7), dpi=140)
    ax_list = [axes[0, 0], axes[0, 1], axes[0, 2], axes[1, 0], axes[1, 1]]
    axes[1, 2].axis("off")

    fig.suptitle("d + p observables from Tic-tac Faddeev U-matrix", fontsize=14, fontweight="bold")

    colors = ["#dc2626", "#1d4ed8", "#059669", "#a16207", "#7c3aed"]
    for idx, (target_tlab, solver_tlab, model) in enumerate(all_model):
        color = colors[idx % len(colors)]
        label = f"target {target_tlab:.1f}, solver {solver_tlab:.3f} MeV"
        theta = model["theta_deg"]
        for ax, key in zip(ax_list, OBS_KEYS):
            ax.plot(theta, model[key], color=color, lw=1.8, label=label)

    overview_limits: Dict[str, Tuple[float, float]] = {}
    for key in OBS_KEYS:
        if key == "dSigma_dOmega":
            continue
        y_values: List[float] = []
        for _, _, model in all_model:
            y_values.extend(_finite_values(model[key]))
        overview_limits[key] = _auto_linear_ylim(y_values)

    for ax, key in zip(ax_list, OBS_KEYS):
        if key == "dSigma_dOmega":
            ax.set_title(f"{OBS_LABELS[key]} [{dsigma_unit_text}]")
        else:
            ax.set_title(OBS_LABELS[key])
        ax.set_xlabel(r"$\theta_{\mathrm{c.m.}}$ (deg)")
        ax.grid(True, alpha=0.25)
        if key == "dSigma_dOmega":
            ax.set_yscale("log")
        else:
            ax.set_ylim(*overview_limits[key])
        ax.legend(loc="best", fontsize=8)

    fig.tight_layout(rect=(0, 0, 1, 0.94))
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_file)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot dpol-p observables from run_dpol_p_observables outputs")
    parser.add_argument("--work-dir", default="output/dpol_p_observables", help="Output root from run_dpol_p_observables.py")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    work_dir = (root / args.work_dir).resolve()
    analysis_dir = work_dir / "analysis"
    figures_dir = work_dir / "figures"
    summary_path = analysis_dir / "summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(f"Missing summary file: {summary_path}")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    energy_entries = summary.get("energies", [])
    dsigma_unit = str(summary.get("units", {}).get("dSigma_dOmega", "mb/sr"))
    dsigma_unit_text = "fm^2/sr" if dsigma_unit == "fm2/sr" else dsigma_unit
    if not energy_entries:
        raise RuntimeError("No energy entries in summary.json")

    all_model: List[Tuple[float, float, Dict[str, List[float]]]] = []

    for entry in energy_entries:
        target_tlab = float(entry["target_tlab_mev_per_u"])
        solver_tlab = float(entry["selected_solver_tlab_mev"])
        energy_dir = Path(str(entry["analysis_dir"]))
        model_csv = energy_dir / "observables_model.csv"
        model = read_model_csv(model_csv)
        all_model.append((target_tlab, solver_tlab, model))

        exp_csv = energy_dir / "observables_experiment_190.csv"
        exp_data = read_experiment_csv(exp_csv) if exp_csv.exists() else None

        out_file = figures_dir / f"{energy_dir.name}_observables.png"
        plot_single_energy(
            model=model,
            exp=exp_data,
            target_tlab=target_tlab,
            solver_tlab=solver_tlab,
            out_file=out_file,
            dsigma_unit_text=dsigma_unit_text,
        )

    overview_file = figures_dir / "overview_observables_multi_energy.png"
    plot_overview(all_model=all_model, out_file=overview_file, dsigma_unit_text=dsigma_unit_text)

    print(f"summary: {summary_path}")
    for entry in energy_entries:
        e_dir = Path(str(entry["analysis_dir"]))
        print(f"energy_figure: {figures_dir / (e_dir.name + '_observables.png')}")
    print(f"overview_figure: {overview_file}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
