#!/usr/bin/env python3
"""Plot experiment vs model validation curves using matplotlib."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_curve_csv(path: Path) -> tuple[list[float], list[float], list[float], str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if len(header) != 3:
            raise ValueError(f"Unexpected header in {path}: {header}")
        theta_col, exp_col, model_col = header

        theta: list[float] = []
        exp_vals: list[float] = []
        model_vals: list[float] = []
        for row in reader:
            if len(row) != 3:
                continue
            theta.append(float(row[0]))
            exp_vals.append(float(row[1]))
            model_vals.append(float(row[2]))

    if not theta:
        raise ValueError(f"No data rows in {path}")

    return theta, exp_vals, model_vals, exp_col, model_col


def make_plot(
    *,
    theta: list[float],
    exp_vals: list[float],
    model_vals: list[float],
    exp_label: str,
    model_label: str,
    title: str,
    y_label: str,
    out_file: Path,
    log_y: bool,
) -> None:
    fig = plt.figure(figsize=(9.5, 6.2), dpi=140)
    ax = fig.add_subplot(1, 1, 1)

    ax.scatter(theta, exp_vals, s=24, color="#111827", label=f"Experiment ({exp_label})", zorder=3)
    ax.plot(theta, model_vals, color="#dc2626", linewidth=2.2, label=f"Model ({model_label})", zorder=2)

    ax.set_title(title)
    ax.set_xlabel("theta_cm (deg)")
    ax.set_ylabel(y_label)
    ax.grid(True, alpha=0.25)
    if log_y:
        ax.set_yscale("log")
    ax.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_file)
    plt.close(fig)


def make_combined_plot(
    *,
    it11_theta: list[float],
    it11_exp: list[float],
    it11_model: list[float],
    ds_theta: list[float],
    ds_exp: list[float],
    ds_model: list[float],
    out_file: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.4), dpi=140)

    ax0 = axes[0]
    ax0.scatter(it11_theta, it11_exp, s=20, color="#111827", label="Experiment", zorder=3)
    ax0.plot(it11_theta, it11_model, color="#dc2626", linewidth=2.0, label="Model", zorder=2)
    ax0.set_title("iT11 Comparison")
    ax0.set_xlabel("theta_cm (deg)")
    ax0.set_ylabel("iT11")
    ax0.grid(True, alpha=0.25)
    ax0.legend(loc="best")

    ax1 = axes[1]
    ax1.scatter(ds_theta, ds_exp, s=20, color="#111827", label="Experiment", zorder=3)
    ax1.plot(ds_theta, ds_model, color="#dc2626", linewidth=2.0, label="Model", zorder=2)
    ax1.set_yscale("log")
    ax1.set_title("dSigma/dOmega Comparison")
    ax1.set_xlabel("theta_cm (deg)")
    ax1.set_ylabel("dSigma/dOmega")
    ax1.grid(True, alpha=0.25)
    ax1.legend(loc="best")

    fig.tight_layout()
    fig.savefig(out_file)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Plot model-vs-experiment curves from validation CSV outputs")
    parser.add_argument("--work-dir", default="output/deuteron_proton_Ay", help="Validation output directory")
    parser.add_argument("--it11-csv", default="best_energy_iT11_curve.csv")
    parser.add_argument("--dsigma-csv", default="best_energy_dsigma_curve.csv")
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    work_dir = (root / args.work_dir).resolve()

    it11_csv = work_dir / args.it11_csv
    dsigma_csv = work_dir / args.dsigma_csv

    it11_theta, it11_exp, it11_model, it11_exp_label, it11_model_label = read_curve_csv(it11_csv)
    ds_theta, ds_exp, ds_model, ds_exp_label, ds_model_label = read_curve_csv(dsigma_csv)

    it11_png = work_dir / "best_energy_iT11_comparison.png"
    dsigma_png = work_dir / "best_energy_dsigma_comparison.png"
    combined_png = work_dir / "best_energy_comparison.png"

    make_plot(
        theta=it11_theta,
        exp_vals=it11_exp,
        model_vals=it11_model,
        exp_label=it11_exp_label,
        model_label=it11_model_label,
        title="190 MeV/u dpol-p: iT11 (Experiment vs Model)",
        y_label="iT11",
        out_file=it11_png,
        log_y=False,
    )
    make_plot(
        theta=ds_theta,
        exp_vals=ds_exp,
        model_vals=ds_model,
        exp_label=ds_exp_label,
        model_label=ds_model_label,
        title="190 MeV/u dpol-p: dSigma/dOmega (Experiment vs Model)",
        y_label="dSigma/dOmega",
        out_file=dsigma_png,
        log_y=True,
    )
    make_combined_plot(
        it11_theta=it11_theta,
        it11_exp=it11_exp,
        it11_model=it11_model,
        ds_theta=ds_theta,
        ds_exp=ds_exp,
        ds_model=ds_model,
        out_file=combined_png,
    )

    print(f"it11_png: {it11_png}")
    print(f"dsigma_png: {dsigma_png}")
    print(f"combined_png: {combined_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
