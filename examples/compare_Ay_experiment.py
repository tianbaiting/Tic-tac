#!/usr/bin/env python3
"""
190 MeV dpol-p solver-vs-experiment comparison.

Comparison is based on solver-produced U-matrix elements from Tic-tac.
No experimental-data interpolation is used in the model prediction path.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional


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


def _residual_metrics(pred: float, obs: List[float]) -> Dict[str, float]:
    residuals = [pred - val for val in obs]
    n = len(residuals)
    mae = sum(abs(val) for val in residuals) / n if n else math.nan
    rmse = math.sqrt(sum(val * val for val in residuals) / n) if n else math.nan
    max_abs = max(abs(val) for val in residuals) if residuals else math.nan
    return {
        "count": float(n),
        "mae": mae,
        "rmse": rmse,
        "max_abs_error": max_abs,
    }


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
        if dsigma_proxy > 1e-15:
            ay_proxy = (f_no_flip.conjugate() * f_flip).imag / dsigma_proxy
        else:
            ay_proxy = 0.0

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


def combine_channels_by_energy(points: List[SolverChannelPoint]) -> List[Dict[str, float]]:
    grouped: Dict[float, List[SolverChannelPoint]] = {}
    for item in points:
        key = round(item.tlab, 6)
        grouped.setdefault(key, []).append(item)

    combined: List[Dict[str, float]] = []
    for key in sorted(grouped.keys()):
        block = grouped[key]
        total_w = sum(item.dsigma_proxy for item in block)
        if total_w <= 1e-20:
            continue

        ay_combined = sum(item.ay_proxy * item.dsigma_proxy for item in block) / total_w
        ecm_combined = sum(item.ecm * item.dsigma_proxy for item in block) / total_w

        combined.append(
            {
                "tlab": key,
                "ecm_weighted": ecm_combined,
                "ay_proxy_combined": ay_combined,
                "dsigma_proxy_combined": total_w,
                "num_channels": float(len(block)),
            }
        )

    return combined


def run_solver_if_requested(work_dir: Path, regenerate: bool) -> None:
    if not regenerate:
        return

    cmd = [
        "python3",
        "examples/deuteron_proton_Ay.py",
        "--work-dir",
        str(work_dir),
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
    lines: List[str] = []
    lines.append("190 MeV dpol-p: Tic-tac solver vs experiment")
    lines.append("============================================")
    lines.append(f"status: {summary['status']}")
    lines.append(f"reason: {summary['status_reason']}")
    lines.append(f"target_tlab: {summary['target_tlab']:.3f} MeV")
    lines.append(f"best_tlab: {summary['best_energy']['tlab']:.3f} MeV")
    lines.append(f"abs_delta_tlab: {summary['best_energy']['delta_tlab']:.3f} MeV")
    lines.append("")

    lines.append("Best-energy proxy metrics:")
    lines.append(f"  ay_proxy_combined: {summary['best_energy']['ay_proxy_combined']:.6f}")
    lines.append(f"  dsigma_proxy_combined: {summary['best_energy']['dsigma_proxy_combined']:.6f}")
    lines.append(f"  iT11 mae: {summary['best_energy']['it11_mae']:.6f}")
    lines.append(f"  iT11 rmse: {summary['best_energy']['it11_rmse']:.6f}")
    lines.append(f"  iT11 max_abs_error: {summary['best_energy']['it11_max_abs_error']:.6f}")
    lines.append(f"  dSigma raw rmse: {summary['best_energy']['dsigma_raw_rmse']:.6f}")
    lines.append(f"  dSigma norm-shape rmse: {summary['best_energy']['dsigma_shape_rmse']:.6f}")
    lines.append("")

    lines.append("Thresholds:")
    lines.append(f"  ay_rmse_pass: {summary['thresholds']['ay_rmse_pass']:.6f}")
    lines.append(f"  ay_mae_pass: {summary['thresholds']['ay_mae_pass']:.6f}")
    lines.append(f"  max_energy_delta_pass: {summary['thresholds']['energy_delta_pass']:.6f}")
    lines.append("")

    lines.append("All available solver energies:")
    for item in summary["energies"]:
        lines.append(
            "  "
            f"Tlab={item['tlab']:.3f} MeV, "
            f"ay_proxy={item['ay_proxy_combined']:.6f}, "
            f"iT11_rmse={item['it11_rmse']:.6f}, "
            f"dSigma_shape_rmse={item['dsigma_shape_rmse']:.6f}"
        )

    lines.append("")
    lines.append("Note:")
    lines.append("  This validation uses solver-produced U-matrix elements only.")
    lines.append("  No interpolation of experimental observable curves is used in model prediction.")

    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Validate Tic-tac solver output against 190MeV dpol-p data")
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
    parser.add_argument("--target-tlab", type=float, default=135.6)
    parser.add_argument("--regenerate", action="store_true", help="Run solver before comparison")
    parser.add_argument("--ay-rmse-pass", type=float, default=0.12)
    parser.add_argument("--ay-mae-pass", type=float, default=0.10)
    parser.add_argument("--energy-delta-pass", type=float, default=40.0)
    return parser


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]

    work_dir = (root / args.work_dir).resolve()
    solver_out_dir = (root / args.solver_out_dir).resolve()
    tensor_data = (root / args.tensor_data).resolve()
    dsigma_data = (root / args.dsigma_data).resolve()

    work_dir.mkdir(parents=True, exist_ok=True)
    run_solver_if_requested(work_dir, args.regenerate)

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

    exp_it11_values = exp_it11["values"]
    exp_dsigma_values = exp_dsigma["values"]
    exp_dsigma_mean = sum(exp_dsigma_values) / len(exp_dsigma_values)
    exp_dsigma_norm = [val / exp_dsigma_mean for val in exp_dsigma_values]

    enriched: List[Dict[str, float]] = []
    for item in combined:
        it11_metrics = _residual_metrics(item["ay_proxy_combined"], exp_it11_values)

        dsigma_metrics_raw = _residual_metrics(item["dsigma_proxy_combined"], exp_dsigma_values)
        dsigma_shape_rmse = math.sqrt(
            sum((1.0 - val) * (1.0 - val) for val in exp_dsigma_norm) / len(exp_dsigma_norm)
        )

        enriched_item = dict(item)
        enriched_item["delta_tlab"] = abs(item["tlab"] - args.target_tlab)
        enriched_item["it11_mae"] = it11_metrics["mae"]
        enriched_item["it11_rmse"] = it11_metrics["rmse"]
        enriched_item["it11_max_abs_error"] = it11_metrics["max_abs_error"]
        enriched_item["dsigma_raw_rmse"] = dsigma_metrics_raw["rmse"]
        enriched_item["dsigma_shape_rmse"] = dsigma_shape_rmse
        enriched.append(enriched_item)

    best = min(enriched, key=lambda item: item["delta_tlab"])

    pass_flag = (
        best["it11_rmse"] <= args.ay_rmse_pass
        and best["it11_mae"] <= args.ay_mae_pass
        and best["delta_tlab"] <= args.energy_delta_pass
    )

    summary = {
        "status": "pass" if pass_flag else "fail",
        "status_reason": (
            "best solver proxy is within thresholds"
            if pass_flag
            else "best solver proxy exceeds thresholds"
        ),
        "target_tlab": args.target_tlab,
        "thresholds": {
            "ay_rmse_pass": args.ay_rmse_pass,
            "ay_mae_pass": args.ay_mae_pass,
            "energy_delta_pass": args.energy_delta_pass,
        },
        "inputs": {
            "solver_out_dir": str(solver_out_dir),
            "u_files": [str(p) for p in u_files],
            "tensor_data": str(tensor_data),
            "dsigma_data": str(dsigma_data),
        },
        "best_energy": best,
        "energies": enriched,
    }

    json_file = work_dir / "solver_validation_190MeV.json"
    txt_file = work_dir / "solver_validation_190MeV.txt"

    json_file.write_text(json.dumps(summary, indent=2, ensure_ascii=False), encoding="utf-8")
    txt_file.write_text(build_report_text(summary), encoding="utf-8")

    print(build_report_text(summary), end="")
    print(f"json_report: {json_file}")
    print(f"text_report: {txt_file}")

    return 0 if pass_flag else 2


if __name__ == "__main__":
    raise SystemExit(main())
