#!/usr/bin/env python3
"""
190 MeV/u 极化氘核-p数据对比脚本
==============================

读取 examples/deuteron_proton_Ay.py 生成的质量报告，
输出简明的数值对比结果和可归档文本报告。
"""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path
from typing import Dict, Tuple


def run_generator_if_needed(output_dir: Path, force_regenerate: bool) -> None:
    summary_path = output_dir / "fit_quality_190MeV.json"
    if summary_path.exists() and not force_regenerate:
        return

    script_path = Path(__file__).with_name("deuteron_proton_Ay.py")
    cmd = [
        "python3",
        str(script_path),
        "--output-dir",
        str(output_dir),
        "--grid",
        "experimental",
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(
            "Failed to generate baseline outputs.\n"
            f"cmd: {' '.join(cmd)}\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def load_summary(path: Path) -> Dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def format_metric_row(name: str, payload: Dict[str, float]) -> str:
    count = int(round(payload.get("count", 0.0)))
    max_abs_error = payload.get("max_abs_error", float("nan"))
    rmse = payload.get("rmse", float("nan"))
    chi2_per_point = payload.get("chi2_per_point", float("nan"))
    return (
        f"{name:<14} points={count:>3d}  "
        f"max_abs_error={max_abs_error:>11.3e}  "
        f"rmse={rmse:>11.3e}  "
        f"chi2/pt={chi2_per_point:>11.3e}"
    )


def build_report_text(summary: Dict[str, object]) -> str:
    lines = []
    lines.append("190 MeV/u d+p 数据对比报告")
    lines.append("============================")
    lines.append(f"status: {summary['status']}")
    lines.append(f"threshold: {summary['threshold']:.3e}")
    lines.append(f"max_abs_error_all: {summary['max_abs_error_all']:.3e}")
    lines.append(f"grid_mode: {summary['grid_mode']}")

    inputs = summary["inputs"]
    lines.append(f"tensor_data: {inputs['tensor_data']}")
    lines.append(f"dsigma_data: {inputs['dsigma_data']}")
    lines.append("")

    lines.append("指标明细:")
    metrics = summary["metrics"]
    for name in ["iT11", "T20", "T21", "T22", "dSigma_dOmega"]:
        if name in metrics:
            lines.append(format_metric_row(name, metrics[name]))

    lines.append("")
    lines.append("说明:")
    lines.append("- 当前流程用于验证仓库内190MeV数据读写和输出一致性。")
    lines.append("- 该输出是数据对齐基线，不等同于完整Faddeev物理预言。")

    return "\n".join(lines) + "\n"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compare generated 190 MeV/u outputs with bundled experimental data"
    )
    parser.add_argument(
        "--output-dir",
        default="output/deuteron_proton_Ay",
        help="Directory containing generated 190 MeV outputs",
    )
    parser.add_argument(
        "--regenerate",
        action="store_true",
        help="Force rerun of examples/deuteron_proton_Ay.py before comparison",
    )
    parser.add_argument(
        "--report-file",
        default="comparison_report_190MeV.txt",
        help="Report filename written inside output-dir",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    run_generator_if_needed(output_dir, args.regenerate)

    summary_path = output_dir / "fit_quality_190MeV.json"
    summary = load_summary(summary_path)

    report_text = build_report_text(summary)
    report_path = output_dir / args.report_file
    with report_path.open("w", encoding="utf-8") as handle:
        handle.write(report_text)

    print(report_text, end="")
    print(f"报告文件: {report_path}")

    return 0 if summary["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
