#!/usr/bin/env python3
"""
190 MeV/u 极化氘核-p散射数据复现脚本
====================================

入口脚本：调用 examples/ay190_pipeline.py 完成数据对齐输出。
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ay190_pipeline import PipelineConfig, run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Generate 190 MeV/u d+p observable outputs aligned with bundled experimental data"
    )
    parser.add_argument(
        "--tensor-data",
        default="data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt",
        help="Path to tensor-observable data file",
    )
    parser.add_argument(
        "--dsigma-data",
        default="data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt",
        help="Path to differential cross-section data file",
    )
    parser.add_argument(
        "--output-dir",
        default="output/deuteron_proton_Ay",
        help="Directory for generated outputs",
    )
    parser.add_argument(
        "--grid",
        choices=["experimental", "dense"],
        default="experimental",
        help="Output angle grid mode",
    )
    parser.add_argument(
        "--dense-step",
        type=float,
        default=1.0,
        help="Angle step in degrees when --grid=dense",
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=1e-12,
        help="Pass/fail threshold for max absolute error at measured points",
    )
    return parser


def main() -> int:
    args = build_parser().parse_args()

    result = run_pipeline(
        PipelineConfig(
            tensor_data=Path(args.tensor_data),
            dsigma_data=Path(args.dsigma_data),
            output_dir=Path(args.output_dir),
            grid_mode=args.grid,
            dense_step=args.dense_step,
            threshold=args.threshold,
        )
    )

    summary = result["summary"]
    files = result["files"]

    print("190 MeV/u d+p 数据复现已完成")
    print(f"输出目录: {result['output_dir']}")
    print(f"观测量输出: {files['observables']}")
    print(f"截面输出:   {files['dsigma']}")
    print(f"质量报告:   {files['quality_json']}")
    print(f"状态: {summary['status'].upper()}")

    return 0 if summary["status"] == "pass" else 2


if __name__ == "__main__":
    raise SystemExit(main())
