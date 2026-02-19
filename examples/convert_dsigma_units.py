#!/usr/bin/env python3
"""Convert differential cross section table units between mb/sr and fm^2/sr."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

from observable_units import (
    SUPPORTED_DSIGMA_UNITS,
    UNIT_MB_PER_SR,
    convert_dsigma_value,
    infer_dsigma_unit_from_lines,
    normalize_dsigma_unit,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Convert dSigma/dOmega data units")
    parser.add_argument(
        "--input",
        default="data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt",
        help="Input text file with angle and dSigma columns",
    )
    parser.add_argument(
        "--output",
        default="data/DataOfCrosssectionAndPol/DSigamaOverDOmega_fm2_per_sr.txt",
        help="Output text file path",
    )
    parser.add_argument(
        "--to-unit",
        default="fm2/sr",
        choices=list(SUPPORTED_DSIGMA_UNITS),
        help="Target dSigma unit",
    )
    parser.add_argument(
        "--from-unit",
        default="",
        help="Optional source unit override (auto-detected by default)",
    )
    return parser


def _replace_unit_marker(line: str, target_unit: str) -> str:
    lower = line.lower()
    if "unit" not in lower:
        return line
    unit_print = "fm^2/sr" if target_unit == "fm2/sr" else "mb/sr"
    if "unit " in lower:
        prefix = line[: lower.index("unit ") + len("unit ")]
        return f"{prefix}{unit_print}"
    return line


def convert_file(input_path: Path, output_path: Path, from_unit: str, to_unit: str) -> int:
    lines = input_path.read_text(encoding="utf-8").splitlines()
    out_lines: List[str] = []
    out_lines.append(
        f"# Converted with examples/convert_dsigma_units.py: {from_unit} -> {to_unit}"
    )

    converted_rows = 0
    for raw_line in lines:
        line = raw_line.rstrip("\n")
        stripped = line.strip()
        if not stripped:
            out_lines.append(line)
            continue

        parts = stripped.split()
        if len(parts) >= 2:
            try:
                angle = float(parts[0])
                value = float(parts[1])
                converted = convert_dsigma_value(value, from_unit, to_unit)
                out_lines.append(f"{angle:.14e}    {converted:.14e}")
                converted_rows += 1
                continue
            except ValueError:
                pass

        out_lines.append(_replace_unit_marker(line, to_unit))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")
    return converted_rows


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    input_path = (root / args.input).resolve()
    output_path = (root / args.output).resolve()

    to_unit = normalize_dsigma_unit(args.to_unit)
    source_lines = input_path.read_text(encoding="utf-8").splitlines()
    from_unit = (
        normalize_dsigma_unit(args.from_unit)
        if args.from_unit.strip()
        else infer_dsigma_unit_from_lines(source_lines, fallback_unit=UNIT_MB_PER_SR)
    )

    count = convert_file(input_path, output_path, from_unit, to_unit)
    print(f"input: {input_path}")
    print(f"output: {output_path}")
    print(f"source_unit: {from_unit}")
    print(f"target_unit: {to_unit}")
    print(f"converted_rows: {count}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
