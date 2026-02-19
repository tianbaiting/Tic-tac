#!/usr/bin/env python3
"""
Purpose:
  Shared unit helpers for dSigma/dOmega tables used by example workflows.

Data flow (left -> right):
  raw unit token or file header lines
    -> normalized canonical unit (`mb/sr` or `fm2/sr`)
    -> scale factor
    -> converted scalar or converted sequence

Called by:
  - `convert_dsigma_units.py`
  - `compare_Ay_experiment.py`
  - `run_dpol_p_observables.py`

Usage:
  Import as a utility module:
    from observable_units import convert_dsigma_value, normalize_dsigma_unit
"""

from __future__ import annotations

import re
from typing import Iterable, List


UNIT_MB_PER_SR = "mb/sr"
UNIT_FM2_PER_SR = "fm2/sr"
SUPPORTED_DSIGMA_UNITS = (UNIT_MB_PER_SR, UNIT_FM2_PER_SR)

_MB_TO_FM2 = 0.1


def normalize_dsigma_unit(unit: str) -> str:
    token = unit.strip().lower().replace(" ", "")
    aliases = {
        "mb": UNIT_MB_PER_SR,
        "mb/sr": UNIT_MB_PER_SR,
        "mbsr": UNIT_MB_PER_SR,
        "mb_per_sr": UNIT_MB_PER_SR,
        "fm2": UNIT_FM2_PER_SR,
        "fm^2": UNIT_FM2_PER_SR,
        "fm2/sr": UNIT_FM2_PER_SR,
        "fm^2/sr": UNIT_FM2_PER_SR,
        "fm2sr": UNIT_FM2_PER_SR,
        "fm^2sr": UNIT_FM2_PER_SR,
        "fm2_per_sr": UNIT_FM2_PER_SR,
    }
    if token not in aliases:
        raise ValueError(f"Unsupported dSigma unit: {unit}")
    return aliases[token]


def dsigma_unit_factor(from_unit: str, to_unit: str) -> float:
    src = normalize_dsigma_unit(from_unit)
    dst = normalize_dsigma_unit(to_unit)
    if src == dst:
        return 1.0
    if src == UNIT_MB_PER_SR and dst == UNIT_FM2_PER_SR:
        return _MB_TO_FM2
    if src == UNIT_FM2_PER_SR and dst == UNIT_MB_PER_SR:
        return 1.0 / _MB_TO_FM2
    raise ValueError(f"Unsupported dSigma conversion: {from_unit} -> {to_unit}")


def convert_dsigma_value(value: float, from_unit: str, to_unit: str) -> float:
    return value * dsigma_unit_factor(from_unit, to_unit)


def convert_dsigma_series(values: Iterable[float], from_unit: str, to_unit: str) -> List[float]:
    factor = dsigma_unit_factor(from_unit, to_unit)
    return [factor * value for value in values]


def infer_dsigma_unit_from_lines(lines: Iterable[str], fallback_unit: str = UNIT_MB_PER_SR) -> str:
    pattern = re.compile(r"unit\s+([A-Za-z0-9^/_]+)", flags=re.IGNORECASE)
    for raw_line in lines:
        line = raw_line.strip()
        if not line:
            continue
        lower = line.lower()
        if "unit" not in lower:
            continue
        match = pattern.search(line)
        if match:
            try:
                return normalize_dsigma_unit(match.group(1))
            except ValueError:
                continue
        if "mb/sr" in lower:
            return UNIT_MB_PER_SR
        if "fm^2/sr" in lower or "fm2/sr" in lower:
            return UNIT_FM2_PER_SR
    return normalize_dsigma_unit(fallback_unit)
