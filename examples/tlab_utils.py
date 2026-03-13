#!/usr/bin/env python3

from __future__ import annotations

import math


def format_tlab_label(tlab_mev: float) -> str:
    nearest_int = round(tlab_mev)
    if math.isclose(tlab_mev, float(nearest_int), rel_tol=0.0, abs_tol=1e-9):
        return f"{int(nearest_int)}MeV"

    token = f"{tlab_mev:.6f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"{token}MeV"


def format_tlab_dir_label(tlab_mev: float) -> str:
    nearest_int = round(tlab_mev)
    if math.isclose(tlab_mev, float(nearest_int), rel_tol=0.0, abs_tol=1e-9):
        return f"{int(nearest_int):03d}MeV"

    token = f"{tlab_mev:.6f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"{token}MeV"


def format_tlab_dir_name(tlab_mev: float) -> str:
    return f"tlab_{format_tlab_dir_label(tlab_mev)}"
