#!/usr/bin/env python3
"""
Purpose:
  Utility helpers for parsing/choosing Tic-tac `U_PW_elements_*.txt` files.

Data flow (left -> right):
  filename string
    -> parsed metadata (Np, Nq, JP, Jmax, PSI)
    -> filtering/grouping logic
    -> selected file family for downstream physics scripts

Called by:
  - `compare_Ay_experiment.py`
  - `run_dpol_p_observables.py`

Usage:
  Import as a helper module:
    from solver_u_file_utils import select_latest_u_file_family, detect_parity_symbol
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple


_U_FILE_RE = re.compile(
    r"^U_PW_elements_Np_(?P<np>\d+)_Nq_(?P<nq>\d+)_JP_(?P<two_j>-?\d+)_(?P<parity>-?1)_Jmax_(?P<jmax>\d+)_PSI_(?P<psi>-?\d+)\.txt$"
)


@dataclass(frozen=True)
class UFileMeta:
    np_wp: int
    nq_wp: int
    two_j: int
    parity: int
    j2max: int
    psi: int


def parse_u_file_meta(path: Path) -> Optional[UFileMeta]:
    match = _U_FILE_RE.match(path.name)
    if match is None:
        return None
    return UFileMeta(
        np_wp=int(match.group("np")),
        nq_wp=int(match.group("nq")),
        two_j=int(match.group("two_j")),
        parity=int(match.group("parity")),
        j2max=int(match.group("jmax")),
        psi=int(match.group("psi")),
    )


def detect_parity_symbol(path: Path) -> str:
    meta = parse_u_file_meta(path)
    if meta is None:
        return "?"
    if meta.parity == 1:
        return "+"
    if meta.parity == -1:
        return "-"
    return "?"


def detect_two_j(path: Path) -> Optional[int]:
    meta = parse_u_file_meta(path)
    if meta is None:
        return None
    return meta.two_j


def select_latest_u_file_family(solver_out_dir: Path) -> List[Path]:
    # "Family" means same (Np, Nq, Jmax/J2max tag, PSI), differing only in JP/parity.
    candidates = list(solver_out_dir.glob("U_PW_elements_*.txt"))
    parsed: List[Tuple[Path, UFileMeta]] = []
    for path in candidates:
        meta = parse_u_file_meta(path)
        if meta is not None:
            parsed.append((path, meta))
    if not parsed:
        return []

    latest_path, latest_meta = max(parsed, key=lambda item: item[0].stat().st_mtime)
    family_key = (latest_meta.np_wp, latest_meta.nq_wp, latest_meta.j2max, latest_meta.psi)

    selected: List[Tuple[Path, UFileMeta]] = []
    for path, meta in parsed:
        key = (meta.np_wp, meta.nq_wp, meta.j2max, meta.psi)
        if key == family_key:
            selected.append((path, meta))

    selected.sort(key=lambda item: (item[1].two_j, item[1].parity, item[0].name))
    return [path for path, _ in selected]


def required_p123_sparse_names(np_wp: int, nq_wp: int, j2max: int, two_j_3n_max: int) -> List[str]:
    # Build exact sparse P123 filenames expected by solver when reusing P123.
    names: List[str] = []
    for two_j in range(1, two_j_3n_max + 1, 2):
        for parity in (1, -1):
            names.append(f"P123_sparse_JP_{two_j}_{parity}_Np_{np_wp}_Nq_{nq_wp}_J2max_{j2max}.h5")
    return names
