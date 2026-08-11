#!/usr/bin/env python3
"""Reconstruct U_PW_elements_*.txt files from stored Neumann-series terms
without waiting for the solver's full Padé loop to converge.

The solver writes `neumann_terms_*.txt` after every Neumann step n=k, giving
us the on-shell sequence {a_0, a_1, ..., a_N} per (alpha', alpha, q_idx).
We replicate the Padé [NM/NM](z=1) resummation and convergence logic from
`src/core/faddeev_solver/solve_faddeev.cpp::pade_approximant` and write a
U_PW file in the format consumed by `examples/pw_amplitudes.py`.

Usage:
    python3 examples/reconstruct_u_from_neumann.py \\
        --neumann-dir CPP/Output/labenpg_3NF_J3N9_par_chn4 \\
        --template     CPP/Output/miller_gate2_v4_J3N9_J2N3/U_PW_elements_Np_20_Nq_20_JP_5_1_Jmax_3_PSI_0.txt \\
        --output-dir   CPP/Output/labenpg_3NF_J3N9_partial

The template gives the U_PW header (channel labels l', 2j', l, 2j) which is
identical between 2NF and 2NF+3NF runs at matched truncation.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

REPO = Path(__file__).resolve().parents[1]

_HEADER_ROW_RE = re.compile(
    r"^\s*U(?P<i>\d+)(?P<j>\d+)\s+(?P<row>\d+)\s+(?P<col>\d+)\s+"
    r"(?P<lp>\d+)\s+(?P<two_jp>\d+)\s+(?P<l>\d+)\s+(?P<two_j>\d+)\s*$"
)

_NEUMANN_LINE_RE = re.compile(
    r"^(?P<re>[+\-0-9.eE]+)\s+(?P<im>[+\-0-9.eE]+)\s*#.*"
    r"alpha'-idx=(?P<ap>\d+)\s+alpha-idx=(?P<a>\d+)\s+q-idx=(?P<q>\d+)"
)


def parse_u_template(template_path: Path) -> Tuple[List[str], List[Tuple[int, int, int, int]]]:
    """Return (header_lines, channel_labels).

    channel_labels[i] = (lp, two_jp, l, two_j) for row i (i = Uii row-idx).
    We also return the raw header lines up to (but not including) the data rows
    so we can replay them when writing the reconstructed U_PW file.
    """
    text = template_path.read_text()
    lines = text.splitlines()
    header_lines: List[str] = []
    channels: Dict[int, Tuple[int, int, int, int]] = {}
    in_header = True
    for raw in lines:
        header_match = _HEADER_ROW_RE.match(raw)
        if header_match and header_match.group("i") == header_match.group("j"):
            i = int(header_match.group("i"))
            channels[i] = (
                int(header_match.group("lp")),
                int(header_match.group("two_jp")),
                int(header_match.group("l")),
                int(header_match.group("two_j")),
            )
        if in_header:
            header_lines.append(raw)
            # Heuristic: header ends at the last line starting with '#'
            if raw.lstrip().startswith("#") and "U00 [MeV]" in raw:
                in_header = False
    n_chan = max(channels) + 1
    chan_tuples = [
        (channels[i][0], channels[i][1], channels[i][2], channels[i][3])
        for i in range(n_chan)
    ]
    return header_lines, chan_tuples


def parse_neumann_file(path: Path) -> Dict[Tuple[int, int, int], List[complex]]:
    """Parse neumann_terms file into {(alpha', alpha, q_idx): [a_0, a_1, ...]}."""
    sequence: Dict[Tuple[int, int, int], List[complex]] = {}
    current_n = -1
    text = path.read_text()
    for raw in text.splitlines():
        stripped = raw.strip()
        if not stripped:
            continue
        m_n = re.match(r"^#\s*n\s*=\s*(\d+)", stripped)
        if m_n:
            current_n = int(m_n.group(1))
            continue
        m_data = _NEUMANN_LINE_RE.match(raw)
        if m_data and current_n >= 0:
            re_part = float(m_data.group("re"))
            im_part = float(m_data.group("im"))
            ap = int(m_data.group("ap"))
            a = int(m_data.group("a"))
            q = int(m_data.group("q"))
            key = (ap, a, q)
            sequence.setdefault(key, [])
            # If we ever see gaps in n, we'll just extend with zeros.
            lst = sequence[key]
            while len(lst) <= current_n:
                lst.append(0j)
            lst[current_n] = complex(re_part, im_part)
    return sequence


def pade_approximant(a: List[complex], N: int, M: int, z: complex = 1 + 0j) -> complex:
    """Mirror of pade_approximant() in solve_faddeev.cpp lines 704-729.

    Requires len(a) >= N+M+1. Returns det(P)/det(Q) for the [N/M] Padé.
    """
    if len(a) < N + M + 1:
        return float("nan") + 0j
    # P and Q are (M+1) x (M+1)
    size = M + 1
    P = np.zeros((size, size), dtype=complex)
    Q = np.zeros((size, size), dtype=complex)
    for row in range(M):
        for col in range(M + 1):
            val = a[N - M + 1 + row + col]
            P[row, col] = val
            Q[row, col] = val
    for col in range(M + 1):
        Q[M, col] = z ** (M - col)
        s = 0 + 0j
        for j in range(M - col, N + 1):
            s += a[j - (M - col)] * (z ** j)
        P[M, col] = s
    det_P = np.linalg.det(P)
    det_Q = np.linalg.det(Q)
    return det_P / det_Q


def compute_best_pade(a: List[complex], NM_max: int = 14) -> Tuple[complex, int, bool]:
    """Replicate the solver's final-tail convergence policy for one element.

    Returns (best_PA_value, idx_best_PA, genuinely_converged).
    """
    if not a:
        return 0 + 0j, 0, False
    n_avail = len(a)  # number of Neumann terms stored: a_0..a_{n_avail-1}
    # NM_max feasible from available data:
    NM_feasible_max = min(NM_max, (n_avail - 1) // 2)
    if NM_feasible_max < 0:
        return a[0], 0, False

    PAs: List[complex] = []

    for NM in range(NM_feasible_max + 1):
        # Need a_0..a_{2*NM}
        if 2 * NM + 1 > n_avail:
            break
        try:
            PA = pade_approximant(a, NM, NM, z=1 + 0j)
        except Exception:
            PA = float("nan") + 0j
        PAs.append(PA)

    if not PAs:
        return 0 + 0j, 0, False

    finite_orders = [idx for idx, value in enumerate(PAs)
                     if np.isfinite(value.real) and np.isfinite(value.imag)]
    if not finite_orders:
        return PAs[0], 0, False
    latest_finite_order = finite_orders[-1]

    # Production only certifies convergence at the configured final order and
    # only if its last three consecutive updates satisfy abs+relative tolerance.
    genuinely_converged = (
        NM_feasible_max == NM_max
        and latest_finite_order == NM_max
        and NM_max >= 3
    )
    if genuinely_converged:
        for upper in range(NM_max, NM_max - 3, -1):
            current, previous = PAs[upper], PAs[upper - 1]
            if not (np.isfinite(current.real) and np.isfinite(current.imag)
                    and np.isfinite(previous.real) and np.isfinite(previous.imag)):
                genuinely_converged = False
                break
            tolerance = 1e-7 + 1e-5 * max(abs(current), abs(previous))
            if abs(current - previous) > tolerance:
                genuinely_converged = False
                break

    if genuinely_converged:
        idx_best_PA = NM_max
    else:
        finite_pairs = [upper for upper in range(1, len(PAs))
                        if (np.isfinite(PAs[upper].real)
                            and np.isfinite(PAs[upper].imag)
                            and np.isfinite(PAs[upper - 1].real)
                            and np.isfinite(PAs[upper - 1].imag))]
        idx_best_PA = (min(finite_pairs,
                           key=lambda upper: abs(PAs[upper] - PAs[upper - 1]))
                       if finite_pairs else latest_finite_order)

    return PAs[idx_best_PA], idx_best_PA, genuinely_converged


def reconstruct_u_from_neumann(
    neumann_dir: Path,
    template_path: Path,
    output_dir: Path,
) -> Path:
    """Reconstruct one U_PW file from neumann data + 2NF template."""
    template_name = template_path.name
    # template_name looks like: U_PW_elements_Np_20_Nq_20_JP_5_1_Jmax_3_PSI_0.txt
    m = re.search(r"JP_(\d+)_(\-?\d+)_Jmax_(\d+)", template_name)
    if not m:
        raise ValueError(f"Cannot parse JP from {template_name}")
    two_J_str, parity_str, J2max = m.group(1), m.group(2), m.group(3)
    parity_int = int(parity_str)
    if parity_int not in (1, -1):
        # file may have "_-1" or "_1" but also possibly "_m1" -- handle here.
        raise ValueError(f"parity parse: {parity_str!r}")
    neumann_glob = f"neumann_terms_Np_*_Nq_*_JP_{two_J_str}_{parity_str}_Jmax_{J2max}.txt"
    candidates = list(neumann_dir.glob(neumann_glob))
    if not candidates:
        raise FileNotFoundError(
            f"No neumann file matching {neumann_glob} in {neumann_dir}")
    neumann_path = candidates[0]

    print(f"  neumann: {neumann_path}", file=sys.stderr)
    print(f"  template: {template_path}", file=sys.stderr)

    header_lines, chan_tuples = parse_u_template(template_path)
    n_chan = len(chan_tuples)
    print(f"  channels (rows): {n_chan}", file=sys.stderr)

    neumann_data = parse_neumann_file(neumann_path)
    print(f"  neumann on-shell tuples: {len(neumann_data)}", file=sys.stderr)

    # Discover alpha indices appearing in the neumann file (should be n_chan of them).
    alpha_indices = sorted({ap for (ap, _, _) in neumann_data.keys()})
    print(f"  alpha'-idx values: {alpha_indices}", file=sys.stderr)
    if len(alpha_indices) != n_chan:
        print(
            f"WARN: template has {n_chan} channels but neumann has "
            f"{len(alpha_indices)} alphas. Will match by order.",
            file=sys.stderr,
        )
    # Build mapping row_idx -> alpha_idx (assume ascending order).
    row_to_alpha = {r: a for r, a in enumerate(alpha_indices)}

    # Collect all q_idx values that appear.
    q_indices = sorted({q for (_, _, q) in neumann_data.keys()})
    print(f"  q_idx values: {q_indices}", file=sys.stderr)

    # Read Tlab/Ecm per q_idx from the template (data rows).
    text = template_path.read_text()
    q_to_tlab_ecm: Dict[int, Tuple[float, float]] = {}
    for raw in text.splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        if _HEADER_ROW_RE.match(raw):
            continue
        parts = s.split()
        if len(parts) < 3:
            continue
        try:
            tlab = float(parts[0]); ecm = float(parts[1]); q_idx = int(parts[2])
        except ValueError:
            continue
        q_to_tlab_ecm[q_idx] = (tlab, ecm)

    # Compute best-PA U values per (alpha'_idx, alpha_idx, q_idx)
    U: Dict[Tuple[int, int, int], complex] = {}
    convergence_stats = {"genuinely_converged": 0, "maxiter_truncated": 0,
                         "low_data": 0, "total": 0}
    for key, a_seq in neumann_data.items():
        best_PA, idx_best, conv = compute_best_pade(a_seq, NM_max=14)
        U[key] = best_PA
        convergence_stats["total"] += 1
        if len(a_seq) < 5:
            convergence_stats["low_data"] += 1
        elif conv:
            convergence_stats["genuinely_converged"] += 1
        else:
            convergence_stats["maxiter_truncated"] += 1
    print(
        f"  convergence: {convergence_stats['genuinely_converged']} converged, "
        f"{convergence_stats['maxiter_truncated']} maxiter, "
        f"{convergence_stats['low_data']} low-data",
        file=sys.stderr,
    )

    # Write U_PW file using template header, but data rows filled from U values.
    output_dir.mkdir(parents=True, exist_ok=True)
    out_path = output_dir / template_path.name.replace(
        "_PSI_0", "_PSI_0_reconstructed"
    )
    with out_path.open("w") as f:
        # Copy header lines verbatim up to and including the column header row.
        for raw in header_lines:
            f.write(raw + "\n")
        # Data rows in the same q_idx order as template
        for q_idx in q_indices:
            if q_idx not in q_to_tlab_ecm:
                continue
            tlab, ecm = q_to_tlab_ecm[q_idx]
            tokens: List[str] = [f"{tlab:.16e}", f"{ecm:.16e}", f"{q_idx:>4d}"]
            for r in range(n_chan):
                for c in range(n_chan):
                    alpha_p = row_to_alpha.get(r)
                    alpha = row_to_alpha.get(c)
                    if alpha_p is None or alpha is None:
                        tokens.append(f"{0.0:+.16e}{0.0:+.16e}j")
                        continue
                    val = U.get((alpha_p, alpha, q_idx), 0 + 0j)
                    re_s = f"{val.real:+.16e}"
                    im_s = f"{val.imag:+.16e}"
                    tokens.append(f"{re_s}{im_s}j")
            f.write("   " + "   ".join(tokens) + "\n")
        f.write("# ################################################################################################### \n")
    print(f"  wrote: {out_path}", file=sys.stderr)
    return out_path


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--neumann-dir", required=True, type=Path)
    ap.add_argument("--template", required=True, type=Path)
    ap.add_argument("--output-dir", required=True, type=Path)
    args = ap.parse_args()
    reconstruct_u_from_neumann(args.neumann_dir, args.template, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
