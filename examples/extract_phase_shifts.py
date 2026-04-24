#!/usr/bin/env python3
"""
Extract nd elastic-scattering partial-wave phase shifts from Tic-tac
U_PW_elements_*.txt output. Implements Miller, Ekström, Hebeler
PRC 106, 024001 (2022), Appendix D/E.

Conversion path (Miller eqs. D3, D5):

    U^plane(q0)  =  U^WP  *  f̄²(q0) / (q0² * ΔE_j)
                 =  U^WP  / (q0 * μ1 * ΔE_j)           (energy WPs: f̄² = q/μ1)

    S           =  I  −  2πi · q0 · m_N · U^plane(q0)

    δ_diag      =  ½ arg(S_diag)              (quick diagonal readout)

    S = U exp(2iδ) U^T     (full eigenphase analysis, Blatt–Biedenharn/E2)

For J^Π = 1/2+, the 2-channel doublet contains (l=0, j=1/2) ²S₁/₂ and
(l=2, j=3/2). Miller Fig. 1 at E_lab = 13 MeV (Nijmegen-I) gives the dominant
doublet phase: Re δ ≈ 105°, Im δ ≈ 15°.

Usage:
  python extract_phase_shifts.py --solver-out-dir <dir> [--tlab 13.0]
"""

import argparse
import math
import re
from pathlib import Path

import numpy as np


HBARC_MEV_FM = 197.327
M_PROTON_MEV = 938.272
M_NEUTRON_MEV = 939.565
M_DEUTERON_MEV = 1875.61

# Miller's m_N for the S-matrix prefactor (PRC 106, text after Eq. D3).
MN_MILLER = 2.0 * M_PROTON_MEV * M_NEUTRON_MEV / (M_PROTON_MEV + M_NEUTRON_MEV)

# Reduced mass of the N-d system (spectator-pair system) for f̄²(q) = q/μ1 and
# ΔE_j = (q_{j+1}² - q_j²)/(2 μ1).
MU1_ND = M_NEUTRON_MEV * M_DEUTERON_MEV / (M_NEUTRON_MEV + M_DEUTERON_MEV)


# --------------------------- parsing helpers --------------------------- #

_COMPLEX_RE = re.compile(r"([+\-]\d\.\d+e[+\-]\d+)([+\-]\d\.\d+e[+\-]\d+)j")


def parse_complex_token(tok: str) -> complex:
    m = _COMPLEX_RE.match(tok.strip())
    if m is None:
        raise ValueError(f"bad complex token: {tok!r}")
    return complex(float(m.group(1)), float(m.group(2)))


def parse_u_pw_file(path: Path):
    """Return (channels, rows) where:
        channels = list[dict(idx, l, two_j)]
        rows     = list[dict(tlab, ecm, q_idx, U_matrix (ndim x ndim complex))]
    """
    channels = []
    rows = []
    in_table1 = False
    in_table2 = False
    header_cols = None

    with path.open() as f:
        for line in f:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                if "Name" in stripped and "row-idx" in stripped and "l'" in stripped:
                    in_table1 = True
                    in_table2 = False
                    continue
                if "Tlab" in stripped and "U00" in stripped:
                    in_table1 = False
                    in_table2 = True
                    # Determine how many complex U columns are present.
                    header_cols = [c for c in stripped.strip("#").split() if c.startswith("U")]
                    continue
                if "####" in stripped:
                    in_table1 = False
                    # keep in_table2 as is; the footer "####" after data closes it
                    continue
                continue

            if in_table1:
                # label row-idx col-idx l' 2j' l 2j
                toks = stripped.split()
                if len(toks) != 7:
                    continue
                label, r, c, lp, jp, l, j = toks
                channels.append(
                    dict(
                        label=label,
                        row=int(r), col=int(c),
                        lp=int(lp), two_jp=int(jp),
                        l=int(l), two_j=int(j),
                    )
                )
                continue

            if in_table2:
                toks = stripped.split()
                if len(toks) < 4:
                    continue
                try:
                    tlab = float(toks[0])
                    ecm = float(toks[1])
                    q_idx = int(toks[2])
                except ValueError:
                    continue
                u_tokens = toks[3:]  # header slot count is a 3x3 template; data is n²
                n = int(round(math.sqrt(len(u_tokens))))
                if n * n != len(u_tokens):
                    raise ValueError(
                        f"U token count {len(u_tokens)} is not a perfect square in {path}"
                    )
                U = np.zeros((n, n), dtype=complex)
                for k, tok in enumerate(u_tokens):
                    U[k // n, k % n] = parse_complex_token(tok)
                rows.append(dict(tlab=tlab, ecm=ecm, q_idx=q_idx, U=U))

    # channels table lists every (row, col) pair; dedupe by row-index.
    by_row = {}
    for c in channels:
        if c["row"] == c["col"]:
            by_row[c["row"]] = dict(row=c["row"], l=c["lp"], two_j=c["two_jp"])
    channels = [by_row[k] for k in sorted(by_row)]
    return channels, rows


def parse_q_kinematics(path: Path):
    """Parse q_kinematics_Nq_*.txt. Returns (boundaries, midpoints) as
    lists of dicts with keys {q, ecm, tlab}. Boundaries has Nq+1 entries,
    midpoints has Nq entries.
    """
    boundaries = []
    midpoints = []
    section = None  # None, "boundaries", "midpoints"
    with path.open() as f:
        for line in f:
            s = line.strip()
            if not s:
                continue
            if s.startswith("#"):
                if "BOUNDARIES" in s:
                    section = "boundaries"
                    continue
                if "MID-POINTS" in s:
                    section = "midpoints"
                    continue
                continue
            toks = s.split()
            if len(toks) < 4:
                continue
            try:
                idx = int(toks[0])
                q = float(toks[1])
                ecm = float(toks[2])
                tlab = float(toks[3])
            except ValueError:
                continue
            row = dict(idx=idx, q=q, ecm=ecm, tlab=tlab)
            if section == "boundaries":
                boundaries.append(row)
            elif section == "midpoints":
                midpoints.append(row)
    return boundaries, midpoints


# ----------------------- physics conversion ----------------------- #

def wp_to_plane_wave(U_wp: np.ndarray, q0_mev: float, dE_mev: float) -> np.ndarray:
    """Miller Eq. (D5) for energy WPs:
        U^plane(q0) = U^WP * f̄²(q0) / (q0² ΔE_j)
                    = U^WP / (q0 · μ1 · ΔE_j)
    Units: U^plane has [MeV]⁻².
    """
    return U_wp / (q0_mev * MU1_ND * dE_mev)


def plane_wave_to_S(U_plane: np.ndarray, q0_mev: float) -> np.ndarray:
    """Miller Eq. (D3): S = I − 2πi q0 m_N U^plane (dimensionless)."""
    n = U_plane.shape[0]
    return np.eye(n, dtype=complex) - 2.0j * math.pi * q0_mev * MN_MILLER * U_plane


def eigenphases(S: np.ndarray):
    """Blatt–Biedenharn style eigenphases: S = U exp(2iδ) U^T.
    For complex S above threshold we use symmetric decomposition: Takagi's
    factorization of a complex symmetric matrix. For a quick readout we just
    diagonalize S and take δ_k = ½ arg(λ_k). Mixing-angle reconstruction is
    out of scope for the first sanity check.
    """
    eigvals = np.linalg.eigvals(S)
    deltas = 0.5 * np.angle(eigvals)  # radians, principal branch
    # Also return |λ| as "inelasticity" diagnostic (|λ|=1 means pure elastic).
    return eigvals, deltas


def rad_to_deg(x: complex | float) -> float:
    return float(x) * 180.0 / math.pi


def canonical_delta_deg(delta_complex: complex) -> tuple[float, float]:
    """Map (Re δ, Im δ) into the Miller plot convention:
    Re δ ∈ (0°, 180°] by using δ ≡ δ+180° (same S), pick branch with positive Re.
    Returns (Re_deg, Im_deg).
    """
    re = delta_complex.real * 180.0 / math.pi
    im = delta_complex.imag * 180.0 / math.pi
    # shift Re into [0, 180)
    re = re - 180.0 * math.floor(re / 180.0)
    return re, im


# --------------------------- driver --------------------------- #

def extract_for_file(u_path: Path, q_kin_path: Path, tlab_target: float | None,
                     verbose: bool = True):
    channels, rows = parse_u_pw_file(u_path)
    boundaries, midpoints = parse_q_kinematics(q_kin_path)

    if not rows:
        print(f"[warn] no U rows in {u_path.name}")
        return

    # header metadata
    header_text = u_path.read_text().splitlines()[:12]
    meta = {}
    for ln in header_text:
        if "JP:" in ln:
            meta["JP"] = ln.split(":", 1)[1].strip()
        if "Potential:" in ln:
            meta["Potential"] = ln.split(":", 1)[1].strip()

    print(f"\n=== {u_path.name} ===")
    for k, v in meta.items():
        print(f"  {k}: {v}")
    print(f"  channels (row-idx : l,2j): "
          + ", ".join(f"{c['row']}:{c['l']},{c['two_j']}" for c in channels))

    # For each Tlab row, extract phase shifts
    results = []
    for r in rows:
        q_idx = r["q_idx"]
        q0 = midpoints[q_idx]["q"]                       # MeV
        # Bin width in energy (ΔE_j = E_{j+1}^bnd − E_j^bnd for energy WP)
        # boundaries has Nq+1 entries indexed 0..Nq; bin j spans [j, j+1]
        dE = boundaries[q_idx + 1]["ecm"] - boundaries[q_idx]["ecm"]

        U_wp = r["U"]
        # Be defensive: keep only the N=len(channels) sub-block
        n_ch = len(channels)
        U_wp_chn = U_wp[:n_ch, :n_ch]

        U_plane = wp_to_plane_wave(U_wp_chn, q0, dE)
        S = plane_wave_to_S(U_plane, q0)
        eigs, deltas = eigenphases(S)

        # Also directly extract δ from diagonal S (no mixing)
        delta_diag = 0.5 * np.angle(np.diag(S))
        inelast_diag = np.abs(np.diag(S))

        results.append(dict(
            tlab=r["tlab"], ecm=r["ecm"], q0=q0, dE=dE,
            S=S, eigvals=eigs, delta_eig=deltas,
            delta_diag=delta_diag, inelast_diag=inelast_diag,
        ))

    # Pretty-print selected Tlab
    for res in results:
        if tlab_target is not None and abs(res["tlab"] - tlab_target) > 5.0:
            continue
        print(f"\n  Tlab = {res['tlab']:.3f} MeV  |  Ecm = {res['ecm']:.3f} MeV  "
              f"|  q0 = {res['q0']:.3f} MeV  |  ΔE_j = {res['dE']:.3f} MeV")
        # Debug: print S and asymmetry / unitarity diagnostics
        S = res["S"]
        asym = np.linalg.norm(S - S.T)
        unit_def = np.linalg.norm(S @ S.conj().T - np.eye(S.shape[0]))
        print(f"    ||S - S^T|| = {asym:.4e}   ||SS†-1|| = {unit_def:.4e}")
        # Diagonal readout (no recoupling / no eigenphase mix)
        for k in range(len(res["delta_diag"])):
            d = res["delta_diag"][k]
            lab = channels[k]
            re, im = canonical_delta_deg(d)
            print(f"    δ_diag [row={k}, (l,2j)=({lab['l']},{lab['two_j']})]  "
                  f"= {re:+8.3f}° + i·{im:+7.3f}°   "
                  f"|S_kk| = {res['inelast_diag'][k]:.4f}")
        # Eigenphase (2x2 Blatt–Biedenharn — quick)
        for k, (ev, de) in enumerate(zip(res["eigvals"], res["delta_eig"])):
            re, im = canonical_delta_deg(de)
            print(f"    δ_eig  [k={k}]          "
                  f"= {re:+8.3f}° + i·{im:+7.3f}°   "
                  f"|λ_k|   = {abs(ev):.4f}")

    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solver-out-dir", required=True, type=Path,
                    help="directory holding U_PW_elements_*.txt + q_kinematics_Nq_*.txt")
    ap.add_argument("--tlab", type=float, default=None,
                    help="show only rows within ±5 MeV of this Tlab (default: all)")
    args = ap.parse_args()

    u_files = sorted(args.solver_out_dir.glob("U_PW_elements_*.txt"))
    if not u_files:
        raise SystemExit(f"no U_PW_elements_*.txt under {args.solver_out_dir}")

    # Pick matching q_kinematics file by Nq suffix.
    for u_path in u_files:
        m = re.search(r"Nq_(\d+)", u_path.name)
        if not m:
            continue
        nq = m.group(1)
        q_kin = args.solver_out_dir / f"q_kinematics_Nq_{nq}.txt"
        if not q_kin.exists():
            print(f"[skip] no q_kinematics_Nq_{nq}.txt for {u_path.name}")
            continue
        extract_for_file(u_path, q_kin, args.tlab)


if __name__ == "__main__":
    main()
