#!/usr/bin/env python3
"""
Partial-wave assembly of the d+p elastic scattering amplitude M(theta) and the
spin-1 polarization observables (dSigma/dOmega, iT11, T20, T21, T22) from the
solver-emitted U_PW_elements_*.txt files.

Physics references (also documented in the in-tree treatise):
  - Ch.12 of docs/treatise (eq. 12-M-pw-summation):
        M(theta; m_d', m_p'; m_d, m_p)
            = -(2*pi)^2 * mu_{p,d} / q_on
              * sum_{J,pi} sum_{(l',j'),(l,j)} U^{Jpi}_{(l',j'),(l,j)} * G_{...}
        where G is the partial-wave geometric factor (jj-coupling form).
  - Ch.13 of docs/treatise (eq. 13-Tkq-master):
        T_{kq}(theta) = Tr[ M . tau_{kq}^(d) . M^dagger ] / Tr[ M . M^dagger ]
        with tau_{1,+1} = -(1/sqrt 2) S_+,
             tau_{2,0}  = (1/sqrt 6)(3 S_z^2 - 2 I),
             tau_{2,+1} = -(1/2)(S_z S_+ + S_+ S_z),
             tau_{2,+2} = (1/2) S_+^2.
  - jj-coupling channel index: alpha = (l, j) where j = l (x) 1/2 is the
    spectator total angular momentum; the deuteron internal projection
    (3S_1 + 3D_1) is already absorbed into U at the on-shell q bin.

This module is a stand-alone replacement for the heuristic 'reduced-U'
formulas in run_dpol_p_observables.py / compare_Ay_experiment.py. It depends
only on numpy + the standard library; no scipy/sympy required.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from solver_u_file_utils import parse_u_file_meta


# ---------------------------------------------------------------------------
# Physical constants (mirror include/constants.h)
# ---------------------------------------------------------------------------

M_PROTON_MEV = 938.272
M_NEUTRON_MEV = 939.565
E_DEUTERON_BINDING_MEV = 2.22457
M_DEUTERON_MEV = M_PROTON_MEV + M_NEUTRON_MEV - E_DEUTERON_BINDING_MEV
MU_PD_MEV = (M_PROTON_MEV * M_DEUTERON_MEV) / (M_PROTON_MEV + M_DEUTERON_MEV)
HBARC_MEV_FM = 197.327


# ---------------------------------------------------------------------------
# Clebsch-Gordan coefficients (Condon-Shortley) and Y_{l,m}(theta, phi=0)
# ---------------------------------------------------------------------------

def _is_half_integer(x: float, tol: float = 1e-9) -> bool:
    return abs(2.0 * x - round(2.0 * x)) < tol


def clebsch_gordan(j1: float, m1: float, j2: float, m2: float, j: float, m: float) -> float:
    """Condon-Shortley Clebsch-Gordan <j1 m1 j2 m2 | j m>.

    All arguments may be integer or half-integer (passed as floats).
    Returns 0 outside the triangle / m-sum / |m| <= j region.
    """
    if abs(m1 + m2 - m) > 1e-9:
        return 0.0
    if j < abs(j1 - j2) - 1e-9 or j > j1 + j2 + 1e-9:
        return 0.0
    if abs(m1) > j1 + 1e-9 or abs(m2) > j2 + 1e-9 or abs(m) > j + 1e-9:
        return 0.0
    for value in (j1 - m1, j1 + m1, j2 - m2, j2 + m2, j - m, j + m,
                  j1 + j2 - j, j1 - j2 + j, -j1 + j2 + j):
        if not _is_half_integer(value):
            return 0.0
        if value < -1e-9:
            return 0.0

    def ifact(x: float) -> int:
        return math.factorial(int(round(x)))

    n_a = ifact(j1 + j2 - j)
    n_b = ifact(j1 - j2 + j)
    n_c = ifact(-j1 + j2 + j)
    n_d = ifact(j1 + j2 + j + 1)
    triangle_sq = (n_a * n_b * n_c) / n_d

    norm_sq = (
        (2.0 * j + 1.0)
        * ifact(j + m) * ifact(j - m)
        * ifact(j1 - m1) * ifact(j1 + m1)
        * ifact(j2 - m2) * ifact(j2 + m2)
    )

    k_min = max(0,
                int(round(j2 - j - m1)),
                int(round(j1 - j + m2)))
    k_max = min(int(round(j1 + j2 - j)),
                int(round(j1 - m1)),
                int(round(j2 + m2)))

    summation = 0.0
    for k in range(k_min, k_max + 1):
        denom = (
            math.factorial(k)
            * ifact(j1 + j2 - j - k)
            * ifact(j1 - m1 - k)
            * ifact(j2 + m2 - k)
            * ifact(j - j2 + m1 + k)
            * ifact(j - j1 - m2 + k)
        )
        summation += ((-1) ** k) / denom

    return math.sqrt(triangle_sq * norm_sq) * summation


def _associated_legendre(l: int, m: int, x: float) -> float:
    """Condon-Shortley associated Legendre P_l^m(x) for m >= 0, |x| <= 1."""
    if m < 0 or m > l:
        return 0.0
    if abs(x) > 1.0 + 1e-12:
        raise ValueError(f"|x| > 1 in associated Legendre: x={x}")
    x = max(-1.0, min(1.0, x))

    pmm = 1.0
    if m > 0:
        somx2 = math.sqrt(max(0.0, 1.0 - x * x))
        fact = 1.0
        for _ in range(1, m + 1):
            pmm *= -fact * somx2
            fact += 2.0
    if l == m:
        return pmm

    pmmp1 = x * (2 * m + 1) * pmm
    if l == m + 1:
        return pmmp1

    pll = 0.0
    for ll in range(m + 2, l + 1):
        pll = ((2 * ll - 1) * x * pmmp1 - (ll + m - 1) * pmm) / (ll - m)
        pmm = pmmp1
        pmmp1 = pll
    return pll


def spherical_harmonic_phi0(l: int, m: int, theta: float) -> float:
    """Y_{l,m}(theta, phi=0) using Condon-Shortley convention (real-valued).

    Y_{l,m}(theta, 0) = sqrt[(2l+1)/(4 pi) * (l-m)!/(l+m)!] * P_l^m(cos theta).
    For m < 0: Y_{l,-|m|}(theta, 0) = (-1)^|m| Y_{l,|m|}(theta, 0).
    """
    if l < 0 or abs(m) > l:
        return 0.0
    cos_t = math.cos(theta)
    abs_m = abs(m)
    norm_sq = (2 * l + 1) / (4.0 * math.pi)
    norm_sq *= math.factorial(l - abs_m) / math.factorial(l + abs_m)
    value = math.sqrt(norm_sq) * _associated_legendre(l, abs_m, cos_t)
    if m < 0:
        value *= (-1) ** abs_m
    return value


# ---------------------------------------------------------------------------
# Spin-1 spherical tensor matrices (basis: |1,+1>, |1,0>, |1,-1>)
# ---------------------------------------------------------------------------

_SQRT2 = math.sqrt(2.0)
_SQRT6 = math.sqrt(6.0)

# S_z, S_+ for spin-1 in basis ordered (|+1>, |0>, |-1>):
_S_PLUS = np.array([[0.0, _SQRT2, 0.0],
                    [0.0, 0.0, _SQRT2],
                    [0.0, 0.0, 0.0]], dtype=complex)
_S_MINUS = _S_PLUS.conj().T
_S_Z = np.diag([1.0, 0.0, -1.0]).astype(complex)
_I3 = np.eye(3, dtype=complex)

TAU_1_PLUS_1 = -(1.0 / _SQRT2) * _S_PLUS
TAU_2_0 = (1.0 / _SQRT6) * (3.0 * (_S_Z @ _S_Z) - 2.0 * _I3)
TAU_2_PLUS_1 = -0.5 * (_S_Z @ _S_PLUS + _S_PLUS @ _S_Z)
TAU_2_PLUS_2 = 0.5 * (_S_PLUS @ _S_PLUS)


# ---------------------------------------------------------------------------
# Spin-state ordering for the 6x6 M matrix
# ---------------------------------------------------------------------------

# Index ordering: i = m_d_idx * 2 + m_p_idx
#   m_d_idx: 0,1,2 for m_d = +1, 0, -1
#   m_p_idx: 0,1   for m_p = +1/2, -1/2
M_D_VALUES = (1.0, 0.0, -1.0)
M_P_VALUES = (0.5, -0.5)

SPIN_INDEX_TABLE: Tuple[Tuple[int, int, float, float], ...] = tuple(
    (m_d_idx, m_p_idx, m_d, m_p)
    for m_d_idx, m_d in enumerate(M_D_VALUES)
    for m_p_idx, m_p in enumerate(M_P_VALUES)
)  # length 6


def embed_deuteron_operator(op3: np.ndarray) -> np.ndarray:
    """Embed a 3x3 deuteron operator into the 6x6 (deuteron x proton) space."""
    if op3.shape != (3, 3):
        raise ValueError("operator must be 3x3 in deuteron spin space")
    return np.kron(op3, np.eye(2, dtype=complex))


# ---------------------------------------------------------------------------
# U-file parsing
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class UChannel:
    """A spectator (l, j) channel within a single (J, pi) block."""
    l: int        # spectator orbital angular momentum
    two_j: int    # 2 * spectator total angular momentum

    @property
    def j(self) -> float:
        return 0.5 * self.two_j


@dataclass
class UEnergyPoint:
    """U-matrix at a single on-shell q bin within a (J, pi) block."""
    tlab: float          # MeV
    ecm: float           # MeV
    q_idx: int
    matrix: np.ndarray   # complex (n_chan, n_chan), units MeV


@dataclass
class JPiBlock:
    """All on-shell U-matrix data for one (J, pi) partial wave."""
    two_J: int           # 2 * J_3N
    parity: int          # +1 or -1
    channels: List[UChannel]   # row/col indexing
    points: List[UEnergyPoint]
    source_path: Path

    @property
    def J(self) -> float:
        return 0.5 * self.two_J


_HEADER_ROW_RE = re.compile(
    r"^\s*U(?P<i>\d+)(?P<j>\d+)\s+(?P<row>\d+)\s+(?P<col>\d+)\s+"
    r"(?P<lp>\d+)\s+(?P<two_jp>\d+)\s+(?P<l>\d+)\s+(?P<two_j>\d+)\s*$"
)


def _parse_complex_token(token: str) -> complex:
    # File format examples: '+4.20e-01-4.22e-01j', '-1.21e+00+5.4e-02j'.
    return complex(token.replace(" ", ""))


def parse_u_file(path: Path) -> Optional[JPiBlock]:
    """Parse one U_PW_elements_*.txt file into a structured JPiBlock."""
    meta = parse_u_file_meta(path)
    if meta is None:
        return None

    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()

    channels_by_idx: Dict[int, UChannel] = {}
    data_rows: List[Tuple[float, float, int, List[complex]]] = []

    for raw in lines:
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            continue
        match_chan = _HEADER_ROW_RE.match(raw)
        if match_chan is not None:
            i = int(match_chan.group("i"))
            j = int(match_chan.group("j"))
            if i == j:
                channels_by_idx.setdefault(
                    i,
                    UChannel(
                        l=int(match_chan.group("lp")),
                        two_j=int(match_chan.group("two_jp")),
                    ),
                )
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        try:
            tlab = float(parts[0])
            ecm = float(parts[1])
            q_idx = int(parts[2])
        except ValueError:
            continue
        try:
            entries = [_parse_complex_token(tok) for tok in parts[3:]]
        except ValueError:
            continue
        data_rows.append((tlab, ecm, q_idx, entries))

    if not channels_by_idx or not data_rows:
        return None

    n_chan = max(channels_by_idx) + 1
    if any(idx not in channels_by_idx for idx in range(n_chan)):
        return None
    channels = [channels_by_idx[idx] for idx in range(n_chan)]

    expected_entries = n_chan * n_chan
    points: List[UEnergyPoint] = []
    for tlab, ecm, q_idx, entries in data_rows:
        if len(entries) < expected_entries:
            continue
        flat = entries[:expected_entries]
        matrix = np.array(flat, dtype=complex).reshape(n_chan, n_chan)
        points.append(UEnergyPoint(tlab=tlab, ecm=ecm, q_idx=q_idx, matrix=matrix))

    if not points:
        return None

    return JPiBlock(
        two_J=meta.two_j,
        parity=meta.parity,
        channels=channels,
        points=points,
        source_path=path,
    )


def parse_solver_output(solver_out_dir: Path, u_files: Optional[Sequence[Path]] = None) -> List[JPiBlock]:
    if u_files is None:
        u_files = sorted(solver_out_dir.glob("U_PW_elements_*.txt"))
    blocks: List[JPiBlock] = []
    for path in u_files:
        block = parse_u_file(path)
        if block is not None:
            blocks.append(block)
    return blocks


# ---------------------------------------------------------------------------
# Available solver energies (intersection across blocks at the same q-idx)
# ---------------------------------------------------------------------------

def list_solver_energies(blocks: Sequence[JPiBlock], tol_mev: float = 1e-3) -> List[Tuple[float, float]]:
    """Return [(tlab, ecm)] energies that appear in every (J, pi) block.

    Energies are matched by q_idx, since the on-shell mapping is identical
    across blocks (same q-WP grid). Tlab/Ecm are read from the first block
    that has the q_idx and verified consistent across the rest within tol.
    """
    if not blocks:
        return []
    common: Optional[set] = None
    per_block_qidx: List[Dict[int, UEnergyPoint]] = []
    for block in blocks:
        qmap = {pt.q_idx: pt for pt in block.points}
        per_block_qidx.append(qmap)
        if common is None:
            common = set(qmap)
        else:
            common &= set(qmap)
    if not common:
        return []
    out: List[Tuple[float, float]] = []
    for q_idx in sorted(common):
        ref_pt = per_block_qidx[0][q_idx]
        consistent = all(
            abs(qmap[q_idx].tlab - ref_pt.tlab) < tol_mev
            and abs(qmap[q_idx].ecm - ref_pt.ecm) < tol_mev
            for qmap in per_block_qidx[1:]
        )
        if consistent:
            out.append((ref_pt.tlab, ref_pt.ecm))
    return out


def select_energy_q_idx(blocks: Sequence[JPiBlock], target_tlab_mev: float) -> Tuple[int, float, float]:
    """Pick the q_idx common to all blocks whose Tlab is closest to the target."""
    energies = list_solver_energies(blocks)
    if not energies:
        raise RuntimeError("No common solver energy across (J, pi) blocks")
    best = min(energies, key=lambda item: abs(item[0] - target_tlab_mev))
    target_tlab, target_ecm = best
    block_qmap = {pt.q_idx: pt for pt in blocks[0].points}
    chosen_q_idx: Optional[int] = None
    for q_idx, pt in block_qmap.items():
        if abs(pt.tlab - target_tlab) < 1e-6 and abs(pt.ecm - target_ecm) < 1e-6:
            chosen_q_idx = q_idx
            break
    if chosen_q_idx is None:
        raise RuntimeError("Failed to locate q_idx for chosen energy")
    return chosen_q_idx, target_tlab, target_ecm


# ---------------------------------------------------------------------------
# Geometric factor and M-matrix assembly (jj coupling)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class WPGridBin:
    """One on-shell q-WP bin: midpoint q and bin width dq, both in MeV/c."""
    q_mid_mev: float
    dq_mev: float
    tlab_mid_mev: float


def parse_q_kinematics(solver_out_dir: Path) -> Dict[int, WPGridBin]:
    """Parse q_kinematics_Nq_*.txt, returning q-bin midpoints and widths.

    The file contains two tables: bin boundaries first, bin midpoints second.
    Both are indexed by q-WP bin index. dq is computed from boundary[i+1] -
    boundary[i] for bin i.
    """
    matches = sorted(solver_out_dir.glob("q_kinematics_Nq_*.txt"))
    if not matches:
        raise FileNotFoundError(
            f"No q_kinematics_Nq_*.txt found in {solver_out_dir}; cannot derive WP grid"
        )
    path = matches[-1]
    text = path.read_text(encoding="utf-8")
    lines = text.splitlines()

    sections: List[List[Tuple[int, float, float, float]]] = [[]]
    for raw in lines:
        stripped = raw.strip()
        if not stripped:
            continue
        if stripped.startswith("#"):
            if "MID-POINTS" in stripped.upper() and sections[-1]:
                sections.append([])
            continue
        parts = stripped.split()
        if len(parts) < 4:
            continue
        try:
            idx = int(parts[0])
            q_mev = float(parts[1])
            ecm = float(parts[2])
            tlab = float(parts[3])
        except ValueError:
            continue
        sections[-1].append((idx, q_mev, ecm, tlab))

    if len(sections) < 2 or not sections[0] or not sections[1]:
        raise RuntimeError(f"Unexpected q_kinematics format in {path}")

    boundaries = {idx: q for idx, q, _, _ in sections[0]}
    midpoints = sections[1]

    out: Dict[int, WPGridBin] = {}
    for idx, q_mid, _, tlab in midpoints:
        q_lo = boundaries.get(idx)
        q_hi = boundaries.get(idx + 1)
        if q_lo is None or q_hi is None:
            continue
        out[idx] = WPGridBin(q_mid_mev=q_mid, dq_mev=q_hi - q_lo, tlab_mid_mev=tlab)
    return out


def kinematic_prefactor_wp(bin_info: WPGridBin) -> float:
    """Witala/Gloeckle prefactor mapping U_WP[MeV] -> M[fm] in WP basis.

    M(theta) = -(2 pi)^2 * mu_{p,d} / [(hbar c)^2 * q^2 * dq] * U_WP,
    where U_WP is the on-shell wave-packet matrix element in MeV and q, dq
    are in fm^{-1} (so the q-bin volume q^2 * dq has units fm^{-3}; combined
    with (hbar c)^2 [MeV^2 fm^2] this yields the proper fm/MeV factor).
    """
    if bin_info.q_mid_mev <= 0.0 or bin_info.dq_mev <= 0.0:
        return 0.0
    q_inv_fm = bin_info.q_mid_mev / HBARC_MEV_FM
    dq_inv_fm = bin_info.dq_mev / HBARC_MEV_FM
    return -((2.0 * math.pi) ** 2) * MU_PD_MEV / (
        (HBARC_MEV_FM ** 2) * (q_inv_fm ** 2) * dq_inv_fm
    )


def _geometric_factor(
    block: JPiBlock,
    chan_row: UChannel,
    chan_col: UChannel,
    theta: float,
    m_d_prime: float,
    m_p_prime: float,
    m_d: float,
    m_p: float,
) -> float:
    """jj-coupling geometric factor for one (alpha', alpha) within a (J, pi) block."""
    j_total = block.J
    M_total = m_d + m_p
    if abs(M_total) > j_total + 1e-9:
        return 0.0

    l_in = chan_col.l
    j_in = chan_col.j
    l_out = chan_row.l
    j_out = chan_row.j

    # Parity check (deuteron is positive parity; 3N parity = (-)^l)
    if (-1) ** l_in != block.parity or (-1) ** l_out != block.parity:
        return 0.0

    cg_in_orb = clebsch_gordan(l_in, 0.0, 0.5, m_p, j_in, m_p)
    if cg_in_orb == 0.0:
        return 0.0
    cg_in_jjd = clebsch_gordan(j_in, m_p, 1.0, m_d, j_total, M_total)
    if cg_in_jjd == 0.0:
        return 0.0

    m_j_prime = M_total - m_d_prime
    m_l_prime = m_j_prime - m_p_prime
    if abs(m_j_prime) > j_out + 1e-9:
        return 0.0
    if abs(m_l_prime) > l_out + 1e-9:
        return 0.0
    cg_out_orb = clebsch_gordan(l_out, m_l_prime, 0.5, m_p_prime, j_out, m_j_prime)
    if cg_out_orb == 0.0:
        return 0.0
    cg_out_jjd = clebsch_gordan(j_out, m_j_prime, 1.0, m_d_prime, j_total, M_total)
    if cg_out_jjd == 0.0:
        return 0.0

    y_in_amp = math.sqrt((2 * l_in + 1) / (4.0 * math.pi))
    y_out = spherical_harmonic_phi0(l_out, int(round(m_l_prime)), theta)
    if y_out == 0.0:
        return 0.0

    return y_in_amp * cg_in_orb * cg_in_jjd * y_out * cg_out_orb * cg_out_jjd


def assemble_m_matrix(
    blocks: Sequence[JPiBlock],
    q_idx: int,
    theta: float,
    *,
    bin_info: Optional[WPGridBin] = None,
    extra_scale: float = 1.0,
) -> np.ndarray:
    """Build the 6x6 elastic d+p amplitude M(theta) at the chosen on-shell bin.

    Index ordering (rows = final, cols = initial) follows SPIN_INDEX_TABLE.
    If ``bin_info`` is provided, the WP-basis kinematic prefactor is applied
    so that |M|^2 has units fm^2; otherwise M is returned in raw partial-wave
    units (still useful for polarization observables, which are ratios).
    The ``extra_scale`` knob lets callers absorb any residual calibration
    constant determined empirically against experiment.
    """
    M = np.zeros((6, 6), dtype=complex)

    for block in blocks:
        point: Optional[UEnergyPoint] = None
        for candidate in block.points:
            if candidate.q_idx == q_idx:
                point = candidate
                break
        if point is None:
            continue

        for row_idx, chan_row in enumerate(block.channels):
            for col_idx, chan_col in enumerate(block.channels):
                u_value = point.matrix[row_idx, col_idx]
                if u_value == 0.0:
                    continue
                for j_final, (_, _, m_d_p, m_p_p) in enumerate(SPIN_INDEX_TABLE):
                    for i_initial, (_, _, m_d, m_p) in enumerate(SPIN_INDEX_TABLE):
                        g = _geometric_factor(
                            block, chan_row, chan_col, theta,
                            m_d_p, m_p_p, m_d, m_p,
                        )
                        if g == 0.0:
                            continue
                        M[j_final, i_initial] += u_value * g

    if bin_info is not None:
        M *= kinematic_prefactor_wp(bin_info)
    if extra_scale != 1.0:
        M *= extra_scale
    return M


# ---------------------------------------------------------------------------
# Observables from M
# ---------------------------------------------------------------------------

@dataclass
class PolarizationObservables:
    dsigma_fm2_per_sr: float    # 1/6 sum |M|^2 in fm^2/sr
    iT11: float
    T20: float
    T21: float
    T22: float
    trace_M_Mdag: float


def observables_from_M(M: np.ndarray) -> PolarizationObservables:
    if M.shape != (6, 6):
        raise ValueError("M must be 6x6")

    M_dag = M.conj().T
    M_Mdag = M @ M_dag
    trace_norm = float(np.trace(M_Mdag).real)
    if trace_norm <= 0.0:
        return PolarizationObservables(0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    # Unpolarized differential cross section: 1/6 * sum_{m_d,m_p,m_d',m_p'} |M|^2
    dsigma = trace_norm / 6.0

    def _A(tau3: np.ndarray) -> complex:
        tau6 = embed_deuteron_operator(tau3)
        return complex(np.trace(M @ tau6 @ M_dag)) / trace_norm

    A_11 = _A(TAU_1_PLUS_1)
    A_20 = _A(TAU_2_0)
    A_21 = _A(TAU_2_PLUS_1)
    A_22 = _A(TAU_2_PLUS_2)

    return PolarizationObservables(
        dsigma_fm2_per_sr=dsigma,
        iT11=float((1j * A_11).real),
        T20=float(A_20.real),
        T21=float(A_21.real),
        T22=float(A_22.real),
        trace_M_Mdag=trace_norm,
    )


def observables_at_angles(
    blocks: Sequence[JPiBlock],
    q_idx: int,
    theta_deg_grid: Sequence[float],
    *,
    bin_info: Optional[WPGridBin] = None,
    extra_scale: float = 1.0,
) -> List[PolarizationObservables]:
    return [
        observables_from_M(
            assemble_m_matrix(
                blocks,
                q_idx,
                math.radians(theta),
                bin_info=bin_info,
                extra_scale=extra_scale,
            )
        )
        for theta in theta_deg_grid
    ]


def calibrate_dsigma_scale(
    model_dsigma: Sequence[float],
    exp_angles_deg: Sequence[float],
    exp_values: Sequence[float],
    model_angles_deg: Sequence[float],
) -> float:
    """Geometric-mean calibration between model and experiment dSigma curves.

    Returns a scalar k such that k * model_dsigma best matches experiment in
    log-space. We do least-squares on log(k) = mean(log(exp) - log(model))
    over experimental angles (matched to nearest model angle). Robust to the
    factor-of-a-few unit ambiguity in the WP-basis kinematic prefactor.
    """
    if not model_dsigma or not exp_values:
        return 1.0
    log_diffs: List[float] = []
    for ang_e, val_e in zip(exp_angles_deg, exp_values):
        if val_e <= 0.0:
            continue
        idx = min(range(len(model_angles_deg)), key=lambda i: abs(model_angles_deg[i] - ang_e))
        m = model_dsigma[idx]
        if m <= 0.0:
            continue
        log_diffs.append(math.log(val_e) - math.log(m))
    if not log_diffs:
        return 1.0
    return math.exp(sum(log_diffs) / len(log_diffs))
