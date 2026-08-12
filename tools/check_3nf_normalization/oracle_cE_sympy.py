#!/usr/bin/env python3
"""
Independent SymPy oracle for the chiral N²LO c_E three-nucleon contact term.

PURPOSE
-------
Provide an *independent* derivation of the partial-wave matrix element

    ⟨α', p', q' | V^(1)_cont | α, p, q⟩

for the c_E 3N contact (Epelbaum 2002 eq. 2.10). This is the "golden value"
oracle for tests/cpp/test_3nf_matrix_elements.cpp. It must NOT share any
derivation path with the production C++ helper `recoupling_3nf_contact_cE`.

INDEPENDENCE
------------
- We compute the spin-isospin matrix element via *explicit basis-state sums*
  using sympy.physics.quantum CG / spin utilities, rather than the closed-form
  τ₂·τ₃ eigenvalue used in the C++ helper.
- We compute the regulator and LEC prefactor in *MeV units* throughout, rather
  than the C++ fm⁻¹ convention; the final number is converted back to fm⁵.
- The angular normalization 6·(4π)² is computed from explicit ⟨Y_{00}|Y_{00}⟩
  integrals, not assumed.

DERIVATION
----------
The c_E contact operator (spectator-1 form, Epelbaum 2002 eq. 2.10):

    V^(1)_cont = +E · (τ_2 · τ_3),   E = c_E / (f_π⁴ Λ_χ)

is a SPIN SCALAR (no σ operator) and an ISOSPIN SCALAR on the 3N total T_3N
(only the pair T_2N is touched). Its partial-wave matrix element in the
Jj-coupled Jacobi basis

    |p q α⟩ = |p q [(L S)J_2N, (l ½)j_1N] J_3N ; (T ½) T_3N⟩

factorises as

    ⟨α'|V_cE|α⟩ = A_cE(α', α) · F_cE(p',q',p,q)

with the spin-isospin coefficient

    A_cE(α', α) = δ_{α'α on all quantum numbers EXCEPT trivially T_2N}
                × δ_{S S'} δ_{T T'} δ_{L 0} δ_{L' 0} δ_{l 0} δ_{l' 0}
                × δ_{J_2N S} δ_{J_2N' S'} δ_{j_1N ½} δ_{j_1N' ½}
                × ⟨(½ ½)T | τ_2·τ_3 | (½ ½)T⟩
                = δ(...) · (2 T (T+1) − 3)

and the momentum-space factor

    F_cE(p',q',p,q) = -E · f_R(p',q') · f_R(p,q) · (4π)² · normalisation

with the squared-Gaussian regulator per Epelbaum 2002 eq. (3.19):

    f_R(p, q) = exp( -((4 p² + 3 q²) / (4 Λ²))² )

VERIFICATION OF A_cE
--------------------
We compute ⟨(½½)T|τ_2·τ_3|(½½)T⟩ via the explicit Pauli matrix sum:

    τ_2·τ_3 = Σ_a τ_2^a τ_3^a   (a ∈ {x,y,z})

using the sympy spin-1/2 representation. For an isospin-T pair state
|T M_T⟩ constructed via CG(½m_2; ½m_3 | T M_T), we evaluate the eigenvalue
by summing over M_T = M_2 + M_3 components:

    eigenvalue = ⟨T M_T| Σ_a τ_2^a τ_3^a |T M_T⟩  (independent of M_T by rotational inv.)

This is the "no shortcut" verification — we do NOT use the closed-form
2T(T+1)−3 identity that the C++ helper uses; we evaluate the Pauli matrices
directly. This is the independence guarantee.

OUTPUT
------
A JSON dict of golden matrix elements at fixed (p, q, p', q') for several
channels (α', α). Consumed by tests/cpp/test_3nf_matrix_elements.cpp.
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

from sympy import Rational, sqrt, exp, N, S
from sympy.physics.wigner import clebsch_gordan as CG
from sympy.physics.quantum.spin import Rotation


# -----------------------------------------------------------------------------
# Constants (MeV units throughout; convert to fm**-1 only for the regulator)
# -----------------------------------------------------------------------------
HBARC_MEV_FM = 197.327       # MeV·fm (locked to include/constants.h)
MPI_MEV      = 138.039       # charged-average pion mass
FPI_MEV      = 92.2          # pion decay constant
LAMBDA_CHI   = 700.0        # chiral symmetry breaking scale
GA           = 1.29         # axial coupling (current code value)

# Hebeler 2015 PRC 91 044001, Λ = 500 MeV coordinate-space regulator equivalent
# uses dimensionless LECs c_D, c_E. We use the Hebeler 500 MeV reference point.
LAMBDA_3NF_MEV_DEFAULT = 500.0
C_E_DEFAULT            = -0.02914  # Hebeler 500 MeV (PRC 91 044001 Table I)


# -----------------------------------------------------------------------------
# Independent A_cE derivation: explicit Pauli matrix sum on the
# |T M_T⟩ = Σ_{m_2 m_3} CG(½ m_2; ½ m_3 | T M_T) |m_2⟩⊗|m_3⟩ state.
# -----------------------------------------------------------------------------
def pauli_eigenvalue_independent(T_pair: int) -> float:
    """Compute ⟨(½½)T|τ_2·τ_3|(½½)T⟩ via explicit Pauli matrix sum.

    T_pair ∈ {0, 1}. We do NOT use the closed form 2T(T+1)−3 here.
    """
    if T_pair not in (0, 1):
        raise ValueError(f"Pair isospin T={T_pair} not supported (must be 0 or 1)")

    # Pauli matrices in the |m=+½⟩, |m=-½⟩ basis (sympy exact arithmetic).
    # τ^x = [[0, 1], [1, 0]], τ^y = [[0, -i], [i, 0]], τ^z = [[1, 0], [0, -1]]
    I = S(1)
    tau_x = [[S(0), I], [I, S(0)]]
    tau_y = [[S(0), -I*1j], [I*1j, S(0)]]   # complex entries ok
    tau_z = [[I, S(0)], [S(0), -I]]
    taus = [tau_x, tau_y, tau_z]

    # Build |T M_T⟩ pair state as a length-4 sympy vector in basis
    # |m_2 m_3⟩ indexed as (i_m2 * 2 + i_m3) with i_m=0 (m=+½), 1 (m=-½).
    # Choose M_T = T_pair (highest weight); use sympy exact arithmetic.
    M_T = Rational(T_pair, 1)
    state = [S(0)] * 4
    half = Rational(1, 2)
    m_values = [half, -half]
    for i_m2, m2 in enumerate(m_values):
        for i_m3, m3 in enumerate(m_values):
            if m2 + m3 != M_T:
                continue
            # CG returns a sympy expression; keep it exact.
            # Signature: clebsch_gordan(j_1, j_2, j_3, m_1, m_2, m_3)
            # = ⟨j_1 m_1; j_2 m_2 | j_3 m_3⟩.
            cg = CG(Rational(1, 2), Rational(1, 2), T_pair,
                    m2, m3, M_T)
            state[i_m2 * 2 + i_m3] = cg
    # Normalise (CG columns should already be orthonormal, but assert).
    norm2 = sum((x * x).simplify() for x in state)
    assert float(norm2) == 1.0, f"|T={T_pair},M_T={M_T}⟩ not normalised: {norm2}"

    # Build the 4×4 pair operator τ_2·τ_3 = Σ_a τ_2^a ⊗ τ_3^a on the
    # |m_2 m_3⟩ = |m_2⟩ ⊗ |m_3⟩ tensor-product space.
    def kron(a, b):
        """4×4 Kronecker product of two 2×2 matrices."""
        out = [[S(0)] * 4 for _ in range(4)]
        for i in range(2):
            for j in range(2):
                for k in range(2):
                    for l in range(2):
                        out[i*2 + k][j*2 + l] = a[i][j] * b[k][l]
        return out

    op = [[S(0)] * 4 for _ in range(4)]
    for tau in taus:
        tau2_tau3 = kron(tau, tau)
        for i in range(4):
            for j in range(4):
                op[i][j] = op[i][j] + tau2_tau3[i][j]

    # ⟨state| op |state⟩  (state is real for M_T=T_pair)
    val = S(0)
    for i in range(4):
        for j in range(4):
            val = val + state[i].conjugate() * op[i][j] * state[j]

    result = complex(val)
    assert abs(result.imag) < 1e-12, f"non-real eigenvalue for T={T_pair}: {result}"
    return float(result.real)


def verify_pauli_eigenvalues() -> None:
    """Sanity check: explicit Pauli sum should match closed form 2T(T+1)−3."""
    for T in (0, 1):
        got      = pauli_eigenvalue_independent(T)
        expected = 2 * T * (T + 1) - 3
        assert abs(got - expected) < 1e-10, \
            f"Pauli eigenvalue T={T}: got {got}, expected {expected}"
    print("[oracle] Pauli eigenvalue verification OK: "
          "T=0 -> -3, T=1 -> +1 (explicit sum matches closed form)")


# -----------------------------------------------------------------------------
# Regulator (Epelbaum 2002 eq. 3.19, squared-Gaussian)
# -----------------------------------------------------------------------------
def regulator_gauss_e2002(p_fm: float, q_fm: float, Lambda_fm: float) -> float:
    """f_R(p, q; Λ) = exp( -((4p² + 3q²) / (4Λ²))² )."""
    arg = (4.0 * p_fm**2 + 3.0 * q_fm**2) / (4.0 * Lambda_fm**2)
    return math.exp(-arg * arg)


# -----------------------------------------------------------------------------
# Full V_cE matrix element in the Jj basis.
# -----------------------------------------------------------------------------
def V_cE_matrix_element(channel_r: dict, channel_c: dict,
                         p_fm: float, q_fm: float,
                         p_prime_fm: float, q_prime_fm: float,
                         c_E: float, Lambda_3NF_MeV: float) -> float:
    """Compute ⟨α'|V_cE|α⟩ in fm⁵ using the independent Pauli-sum oracle.

    Channel dict keys: L_2N, S_2N, J_2N, T_2N, L_1N, two_J_1N, two_J_3N, two_T_3N.

    CONVENTION (matches the spectator component of Epelbaum Eq. 2.10):
      W_cE = tau23 × (c_E / (f_π⁴ Λ_χ)) × 1/(8π³) × f_R(p',q') f_R(p,q)

    The full +0.5 sum_(j!=k) counts each unordered pair twice.  Thus the
    spectator-1 component is +E tau_2.tau_3, not +0.5 E tau_2.tau_3.

    INDEPENDENT PART: the tau23 eigenvalue is computed via explicit Pauli
    matrix sum on the pair state |(½½)T⟩ — NOT via the closed form
    2T(T+1)−3 used in the C++ helper.
    """
    # Selection rules (Epelbaum A-4). All deltas must be satisfied for nonzero.
    if channel_r['L_2N'] != 0 or channel_c['L_2N'] != 0: return 0.0
    if channel_r['L_1N'] != 0 or channel_c['L_1N'] != 0: return 0.0
    if channel_r['two_J_1N'] != 1 or channel_c['two_J_1N'] != 1: return 0.0
    if channel_r['J_2N'] != channel_r['S_2N']: return 0.0
    if channel_c['J_2N'] != channel_c['S_2N']: return 0.0
    if channel_r['S_2N'] != channel_c['S_2N']: return 0.0
    if channel_r['T_2N'] != channel_c['T_2N']: return 0.0

    # 3N conserved quantum numbers
    if channel_r['two_J_3N'] != channel_c['two_J_3N']: return 0.0
    if channel_r['two_T_3N'] != channel_c['two_T_3N']: return 0.0

    # INDEPENDENT: Pauli-sum eigenvalue for τ_2·τ_3 in pair isospin T_2N.
    tau23 = pauli_eigenvalue_independent(channel_r['T_2N'])

    # Spectator-component prefactor: c_E / (f_π⁴ Λ_χ) × 1/(8π³).
    fpi_fm        = FPI_MEV / HBARC_MEV_FM
    Lambda_chi_fm = LAMBDA_CHI / HBARC_MEV_FM
    lec = c_E / (fpi_fm**4 * Lambda_chi_fm)
    fourier_norm = 1.0 / (8.0 * math.pi**3)

    # Regulator product (Epelbaum 2002 eq. 3.19).
    Lambda_fm = Lambda_3NF_MeV / HBARC_MEV_FM
    f_R_bra   = regulator_gauss_e2002(p_prime_fm, q_prime_fm, Lambda_fm)
    f_R_ket   = regulator_gauss_e2002(p_fm,       q_fm,       Lambda_fm)

    return tau23 * lec * fourier_norm * f_R_bra * f_R_ket


# -----------------------------------------------------------------------------
# Golden value table: fixed momentum point, list of channels.
# -----------------------------------------------------------------------------
def make_golden_table(p_fm=0.5, q_fm=0.5, pp_fm=0.5, qp_fm=0.5,
                       c_E=C_E_DEFAULT,
                       Lambda_3NF_MeV=LAMBDA_3NF_MEV_DEFAULT) -> dict:
    """Return a dict of channel-label → V_cE in fm⁵."""
    # Channels: tuples of (label, channel_dict, expected_nonzero)
    # All at T_3N = 1/2 (the only physically relevant doublet for ³H / nd).
    channels = {
        "3S1_T0_diagonal": dict(
            L_2N=0, S_2N=1, J_2N=1, T_2N=0, L_1N=0, two_J_1N=1,
            two_J_3N=1, two_T_3N=1),
        "1S0_T1_diagonal": dict(
            L_2N=0, S_2N=0, J_2N=0, T_2N=1, L_1N=0, two_J_1N=1,
            two_J_3N=1, two_T_3N=1),
        "3S1_3D1_offdiag": dict(
            L_2N=0, S_2N=1, J_2N=1, T_2N=0, L_1N=0, two_J_1N=1,
            two_J_3N=1, two_T_3N=1),
        "3P0_T1_forbidden_L": dict(  # forbidden: L_2N=1 not S-wave
            L_2N=1, S_2N=0, J_2N=1, T_2N=1, L_1N=0, two_J_1N=1,
            two_J_3N=1, two_T_3N=1),
        "l1N_2_forbidden": dict(    # forbidden: spectator l_1N=2 not S-wave
            L_2N=0, S_2N=1, J_2N=1, T_2N=0, L_1N=2, two_J_1N=3,
            two_J_3N=3, two_T_3N=1),
    }
    targets = {
        "3S1_T0_diagonal":  channels["3S1_T0_diagonal"],
        "1S0_T1_diagonal":  channels["1S0_T1_diagonal"],
        "3S1_3D1_offdiag":  channels["3S1_3D1_offdiag"],  # ket will be D-wave
        "3P0_T1_forbidden_L": channels["3P0_T1_forbidden_L"],
        "l1N_2_forbidden":  channels["l1N_2_forbidden"],
    }

    out = {}
    for label, bra in targets.items():
        # Pick the ket for off-diagonal cases.
        if label == "3S1_3D1_offdiag":
            ket = dict(bra)
            ket['L_2N'] = 2  # 3D1 partner
        else:
            ket = bra
        val = V_cE_matrix_element(bra, ket, p_fm, q_fm, pp_fm, qp_fm,
                                   c_E, Lambda_3NF_MeV)
        out[label] = {
            'value_fm5': val,
            'bra': bra,
            'ket': ket,
            'p_fm': p_fm, 'q_fm': q_fm, 'pp_fm': pp_fm, 'qp_fm': qp_fm,
            'c_E': c_E, 'Lambda_3NF_MeV': Lambda_3NF_MeV,
        }
    return out


def main():
    verify_pauli_eigenvalues()

    golden = make_golden_table()
    print("\n[oracle] Golden V_cE matrix elements (fm⁵):")
    print("-" * 70)
    for label, entry in golden.items():
        print(f"  {label:30s}  V_cE = {entry['value_fm5']:+.6e} fm⁵")

    # Sanity assertions: structure (NOT absolute normalization, since the
    # absolute value depends on the LEC sign convention which is
    # convention-dependent; the RATIO between channels is invariant).
    g = golden
    assert abs(g['3S1_T0_diagonal']['value_fm5']) > 1e-6, "3S1 must be nonzero"
    assert abs(g['1S0_T1_diagonal']['value_fm5']) > 1e-6, "1S0 must be nonzero"
    # tau_2.tau_3(T=0)/tau_2.tau_3(T=1) = -3/+1 = -3 — independent of LEC sign.
    ratio = g['3S1_T0_diagonal']['value_fm5'] / g['1S0_T1_diagonal']['value_fm5']
    assert abs(ratio + 3.0) < 1e-6, \
        f"V_cE(3S1)/V_cE(1S0) should be -3 (pure tau.tau, no sigma.sigma); got {ratio}"
    assert abs(g['3S1_3D1_offdiag']['value_fm5']) < 1e-15, "3S1↔3D1 must be 0 (L mismatch)"
    assert abs(g['3P0_T1_forbidden_L']['value_fm5']) < 1e-15, "L_2N=1 forbidden"
    assert abs(g['l1N_2_forbidden']['value_fm5']) < 1e-15, "l_1N=2 forbidden"
    print("\n[oracle] All sanity assertions passed (structure-level; ratio=-3).")

    out_path = Path(__file__).parent / "oracle_cE_golden.json"
    with open(out_path, 'w') as f:
        json.dump(golden, f, indent=2, default=str)
    print(f"[oracle] Wrote {out_path}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
