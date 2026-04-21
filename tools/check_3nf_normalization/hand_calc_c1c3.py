#!/usr/bin/env python3
"""
Hand-calculation oracle for the c_1 and c_3 2PE pieces of the chiral N2LO 3NF
partial-wave matrix element.

Computes the RANK-0 contribution for the diagonal 3S1 channel at
p=q=p'=q'=0.5 fm^-1 using Witala LECs. This oracle exercises:
  - lec_bracket = (-4 c_1 m_pi^2 + 2 c_3 q2q3) / fpi^2
  - Rank-0 scalar reduction: 1/3 * (sigma_2.sigma_3) * (tau_2.tau_3)
  - 24-pt Gauss-Legendre quadrature over x = cos(q.q')
  - Squared-Gaussian regulator per E2002 eq. (3.19)
  - Overall (gA/2fpi)^2 prefactor

Expected c_1 rank-0: approximately -0.8 fm^5 (see formula_reference.md §3.5)
Expected c_3 rank-0: approximately +1.2 fm^5 (see formula_reference.md §4.5)

Operator form ([G2010] eq. 18):
    V^(1)_2pi = (gA/2fpi)^2 * (tau_2.tau_3) * (sigma_2.q2)(sigma_3.q3)
                * (-4c1*mpi^2/fpi^2 + 2*c3*(q2.q3)/fpi^2) / [(q2^2+mpi^2)(q3^2+mpi^2)]

Rank-0 decomposition:
    (sigma_2.q2)(sigma_3.q3) -> 1/3 * (sigma_2.sigma_3) * (q2.q3)    [rank-0 piece]

Momentum-transfer kinematics (E2002 azimuthal-averaged convention):
    dp2 = p^2 + p'^2        (|Delta_p|^2, azimuthal averaged)
    dq2 = q^2 + q'^2 - 2*q*q'*x   (|Delta_q|^2(x))
    q2sq = dp2 + dq2/4      (<q2^2> = <q3^2>)
    q2q3 = dp2 - dq2/4      (<q2.q3>)
"""

import math
import numpy as np
from numpy.polynomial.legendre import leggauss

def main():
    # ----------------------------- Physical constants --------------------------
    hbarc        = 197.327          # MeV fm
    fpi_MeV      = 92.2             # pion decay constant, MeV (constants.h value)
    Lambda_chi_MeV = 700.0          # chiral breaking scale
    Lambda_3NF_MeV = 500.0          # 3NF regulator cutoff
    mpi_MeV      = 138.0            # pion mass MeV
    gA_val       = 1.29             # axial coupling (constants.h value)

    # Witala LECs (c_1, c_3 in GeV^-1 -> fm)
    c1_GeV_inv = -0.81              # GeV^-1
    c3_GeV_inv = -3.20              # GeV^-1
    c1_fm = c1_GeV_inv * hbarc / 1000.0   # fm
    c3_fm = c3_GeV_inv * hbarc / 1000.0   # fm

    fpi_fm       = fpi_MeV   / hbarc
    Lambda_fm    = Lambda_3NF_MeV / hbarc
    mpi_fm       = mpi_MeV   / hbarc
    fpi2_fm      = fpi_fm**2
    fpi4_fm      = fpi_fm**4

    print("=" * 72)
    print("c_1 / c_3 2PE oracle (rank-0, 3S1 diagonal, p=q=p'=q'=0.5 fm^-1)")
    print("=" * 72)
    print(f"  fpi_fm    = {fpi_fm:.6f} fm^-1")
    print(f"  Lambda_fm = {Lambda_fm:.6f} fm^-1")
    print(f"  mpi_fm    = {mpi_fm:.6f} fm^-1")
    print(f"  c1_fm     = {c1_fm:.6f} fm")
    print(f"  c3_fm     = {c3_fm:.6f} fm")
    print()

    # ----------------------------- Channel: 3S1 diagonal ----------------------
    # alpha = {L_2N=0, S_2N=1, J_2N=1, T_2N=0, l_1N=0, 2j_1N=1, 2J_3N=1, 2T_3N=1}
    S_2N = 1; T_2N = 0
    sigma23 = 2.0 * S_2N * (S_2N + 1) - 3.0   # +1 for triplet
    tau23   = 2.0 * T_2N * (T_2N + 1) - 3.0   # -3 for isoscalar pair

    # Rank-0 coefficient: 1/3 * sigma_2.sigma_3 * tau_2.tau_3
    recoup_rank0 = (1.0/3.0) * sigma23 * tau23
    print(f"[Rank-0 recoupling (3S1)]")
    print(f"  sigma_2.sigma_3 = {sigma23:.1f}  (S=1 triplet)")
    print(f"  tau_2.tau_3     = {tau23:.1f}  (T=0 isoscalar)")
    print(f"  1/3 * sig * tau = {recoup_rank0:.6f}")
    print()

    # ----------------------------- Momentum kinematics -----------------------
    p_val = q_val = p_p = q_p = 0.5   # all 0.5 fm^-1
    n_gl = 24
    x_nodes, w_weights = leggauss(n_gl)
    fourier_norm = 1.0 / (8.0 * math.pi**3)

    dp2  = p_val**2 + p_p**2   # = 0.5 fm^-2

    # ----------------------------- c_1-only rank-0 ---------------------------
    integ_c1 = 0.0
    integ_c3 = 0.0
    integ_c1c3 = 0.0
    for xi, wi in zip(x_nodes, w_weights):
        dq2  = q_val**2 + q_p**2 - 2.0*q_val*q_p*xi   # |Delta_q|^2(x)
        q2sq = dp2 + 0.25*dq2                            # <q2^2> = <q3^2>
        q2q3 = dp2 - 0.25*dq2                            # <q2.q3>
        mp2  = mpi_fm**2
        prop = 1.0 / ((q2sq + mp2)**2)                   # 1/(q2^2+mp^2)^2

        lec_bracket_c1   = (-4.0 * c1_fm * mp2 + 0) / fpi2_fm
        lec_bracket_c3   = (0 + 2.0 * c3_fm * q2q3) / fpi2_fm
        lec_bracket_both = (-4.0 * c1_fm * mp2 + 2.0 * c3_fm * q2q3) / fpi2_fm

        integ_c1   += wi * fourier_norm * lec_bracket_c1   * q2q3 * prop
        integ_c3   += wi * fourier_norm * lec_bracket_c3   * q2q3 * prop
        integ_c1c3 += wi * fourier_norm * lec_bracket_both * q2q3 * prop

    print(f"[Gauss-Legendre integral (n={n_gl})]")
    print(f"  integ_rank0(c1 only)   = {integ_c1:.8e}  fm^3")
    print(f"  integ_rank0(c3 only)   = {integ_c3:.8e}  fm^3")
    print(f"  integ_rank0(c1+c3)     = {integ_c1c3:.8e}  fm^3")
    print()

    # ----------------------------- Regulator ---------------------------------
    inv_L2 = 1.0 / Lambda_fm**2
    a = (p_val**2 + 0.75*q_val**2) * inv_L2
    f_R = math.exp(-a**2)
    print(f"[Regulator f_R(p=0.5, q=0.5)]")
    print(f"  a = (p^2 + 0.75*q^2)/Lambda^2 = {a:.6f}")
    print(f"  f_R = exp(-a^2) = {f_R:.8f},  f_R^2 = {f_R**2:.8f}")
    print()

    # ----------------------------- Overall coefficient -----------------------
    # (gA/2fpi)^2 = gA^2 / (4 * fpi^2)
    coeff = gA_val**2 / (4.0 * fpi2_fm)
    print(f"[Overall coefficient (gA/2fpi)^2]")
    print(f"  coeff = gA^2 / (4*fpi^2) = {gA_val:.3f}^2 / (4*{fpi2_fm:.4e}) = {coeff:.6e} fm^2")
    print()

    # ----------------------------- Final assembly ----------------------------
    W1_c1   = coeff * f_R**2 * recoup_rank0 * integ_c1
    W1_c3   = coeff * f_R**2 * recoup_rank0 * integ_c3
    W1_c1c3 = coeff * f_R**2 * recoup_rank0 * integ_c1c3

    print("=" * 72)
    print("Rank-0 oracle values (3S1 diagonal, p=q=p'=q'=0.5 fm^-1):")
    print(f"  W1_2pe rank-0 (c_1 only)   = {W1_c1:.6e} fm^5  (expect ~-0.8 fm^5)")
    print(f"  W1_2pe rank-0 (c_3 only)   = {W1_c3:.6e} fm^5  (expect ~+1.2 fm^5)")
    print(f"  W1_2pe rank-0 (c_1+c_3)    = {W1_c1c3:.6e} fm^5")
    print("=" * 72)
    print()
    print("Note: the C++ W1_2pe also includes rank-2 contributions (off-diagonal alpha).")
    print("The 3S1->3S1 diagonal element should approximately match the rank-0 oracle")
    print("above (rank-2 vanishes for L=0->L=0 by CG selection rule CG(0,0;2,0|0,0)=0).")
    print("The 3S1->3D1 off-diagonal element is the rank-2 tensor contribution.")

    return W1_c1, W1_c3, W1_c1c3


if __name__ == "__main__":
    v_c1, v_c3, v_c1c3 = main()
    print()
    print(f"FINAL: W1_2pe rank-0 (c_1 only, 3S1) = {v_c1:.6e} fm^5")
    print(f"FINAL: W1_2pe rank-0 (c_3 only, 3S1) = {v_c3:.6e} fm^5")
    print(f"FINAL: W1_2pe rank-0 (c_1+c_3, 3S1)  = {v_c1c3:.6e} fm^5")
