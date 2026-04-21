#!/usr/bin/env python3
"""
Hand-calculation oracle for the c_D 1PE-contact piece of the chiral N2LO 3NF
partial-wave matrix element.

For the 3S1 diagonal channel at p=q=p'=q'=0.5 fm^-1:

The c_D operator (E2002 eq. 2.10):
    V^(1)_OPE = -(gA / 8 fpi^2) * D * sum_j [ (tau_1.tau_j) *
                  (sigma_j.q_j)(sigma_1.q_j) / (q_j^2 + mpi^2) ]
    D = c_D / (fpi^2 * Lambda_chi)

After partial-wave projection (l_1N=0, L_2N=0, rank-2 vanishes by CG):

    W1_1pe_contact = coeff * f_R^2 * (1/3 * sig13) * tau13 * integral_rank0

where:
    coeff = -(gA * c_D) / (8 * fpi^4 * Lambda_chi) * 2 (pair j=2,3 sum)
    integral_rank0 = (1/(8pi^3)) * integral_{-1}^{+1} dx / (Q^2(x) + mpi^2)
    Q^2(x) = q^2 + q'^2 - 2*q*q'*x = 0.5*(1-x) at q=q'=0.5

Uses sympy for Wigner 6j/9j symbols.
"""

import math
import numpy as np
from numpy.polynomial.legendre import leggauss

# ------------------------- Wigner symbol implementations ---------------------
# Using closed-form formulas for the specific values needed.
# These are exact for small quantum numbers.

def w6j(a, b, c, d, e, f):
    """Wigner 6j symbol {a b c; d e f} via sympy."""
    from sympy.physics.wigner import wigner_6j
    from sympy import Rational, sqrt
    result = float(wigner_6j(a, b, c, d, e, f))
    return result

def w9j(a, b, c, d, e, f, g, h, i):
    """Wigner 9j symbol {a b c; d e f; g h i} via sympy."""
    from sympy.physics.wigner import wigner_9j
    result = float(wigner_9j(a, b, c, d, e, f, g, h, i))
    return result

def main():
    # ----------------------------- physical constants --------------------------
    hbarc = 197.327                         # MeV fm
    fpi_MeV = 92.2                          # constants.h value
    Lambda_chi_MeV = 700.0
    Lambda_3NF_MeV = 500.0
    c_D = -0.20                             # Witala LEC, dimensionless
    gA_val = 1.29                           # axial coupling

    mpi_MeV = 138.0
    fpi_fm = fpi_MeV / hbarc
    Lambda_chi_fm = Lambda_chi_MeV / hbarc
    Lambda_fm = Lambda_3NF_MeV / hbarc
    mpi_fm = mpi_MeV / hbarc
    fpi4_fm = fpi_fm**4

    print("=" * 72)
    print("c_D 1PE-contact oracle (3S1 diagonal, p=q=p'=q'=0.5 fm^-1)")
    print("=" * 72)
    print(f"  fpi_fm        = {fpi_fm:.6f} fm^-1")
    print(f"  Lambda_chi_fm = {Lambda_chi_fm:.6f} fm^-1")
    print(f"  Lambda_fm     = {Lambda_fm:.6f} fm^-1")
    print(f"  mpi_fm        = {mpi_fm:.6f} fm^-1")
    print(f"  fpi^4         = {fpi4_fm:.6e} fm^-4")
    print()

    # ----------------------------- Fourier normalization -----------------------
    fourier_norm = 1.0 / (8.0 * math.pi**3)

    # ----------------------------- Channel QNs (3S1) ----------------------------
    L = S = 1; J = 1; T = 0; l = 0; two_j = 1; two_J3N = 1; two_T3N = 1
    L_p = S_p = 1; J_p = 1; T_p = 0; l_p = 0; two_j_p = 1

    # ----------------------------- sigma_1.sigma_3 matrix element ---------------
    # Via spin_isospin_algebra.cpp::reduced_me_sigma1_dot_sigma3 formula:
    #
    # sig13 = prefactor * phase * hat_J2N * hat_S * hat_j1N
    #          * w9j * w6j_pair_orbital * w6j_pair_spin * w6j_spec
    # prefactor = -sqrt(3) * sqrt(2*J_3N+1) * 6
    # 2*phase_exp = 2L + 2S' + 2*J_2N + 2S + 2l + 5 + two_j_1N

    w6j_pair_spin    = w6j(1/2, 1, 1/2, 1, 1/2, 1)   # {1/2 S' 1/2; S 1/2 1}
    w6j_pair_orbital = w6j(1, 1, 0, 1, 1, 1)           # {S' J_2N' L; J_2N S 1}
    w6j_spec         = w6j(1/2, 1/2, 0, 1/2, 1/2, 1)  # {1/2 j' l; j 1/2 1}

    # 9j rows: (J_2N, J_2N', 1), (j_1N, j_1N', 1), (J_3N, J_3N, 0)
    w9j_spin = w9j(1, 1, 1, 1/2, 1/2, 1, 1/2, 1/2, 0)

    two_exp = 2*L + 2*S_p + 2*J_p + 2*S + 2*l + 5 + two_j
    phase = 1.0 if (two_exp // 2) % 2 == 0 else -1.0

    hat_J2N = math.sqrt((2*J_p + 1.0) * (2*J + 1.0))
    hat_S   = math.sqrt((2*S_p + 1.0) * (2*S + 1.0))
    hat_j1N = math.sqrt((two_j_p + 1.0) * (two_j + 1.0))
    prefactor_spin = -math.sqrt(3.0) * math.sqrt(two_J3N + 1.0) * 6.0

    sig13 = (prefactor_spin * phase * hat_J2N * hat_S * hat_j1N
             * w9j_spin * w6j_pair_orbital * w6j_pair_spin * w6j_spec)

    print("[sigma_1.sigma_3 matrix element (3S1 diagonal)]")
    print(f"  w6j_pair_spin    = {w6j_pair_spin:.8f}")
    print(f"  w6j_pair_orbital = {w6j_pair_orbital:.8f}")
    print(f"  w6j_spec         = {w6j_spec:.8f}")
    print(f"  w9j              = {w9j_spin:.8f}")
    print(f"  phase_exp = {two_exp//2} -> phase = {phase:+.1f}")
    print(f"  hat_J2N={hat_J2N:.3f}, hat_S={hat_S:.3f}, hat_j1N={hat_j1N:.3f}")
    print(f"  prefactor = {prefactor_spin:.6f}")
    print(f"  sig13 (3S1 diag) = {sig13:.8f}")
    print()

    # ----------------------------- tau_1.tau_3 matrix element ------------------
    # Via spin_isospin_algebra.cpp::tau1_dot_tau3 formula:
    # w6j_tau: {1/2 T' 1/2; T 1/2 1} = {1/2 0 1/2; 0 1/2 1}
    w6j_tau_pair = w6j(1/2, 0, 1/2, 0, 1/2, 1)  # {1/2 T' 1/2; T 1/2 1}
    # 9j: (T, T', 1), (1/2, 1/2, 1), (T_3N, T_3N, 0) = {0 0 1; 1/2 1/2 1; 1/2 1/2 0}
    w9j_tau = w9j(0, 0, 1, 1/2, 1/2, 1, 1/2, 1/2, 0)

    phase_tau = (-1.0)**T    # (-1)^T_2N
    hat_T = math.sqrt((2*T_p + 1.0) * (2*T + 1.0))  # = 1 for T=T'=0
    prefactor_tau = -math.sqrt(3.0) * math.sqrt(two_T3N + 1.0) * 6.0

    tau13 = prefactor_tau * phase_tau * hat_T * w9j_tau * w6j_tau_pair

    print("[tau_1.tau_3 matrix element (3S1 diagonal)]")
    print(f"  w6j_tau_pair = {w6j_tau_pair:.8f}")
    print(f"  w9j_tau      = {w9j_tau:.8f}")
    print(f"  phase_tau    = {phase_tau:+.1f}")
    print(f"  tau13 (3S1 diag) = {tau13:.8f}")
    print()

    # ----------------------------- Rank-0 recoupling ---------------------------
    # Factor of 1/3 from (sigma_1.q)(sigma_3.q) rank-0 decomposition = 1/3 sigma_1.sigma_3
    recoup_scalar = (1.0/3.0) * sig13 * tau13
    print(f"[Rank-0 recoupling coefficient] = (1/3)*sig13*tau13 = {recoup_scalar:.8f}")
    print()

    # ----------------------------- Rank-2 check (l=0 kills it) ----------------
    print("[Rank-2 check: CG(l=0,m=0; k=2,m=0 | l'=0,m=0)]")
    print("  Triangle rule: |0-2|<=0<=0+2 => 2<=0: FALSE => CG=0")
    print("  => rank-2 integral = 0 for l_1N=l_1N'=0 (S-wave spectator)")
    print()

    # ----------------------------- Gauss-Legendre integral ---------------------
    # Q^2(x) = q^2 + q'^2 - 2 q q' x
    q_val = q_p = 0.5   # fm^-1
    p_val = p_p = 0.5

    n_gl = 24
    x_nodes, w_weights = leggauss(n_gl)

    integral_rank0 = 0.0
    for xi, wi in zip(x_nodes, w_weights):
        Q2 = q_val**2 + q_p**2 - 2.0*q_val*q_p*xi
        integral_rank0 += wi * fourier_norm / (Q2 + mpi_fm**2)

    print(f"[Gauss-Legendre integral (n={n_gl})]")
    print(f"  Q^2(x) = {q_val}^2 + {q_p}^2 - 2*{q_val}*{q_p}*x = 0.5*(1-x)")
    print(f"  mpi^2  = {mpi_fm**2:.6f} fm^-2")
    print(f"  integral_rank0 = (1/(8pi^3)) * int dx/(Q^2+mpi^2)")
    print(f"                 = {integral_rank0:.8e}  fm^2")
    print()

    # ----------------------------- Regulator ------------------------------------
    inv_L2 = 1.0 / Lambda_fm**2
    a = (p_val**2 + 0.75*q_val**2) * inv_L2
    f_R = math.exp(-a**2)
    print(f"[Regulator f_R(p=0.5, q=0.5)]")
    print(f"  a = (p^2 + 0.75*q^2)/Lambda^2 = {a:.6f}")
    print(f"  f_R = exp(-a^2) = {f_R:.8f}")
    print(f"  f_R^2 = {f_R**2:.8f}")
    print()

    # ----------------------------- Overall coefficient -------------------------
    # -(gA * c_D) / (8 * fpi^4 * Lambda_chi) * 2 (summing j=2,3)
    coeff = -(gA_val * c_D) / (8.0 * fpi4_fm * Lambda_chi_fm) * 2.0
    print(f"[Overall coefficient]")
    print(f"  coeff = -(gA*c_D)/(8*fpi^4*Lambda_chi) * 2")
    print(f"        = -({gA_val}*{c_D})/(8*{fpi4_fm:.4e}*{Lambda_chi_fm:.4f}) * 2")
    print(f"        = {coeff:.6e}  [dimensionless * fm^5 from the kernel]")
    print()

    # ----------------------------- Final assembly --------------------------------
    W1_cD = coeff * f_R**2 * recoup_scalar * integral_rank0
    print("=" * 72)
    print("W1_1pe_contact(3S1, p=q=p'=q'=0.5 fm^-1)")
    print(f"  = coeff * f_R^2 * recoup_scalar * integral_rank0")
    print(f"  = {coeff:.4e} * {f_R**2:.6f} * {recoup_scalar:.6f} * {integral_rank0:.4e}")
    print(f"  = {W1_cD:.6e}  fm^5")
    print("=" * 72)

    return W1_cD


if __name__ == "__main__":
    val = main()
    print(f"\nFINAL: W1_1pe_contact(3S1, p=q=p'=q'=0.5) = {val:.6e} fm^5")
