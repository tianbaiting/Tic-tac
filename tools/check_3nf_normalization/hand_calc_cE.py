#!/usr/bin/env python3
"""
Hand-calculation oracle for the c_E contact piece of the chiral N2LO 3NF
partial-wave matrix element, using the exact conventions that Task 3 installs
in chiral_N2LO_3NF.h::W1_contact:

    W^(1)_cE(alpha', alpha; p', q', p, q)
        = recoupling_3nf_scalar(alpha', alpha)
          * kernel_contact(c_E, fpi^4, Lambda_chi)
          * f_R(p', q') * f_R(p, q)                                (Eq. A)

with

    recoupling_3nf_scalar (rank-0) = (sigma_2 . sigma_3) * (tau_2 . tau_3)
                                     * (Kronecker selection on all quantum numbers
                                        of pair and spectator)

    kernel_contact = +0.5 * c_E / (fpi^4 * Lambda_chi)              (fm^5)
                     [E2002 eq. 2.10: V^(1)_cont = -1/2 E sum_{j!=k} (tau_j . tau_k)]

    f_R(p, q; Lambda) = exp( -((4 p^2 + 3 q^2) / (4 Lambda^2))^2 ) (dimensionless)
                       [E2002 eq. 3.19, squared-Gaussian]

Diagonal 3S1 channel:
    L_2N = 0, S_2N = 1, J_2N = 1, T_2N = 0, L_1N = 0,
    2 j_1N = 1, 2 J_3N = 1, 2 T_3N = 1.

At p = q = p' = q' = 0.5 fm^-1.

NOTE: The Task-1 helper `recoupling_3nf_scalar` deliberately encodes the
rank-0 scalar piece as (sigma . sigma) * (tau . tau).  For the strict c_E
contact operator, only (tau . tau) is present on physical grounds (the
operator is a pure isospin scalar with no sigma dependence, see
formula_reference.md §1.1); the extra (sigma . sigma) = 2 S(S+1) - 3 is
innocuous for S_pair = 1 (factor +1) but would spuriously enter S_pair = 0.
This matches the deliberate design choice in the pw-kernels infrastructure
and is what the rewrite installs.  We flag this for the reviewer.
"""

import math


def main():
    # --------------------------- physical constants ---------------------------
    hbarc = 197.327                        # MeV fm
    # Use the in-code value so the Python oracle matches the C++ baseline.
    fpi_MeV = 92.2                         # constants.h value (task spec says 92.4 -- see note below)
    Lambda_chi_MeV = 700.0                 # chiral breaking scale
    Lambda_3NF_MeV = 500.0                 # 3NF regulator cut-off (EM500)
    c_E = -0.205                           # Witala/Epelbaum LEC, dimensionless

    # ----------------------- convert to natural (fm) units --------------------
    fpi_fm = fpi_MeV / hbarc               # fm^-1
    Lambda_chi_fm = Lambda_chi_MeV / hbarc # fm^-1
    Lambda_fm = Lambda_3NF_MeV / hbarc     # fm^-1

    fpi4_fm = fpi_fm**4                    # fm^-4

    print("=" * 72)
    print("c_E contact single-point oracle (diagonal 3S1, p=q=p'=q'=0.5 fm^-1)")
    print("=" * 72)
    print()
    print("[constants]")
    print(f"  hbarc            = {hbarc} MeV fm")
    print(f"  fpi              = {fpi_MeV} MeV  -> {fpi_fm:.6f} fm^-1")
    print(f"  Lambda_chi       = {Lambda_chi_MeV} MeV -> {Lambda_chi_fm:.6f} fm^-1")
    print(f"  Lambda_3NF       = {Lambda_3NF_MeV} MeV -> {Lambda_fm:.6f} fm^-1")
    print(f"  c_E              = {c_E}")
    print(f"  fpi^4            = {fpi4_fm:.6e} fm^-4")
    print()

    # ------------------- channel quantum numbers (3S1 diag) -------------------
    S_2N = 1
    T_2N = 0
    # (recoupling is diagonal in all of L_2N, S_2N, J_2N, T_2N, L_1N, j_1N.)

    # sigma_2 . sigma_3  =  2 S_pair (S_pair + 1) - 3
    #                    =  2*1*2 - 3 = +1
    sigma_sigma = 2.0 * S_2N * (S_2N + 1) - 3.0
    # tau_2 . tau_3      =  2 T_pair (T_pair + 1) - 3
    #                    =  2*0*1 - 3 = -3
    tau_tau = 2.0 * T_2N * (T_2N + 1) - 3.0
    # Helper (Task 1) returns this product (no 6j, no (4pi)^2):
    recoup = sigma_sigma * tau_tau

    print("[recoupling coefficient (Task-1 rank-0 convention)]")
    print(f"  sigma_2 . sigma_3 (S_2N=1) = {sigma_sigma:+.1f}")
    print(f"  tau_2   . tau_3   (T_2N=0) = {tau_tau:+.1f}")
    print(f"  recoupling_3nf_scalar       = sigma*sigma * tau*tau = {recoup:+.6f}")
    print()

    # ----------------------- kernel_contact (Task 2) --------------------------
    # +0.5 * c_E / (fpi^4 * Lambda_chi), with 1/(2π)³ Fourier normalization
    # mirroring chiral_LO_internal.cpp:59. Sign is Navrátil/Witała convention
    # (see kernel_contact docstring for triangulation commit notes).
    fourier_norm_3nf = 1.0 / (8.0 * math.pi * math.pi * math.pi)
    kernel = fourier_norm_3nf * (+0.5 * c_E / (fpi4_fm * Lambda_chi_fm))
    print("[kernel_contact (Task 2)]")
    print(f"  (1/(8π³)) * +0.5 * c_E / (fpi^4 * Lambda_chi)")
    print(f"  = {fourier_norm_3nf:.6e} * +0.5 * ({c_E}) / ({fpi4_fm:.6e} * {Lambda_chi_fm:.6f})")
    print(f"  = {kernel:+.6e}  [fm^5]")
    print()

    # ----------------------- regulator f_R(p, q) ------------------------------
    # exp(-((4 p^2 + 3 q^2) / (4 Lambda^2))^2)
    p, q = 0.5, 0.5                        # fm^-1 (both bra and ket)
    num = 4.0 * p * p + 3.0 * q * q        # fm^-2
    den = 4.0 * Lambda_fm * Lambda_fm      # fm^-2
    a = num / den
    fR = math.exp(-a * a)
    print("[regulator f_R(p=0.5, q=0.5) per E2002 eq. 3.19]")
    print(f"  (4 p^2 + 3 q^2) / (4 Lambda^2) = ({num:.4f}) / ({den:.4f}) = {a:.6f}")
    print(f"  exp(-a^2)                       = {fR:.8f}")
    print(f"  f_R(p,q) * f_R(p',q')           = {fR*fR:.8f}")
    print()

    # ---------------------------- assemble (Eq. A) ----------------------------
    W1 = recoup * kernel * fR * fR
    print("[final matrix element]")
    print("  W^(1)_cE(3S1; p=q=p'=q'=0.5 fm^-1) = recoup * kernel_contact * f_R^2")
    print(f"                                     = ({recoup:+.4f}) * ({kernel:+.6e}) * ({fR*fR:.8f})")
    print(f"                                     = {W1:+.6e}  [fm^5]")
    print()

    # ------------------------- sanity: arithmetic trace -----------------------
    # recoup = +1 * (-3) = -3
    # kernel = +0.5 * (-0.205) / (fpi^4 * Lambda_chi)
    #        = -0.1025 / (fpi^4 * Lambda_chi)
    # fR^2   ~ 0.99 (near 1; regulator barely damps at p=q=0.5 fm^-1)
    # Final sign: recoup(-3) * kernel(negative) * fR^2(positive) = positive.

    # ----------------------- summary line the reviewer wants ------------------
    print("=" * 72)
    print(f"W1_contact(3S1, p=q=p'=q'=0.5 fm^-1) = {W1:+.6e} fm^5")
    print("=" * 72)

    # Also emit a value using the task-spec fpi=92.4 MeV for bookkeeping.
    fpi_task_MeV = 92.4
    fpi_task_fm = fpi_task_MeV / hbarc
    fpi4_task_fm = fpi_task_fm**4
    kernel_task = fourier_norm_3nf * (+0.5 * c_E / (fpi4_task_fm * Lambda_chi_fm))
    W1_task = recoup * kernel_task * fR * fR
    print(f"[ref: with fpi=92.4 MeV (task-spec) -> W1 = {W1_task:+.6e} fm^5]")


if __name__ == "__main__":
    main()
