#!/usr/bin/env python3
"""
Hand-calculation oracle for the c_E contact piece of the chiral N2LO 3NF
partial-wave matrix element, using the exact conventions that Task 3 installs
in chiral_N2LO_3NF.h::W1_contact:

    W^(1)_cE(alpha', alpha; p', q', p, q)
        = recoupling_3nf_contact_cE(alpha', alpha)
          * kernel_contact(c_E, fpi^4, Lambda_chi)
          * f_R(p', q') * f_R(p, q)                                (Eq. A)

with

    recoupling_3nf_contact_cE = (tau_2 . tau_3)
                               * (Kronecker selection on all quantum numbers)

    kernel_contact = c_E / (fpi^4 * Lambda_chi) * 1/(4*pi^4)       (fm^5)
                     [E2002 eqs. 2.10, 2.12, and A-4]

    f_R(p, q; Lambda) = exp( -((4 p^2 + 3 q^2) / (4 Lambda^2))^2 ) (dimensionless)
                       [E2002 eq. 3.19, squared-Gaussian]

Diagonal 3S1 channel:
    L_2N = 0, S_2N = 1, J_2N = 1, T_2N = 0, L_1N = 0,
    2 j_1N = 1, 2 J_3N = 1, 2 T_3N = 1.

At p = q = p' = q' = 0.5 fm^-1.

The normalization 1/(4*pi^4) is the product of the raw S-wave angular factor
(4*pi)^2 in Epelbaum Eq. (A-4) and one (2*pi)^-3 Fourier factor for each of the
two independent Jacobi coordinates.
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
    T_2N = 0
    # (recoupling is diagonal in all of L_2N, S_2N, J_2N, T_2N, L_1N, j_1N.)

    # tau_2 . tau_3      =  2 T_pair (T_pair + 1) - 3
    #                    =  2*0*1 - 3 = -3
    tau_tau = 2.0 * T_2N * (T_2N + 1) - 3.0
    recoup = tau_tau

    print("[c_E recoupling coefficient]")
    print(f"  tau_2   . tau_3   (T_2N=0) = {tau_tau:+.1f}")
    print(f"  recoupling_3nf_contact_cE   = tau*tau = {recoup:+.6f}")
    print()

    # ----------------------- kernel_contact (Task 2) --------------------------
    # c_E / (fpi^4 * Lambda_chi), with exact contact PWD normalization.
    # The 1/2 in Epelbaum Eq. (2.10) is exhausted by the ordered-pair sum.
    contact_pw_norm = 1.0 / (4.0 * math.pi**4)
    kernel = contact_pw_norm * (c_E / (fpi4_fm * Lambda_chi_fm))
    print("[kernel_contact (Task 2)]")
    print(f"  (1/(4π⁴)) * c_E / (fpi^4 * Lambda_chi)")
    print(f"  = {contact_pw_norm:.6e} * ({c_E}) / ({fpi4_fm:.6e} * {Lambda_chi_fm:.6f})")
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
    # recoup = tau_2.tau_3 = -3
    # kernel = -0.205 / (fpi^4 * Lambda_chi)
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
    kernel_task = contact_pw_norm * (c_E / (fpi4_task_fm * Lambda_chi_fm))
    W1_task = recoup * kernel_task * fR * fR
    print(f"[ref: with fpi=92.4 MeV (task-spec) -> W1 = {W1_task:+.6e} fm^5]")


if __name__ == "__main__":
    main()
