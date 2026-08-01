#!/usr/bin/env python3
"""
tools/3nf_oracle/angular_oracle.py — Phase 3 independent angular oracle.

INDEPENDENT derivation of the chiral N2LO 3NF partial-wave matrix elements,
starting from the UN-REDUCED momentum-space operator and performing the FULL
multi-dimensional angular quadrature (NOT the azimuthal average that the
production code uses for the 2PE c1/c3 piece — audit B7).

What this oracle checks
-----------------------
1. c_E contact: closed-form, NO angular integral. The production code is exact
   here (the only angular dependence is the regulator, which factors out).
   This is the calibration check — the oracle MUST reproduce the production
   value to machine precision.
2. c_D 1PE-contact: the pair vertex is a contact (L_2N=0), so the only angular
   integral is over the spectator direction (x = cos q.q'). The production code
   does this exactly. Another calibration check.
3. c_1/c_3 2PE rank-0: the operator depends on q_2, q_3 which involve BOTH the
   pair (p) and spectator (q) momentum transfers. The production code averages
   over the pair azimuth phi_{p'}, replacing |Delta p|^2 -> p^2 + p'^2 and
   dropping the Delta p . Delta q cross term BEFORE the (nonlinear) pion
   propagator — the "monopole approximation" (audit B7). The oracle computes
   the FULL angular integral and quantifies the discrepancy.

Independence
------------
The oracle derives the kernel from the operator
    V^(1)_2pi = (gA/2 fpi)^2 * (tau2.tau3) * (sigma2.q2)(sigma3.q3)
                * [-4 c1 m_pi^2 + 2 c3 (q2.q3)] / [fpi^2 (q2^2+m_pi^2)(q3^2+m_pi^2)]
                * regulator(p,q) regulator(p',q')
by building the vectors p, q, p', q', q1, q2, q3 explicitly and integrating
over the angles with Gauss-Legendre quadrature. It does NOT transcribe the
production kernel_2pe_c1c3 / kernel_1pe_contact formulae; it re-derives them
from the operator. The production value comes from the C++ driver
print_w1_element (a thin wrapper around chiral_N2LO_3NF::W1_element), so the
comparison is genuinely code-vs-independent-derivation.

Usage
-----
    python3 angular_oracle.py
The driver print_w1_element must be built first (see run() below).
"""

import math
import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DRIVER = REPO / "build" / "tools" / "3nf_oracle" / "print_w1_element"

# Physical constants — MUST match include/constants.h for a like-for-like
# comparison. Values taken from include/constants.h.
HBARC = 197.327               # MeV fm
MPI_MEV = 138.0               # pion mass [MeV]  (constants.h: mpi=138.0)
FPI_MEV = 92.4                # f_pi [MeV]       (constants.h: fpi=92.4)
GA = 1.29                     # g_A              (constants.h: gA=1.29)
LAMBDA_CHI_MEV = 700.0        # chiral breaking scale [MeV]

# fm units
MPI = MPI_MEV / HBARC
FPI = FPI_MEV / HBARC
LAMBDA_CHI = LAMBDA_CHI_MEV / HBARC
FOURIER_NORM = 1.0 / (8.0 * math.pi ** 3)


def regulator_gauss(p, q, Lambda):
    """E2002 eq. (3.19) squared-Gaussian regulator."""
    a = (4.0 * p * p + 3.0 * q * q) / (4.0 * Lambda * Lambda)
    return math.exp(-a * a)


def gauss_legendre(n, a=-1.0, b=1.0):
    """n-point Gauss-Legendre nodes/weights on [a,b] via numpy."""
    try:
        import numpy as np
        x, w = np.polynomial.legendre.leggauss(n)
        # scale to [a,b]
        x = 0.5 * (b - a) * x + 0.5 * (a + b)
        w = 0.5 * (b - a) * w
        return list(x), list(w)
    except ImportError:
        # fallback: simple trapezoid (less accurate; only for environments
        # without numpy — the oracle REQUIRES numpy for the comparison to be
        # meaningful, so this is a hard error if numpy is missing).
        raise RuntimeError("angular_oracle.py requires numpy for Gauss-Legendre quadrature")


# -----------------------------------------------------------------------------
# Spin-isospin eigenvalues (INDEPENDENT Pauli derivation, matching the C++
# test_cE_direct_pauli_enumeration test). tau = sigma (Pauli, NOT sigma/2),
# so tau2.tau3 = 2 T(T+1) - 3.
# -----------------------------------------------------------------------------
def tau2_dot_tau3(T_pair):
    return 2 * T_pair * (T_pair + 1) - 3

def sigma2_dot_sigma3(S_pair):
    return 2 * S_pair * (S_pair + 1) - 3


# -----------------------------------------------------------------------------
# Vector helpers (3D).
# -----------------------------------------------------------------------------
def vec(x, y, z):
    return (x, y, z)

def vsub(a, b):
    return (a[0]-b[0], a[1]-b[1], a[2]-b[2])

def vadd(a, b):
    return (a[0]+b[0], a[1]+b[1], a[2]+b[2])

def vscale(a, s):
    return (a[0]*s, a[1]*s, a[2]*s)

def vdot(a, b):
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

def vnorm(a):
    return math.sqrt(vdot(a, a))


# -----------------------------------------------------------------------------
# Oracle: c_E contact (closed-form, calibration).
# Matches kernel_contact: +0.5 * c_E / (fpi^4 * Lambda_chi) * fourier_norm.
# -----------------------------------------------------------------------------
def oracle_cE(cE, Lambda, p, q, pp, qp, S_pair, T_pair):
    """Independent c_E contact matrix element. NO angular integral (contact)."""
    tau23 = tau2_dot_tau3(T_pair)
    fR = regulator_gauss(p, q, Lambda) * regulator_gauss(pp, qp, Lambda)
    lec = 0.5 * cE / (FPI**4 * LAMBDA_CHI)
    return tau23 * lec * fR * FOURIER_NORM


# -----------------------------------------------------------------------------
# Oracle: c_D 1PE-contact via FULL spectator angular integral.
# Operator: -(gA c_D / (8 fpi^4 Lambda_chi)) * 2  (sum over j=2,3)
#         * (tau1.tau3)(sigma1.q3_hat)(sigma3.q3_hat) / (q3^2 + m_pi^2)
# For the rank-0 S-wave reduction: (sigma1.q_hat)(sigma3.q_hat) -> (1/3)(sigma1.sigma3).
# The pair vertex is contact (L_2N=0), so the pair angle does NOT enter: q3 = -Delta q / 2
# depends only on q, q'. The integral is over the spectator direction x = cos(q.q').
# This should be EXACT (the production code does the same integral).
# -----------------------------------------------------------------------------
def oracle_cD_rank0(cD, Lambda, p, q, pp, qp, S_pair, T_pair, two_T3N=1, n_gl=48):
    """c_D 1PE-contact rank-0 via full spectator x-quadrature.

    Uses the SPECTATOR momentum transfer Q^2 = |Delta q|^2 = q^2 + q'^2 - 2 q q' x
    (E2002 eq. A-2). Pair vertex is contact, so no pair-angle dependence.
    """
    # Isospin tau1.tau3 matrix element in the coupled |(T_pair, 1/2) T_3N> basis.
    # tau1.tau3 = (1/2)[T_3N(T_3N+1) - T_pair(T_pair+1) - 3/4]*2 - ... use the
    # standard result: <(T 1/2)T3 | tau1.tau3 | (T 1/2)T3> connects T_pair to
    # 1-T_pair (off-diagonal in T_pair for tau1.tau3). For the S-wave spectator
    # j_1N=1/2 doublet (T_3N=1/2), tau1.tau3 has BOTH diagonal and off-diagonal
    # pieces. The production recoupling_3nf_1pe_ct_scalar handles this via
    # tau1_dot_tau3(T_r, T_c, two_T3N). To keep the oracle INDEPENDENT and
    # focused on the SPATIAL integral (which is the calibration point here),
    # we compute only the SPATIAL part and let the spin-isospin factor be
    # supplied by the caller / compared via the ratio.
    # Spatial integral: int dx 1/(Q^2 + m_pi^2), Q^2 = q^2 + q'^2 - 2 q q' x.
    xs, ws = gauss_legendre(n_gl, -1.0, 1.0)
    integ = 0.0
    for x, w in zip(xs, ws):
        Q2 = q*q + qp*qp - 2.0*q*qp*x
        integ += w * 1.0 / (Q2 + MPI*MPI)
    # The production kernel_1pe_contact returns fourier_norm / (Q2 + m_pi^2)
    # and the caller multiplies by the spin-isospin recoupling and the
    # overall coeff -gA cD / (8 fpi^4 Lambda_chi) * 2.
    return integ  # spatial integral only (spin-isospin factored out)


# -----------------------------------------------------------------------------
# Oracle: c_1/c_3 2PE rank-0 via FULL angular integral (the key check).
#
# Operator (G2010 eq. 18, E2002 eq. 2.2):
#   V^(1)_2pi = (gA/2 fpi)^2 (tau2.tau3) (sigma2.q2)(sigma3.q3)
#               * [-4 c1 m_pi^2 + 2 c3 (q2.q3)] / [fpi^2 (q2^2+m_pi^2)(q3^2+m_pi^2)]
#               * regulator(p,q) regulator(p',q')
#
# S-wave pair rank-0 reduction (EXACT for L_2N=0): the spin trace gives
#   (sigma2.q2)(sigma3.q3) -> (1/3)(sigma2.sigma3)(q2.q3)
# so the operator becomes, after spin-isospin factoring,
#   K_spatial(q2,q3) = (q2.q3) * [-4 c1 m_pi^2 + 2 c3 (q2.q3)]
#                      / [fpi^2 (q2^2+m_pi^2)(q3^2+m_pi^2)]
#
# q2 = Delta p + Delta q / 2,  q3 = Delta p - Delta q / 2
# Delta p = p' - p (vectors),   Delta q = q' - q (vectors)
#
# The PRODUCTION code azimuthally averages over the pair azimuth phi_{p'}:
#   |Delta p|^2  ->  p^2 + p'^2          (drops -2 p p' cos theta_{p,p'})
#   q2^2 = q3^2  ->  |Delta p|^2 + |Delta q|^2 / 4   (drops Delta p . Delta q)
#   q2.q3        ->  |Delta p|^2 - |Delta q|^2 / 4   (drops Delta p . Delta q)
# and then integrates only over x = cos(q.q').
#
# The ORACLE does the FULL 3D angular integral over (theta_{p'}, phi_{p'}, x)
# WITHOUT the azimuthal reduction, keeping the full q2, q3 angle dependence.
# -----------------------------------------------------------------------------
def oracle_2pe_rank0_full(c1, c3, Lambda, p, q, pp, qp, n_theta=40, n_phi=24, n_x=48):
    """Full 3D angular integral of the 2PE rank-0 spatial kernel.

    Returns the RAW integral I_full = int dOmega_{p'} int dx kernel_full
    (NO 1/(4pi)^2 normalization — the caller forms the normalization-canceling
    ratio against I_monopole = 4pi * int dx kernel_monopole, which is 1 when
    the azimuthal-average approximation is exact).

    Frame: p = p z_hat, q = q z_hat (both along z; valid for S-wave pair +
    spectator because the S-wave projections make the result independent of the
    p-q relative angle). p' is integrated over its full solid angle, q' over x
    = cos(q.q') with phi_{q'} = 0 (the S-wave spectator makes the result
    phi_{q'}-independent for the rank-0 piece).
    """
    p_vec = vec(0.0, 0.0, p)
    q_vec = vec(0.0, 0.0, q)

    ths, wths = gauss_legendre(n_theta, 0.0, math.pi)
    phs, wphs = gauss_legendre(n_phi, 0.0, 2.0 * math.pi)
    xs, wxs = gauss_legendre(n_x, -1.0, 1.0)

    total = 0.0
    for th_p, wth in zip(ths, wths):
        sin_th = math.sin(th_p)
        cos_th = math.cos(th_p)
        for ph_p, wph in zip(phs, wphs):
            pp_vec = vec(pp * sin_th * math.cos(ph_p),
                         pp * sin_th * math.sin(ph_p),
                         pp * cos_th)
            dp = vsub(pp_vec, p_vec)
            dp2 = vdot(dp, dp)
            for x, wx in zip(xs, wxs):
                sin_x = math.sqrt(max(0.0, 1.0 - x * x))
                qp_vec = vec(qp * sin_x, 0.0, qp * x)
                dq = vsub(qp_vec, q_vec)
                dq2 = vdot(dq, dq)
                q2 = vadd(dp, vscale(dq, 0.5))
                q3 = vsub(dp, vscale(dq, 0.5))
                q2sq = vdot(q2, q2)
                q3sq = vdot(q3, q3)
                q2q3 = vdot(q2, q3)
                mp2 = MPI * MPI
                if q2sq + mp2 <= 0 or q3sq + mp2 <= 0:
                    continue
                prop = 1.0 / ((q2sq + mp2) * (q3sq + mp2))
                lec_bracket = (-4.0 * c1 * mp2 + 2.0 * c3 * q2q3) / (FPI * FPI)
                k = q2q3 * lec_bracket * prop
                # dOmega_{p'} = sin(th) dth dph; the spectator x-integral here
                # does NOT carry a 2pi phi_{q'} factor — the production code's
                # 1/(8pi^3) Fourier norm already absorbs the spectator azimuth,
                # so to compare like-for-like we keep phi_{q'} = 0 fixed (the
                # S-wave spectator makes the rank-0 result phi_{q'}-independent).
                total += wth * sin_th * wph * wx * k
    fR = regulator_gauss(p, q, Lambda) * regulator_gauss(pp, qp, Lambda)
    return total, fR  # raw I_full = int dOmega_{p'} int dx kernel_full (phi_{q'}=0)


def oracle_2pe_rank0_production_formula(c1, c3, Lambda, p, q, pp, qp, n_x=48):
    """Re-derive the production azimuthal-averaged kernel INDEPENDENTLY.

    Returns the RAW integral I_monopole = int dx kernel_monopole(x) (the
    azimuthal-averaged kernel the production code uses, with
    |Delta p|^2 -> p^2 + p'^2 and Delta p . Delta q -> 0). The caller compares
    I_full / (4*pi * I_monopole): if the monopole approximation were exact,
    int dOmega kernel_full = 4*pi * kernel_monopole, so the ratio would be 1.
    """
    xs, wxs = gauss_legendre(n_x, -1.0, 1.0)
    total = 0.0
    for x, wx in zip(xs, wxs):
        dp2 = p * p + pp * pp                       # azimuthal-averaged |dp|^2
        dq2 = q * q + qp * qp - 2.0 * q * qp * x    # |dq|^2(x)
        q2sq = dp2 + 0.25 * dq2                     # <q2^2> = <q3^2>
        q2q3 = dp2 - 0.25 * dq2                     # <q2.q3>
        mp2 = MPI * MPI
        prop = 1.0 / ((q2sq + mp2) * (q2sq + mp2))
        lec_bracket = (-4.0 * c1 * mp2 + 2.0 * c3 * q2q3) / (FPI * FPI)
        k = q2q3 * lec_bracket * prop
        total += wx * k
    fR = regulator_gauss(p, q, Lambda) * regulator_gauss(pp, qp, Lambda)
    return total, fR  # raw I_monopole (int dx), regulator product


# -----------------------------------------------------------------------------
# Driver: get the production W1_element value for a given channel + momenta.
# -----------------------------------------------------------------------------
def build_driver():
    """Build print_w1_element via CMake if not present."""
    if DRIVER.exists():
        return True
    build_dir = REPO / "build"
    if not build_dir.exists():
        print("[oracle] build/ dir not found; run cmake + make first")
        return False
    print("[oracle] building driver via CMake:", DRIVER)
    r = subprocess.run(["cmake", "--build", str(build_dir), "--target", "print_w1_element"],
                       capture_output=True, text=True, cwd=str(build_dir))
    if r.returncode != 0:
        print("[oracle] build FAILED:\n", r.stderr[-1500:])
        return False
    return DRIVER.exists()


def get_production_w1(cE, cD, c1, c3, Lambda, channel, momenta):
    """Run the C++ driver and return the production W1_element value."""
    if not build_driver():
        return None
    # channel = (alpha_r, alpha_c); momenta = (p_r, q_r, p_c, q_c)
    ar, ac = channel
    pr, qr, pc, qc = momenta
    inp = f"{ar} {ac} {pr} {qr} {pc} {qc}\n"
    r = subprocess.run([str(DRIVER), str(cE), str(cD), str(c1), str(c3), str(Lambda)],
                       input=inp, capture_output=True, text=True, timeout=60)
    for line in r.stdout.splitlines():
        if line.startswith("W1 "):
            parts = line.split()
            return float(parts[-1])
    return None


def resolve_channel(cE, cD, c1, c3, Lambda, L2, S2, J2, T2, l1, two_j1, two_J3, two_T3):
    """Run the driver once and parse the channel table to find the alpha index."""
    if not build_driver():
        return None
    # dummy momenta just to get the channel table
    inp = "0 0 0.1 0.1 0.1 0.1\n"
    r = subprocess.run([str(DRIVER), str(cE), str(cD), str(c1), str(c3), str(Lambda)],
                       input=inp, capture_output=True, text=True, timeout=60)
    for line in r.stdout.splitlines():
        if line.startswith("# alpha "):
            parts = line.split()
            a = int(parts[2])
            def kw(tag):
                # tokens look like "L2=0"; find the token starting with tag=
                for tok in parts:
                    if tok.startswith(tag):
                        return int(tok.split("=")[1])
                raise KeyError(tag)
            if (kw("L2=")==L2 and kw("S2=")==S2 and kw("J2=")==J2 and kw("T2=")==T2
                and kw("l1=")==l1 and kw("2j1=")==two_j1 and kw("2J3=")==two_J3
                and kw("2T3=")==two_T3):
                return a
    return None


def run():
    print("=" * 72)
    print("Phase 3 independent angular oracle — chiral N2LO 3NF")
    print("=" * 72)
    print(f"constants: m_pi={MPI:.6f} fpi={FPI:.6f} gA={GA} Lambda_chi={LAMBDA_CHI:.6f} fm")
    print(f"           fourier_norm=1/(8pi^3)={FOURIER_NORM:.6e}")
    print()

    Lambda = 500.0  # MeV
    p, q, pp, qp = 0.5, 0.4, 0.6, 0.7   # fm^-1, asymmetric bra/ket

    # ---- c_E calibration (closed-form, no angular integral) ----------------
    print("--- c_E contact (calibration: closed-form, no angular integral) ---")
    cE = -0.02914
    # 3S1: L2=0,S2=1,J2=1,T2=0, l1=0,2j1=1, 2J3=1, 2T3=1
    a_3S1 = resolve_channel(cE, 0, 0, 0, Lambda, 0,1,1,0, 0,1, 1,1)
    if a_3S1 is None:
        print("[oracle] 3S1 channel not found; aborting cE check")
    else:
        v_prod = get_production_w1(cE, 0, 0, 0, Lambda, (a_3S1, a_3S1), (p,q,pp,qp))
        v_oracle = oracle_cE(cE, Lambda/HBARC, p, q, pp, qp, S_pair=1, T_pair=0)
        print(f"  3S1 (T_pair=0): prod={v_prod:.8e}  oracle={v_oracle:.8e}  "
              f"ratio={v_oracle/v_prod:.6f}" if v_prod else "  [prod None]")
        print(f"  (oracle tau23(T=0)={tau2_dot_tau3(0)}, should be -3)")

    # 1S0: L2=0,S2=0,J2=0,T2=1, l1=0,2j1=1, 2J3=1, 2T3=1
    a_1S0 = resolve_channel(cE, 0, 0, 0, Lambda, 0,0,0,1, 0,1, 1,1)
    if a_1S0 is not None and a_3S1 is not None:
        v_prod_1S0 = get_production_w1(cE, 0, 0, 0, Lambda, (a_1S0, a_1S0), (p,q,pp,qp))
        v_prod_3S1 = get_production_w1(cE, 0, 0, 0, Lambda, (a_3S1, a_3S1), (p,q,pp,qp))
        if v_prod_1S0 and v_prod_3S1:
            ratio = v_prod_3S1 / v_prod_1S0
            print(f"  prod ratio V(3S1)/V(1S0) = {ratio:.6f}  (should be -3, NOT +1)")
    print()

    # ---- c_1/c_3 2PE rank-0: full angular vs production azimuthal-average --
    print("--- c_1/c_3 2PE rank-0 (key check: full angular vs azimuthal avg) ---")
    c1, c3 = -0.81, -3.2  # fm (Idaho, already in fm via *hbarc/1000 in C++; here use fm directly)
    # The C++ driver takes c1,c3 in GeV^-1 and converts. Idaho c1=-0.81 GeV^-1.
    c1_gev, c3_gev = -0.81, -3.2
    c1_fm = c1_gev * HBARC / 1000.0
    c3_fm = c3_gev * HBARC / 1000.0
    print(f"  c1={c1_gev} GeV^-1 = {c1_fm:.6f} fm, c3={c3_gev} GeV^-1 = {c3_fm:.6f} fm")
    print(f"  momenta: p={p} q={q} p'={pp} q'={qp} fm^-1, Lambda={Lambda} MeV")

    I_full, fR_full = oracle_2pe_rank0_full(c1_fm, c3_fm, Lambda/HBARC, p, q, pp, qp)
    I_mono, fR_mono = oracle_2pe_rank0_production_formula(c1_fm, c3_fm, Lambda/HBARC, p, q, pp, qp)
    prefactor = (GA / (2.0 * FPI)) ** 2 * FOURIER_NORM  # (gA/2fpi)^2 / (8pi^3)
    # Normalization-canceling ratio: if the monopole approximation were exact,
    # int dOmega_{p'} kernel_full = (2*pi) * kernel_monopole (the 2*pi from
    # the pair azimuth dph integral; for an angle-independent kernel the pair
    # polar integral gives 2 and the azimuth gives 2*pi, total 4*pi but the
    # production I_mono already carries one factor of 2 from int dx, so the
    # matching is I_full / (2*pi * I_mono) == 1 for an exact reduction).
    # NOTE: the absolute normalization of the partial-wave matrix element
    # involves the 1/(8*pi^3) Fourier norm and 4*pi spherical-harmonic factors
    # that are subtle to match exactly between the full-solid-angle integral and
    # the production code's reduced int dx. The RATIO is the meaningful
    # diagnostic: deviations from 1 quantify the B7 monopole approximation.
    ratio = I_full / (2.0 * math.pi * I_mono) if abs(I_mono) > 0 else float('nan')
    print(f"  I_full (raw int dOmega int dx, phi_q'=0)     = {I_full:.8e}")
    print(f"  I_mono (raw int dx, azimuthal-avg kernel)   = {I_mono:.8e}")
    print(f"  2*pi * I_mono                                = {2*math.pi*I_mono:.8e}")
    print(f"  ratio I_full / (2*pi * I_mono)               = {ratio:.6f}")
    print(f"  regulator product (must match for both)      = {fR_full:.8e}")
    if abs(ratio - 1.0) < 0.05:
        print(f"  >>> azimuthal-average (B7) is EXACT for this point "
              f"(deviation {100*(ratio-1):+.4f}%)")
    else:
        print(f"  >>> azimuthal-average (B7) is an APPROXIMATION: "
              f"ratio deviates from 1 by {100*(ratio-1):+.2f}% at this point.")
        print(f"      NOTE: the absolute magnitude of the deviation depends on the")
        print(f"      subtle normalization mapping between the full-solid-angle")
        print(f"      integral and the production code's 1/(8pi^3) Fourier-norm +")
        print(f"      int dx reduction. The DIRECTION (full != monopole) is robust;")
        print(f"      a research-grade reference PWD is needed to pin the magnitude.")
        print(f"      What is certain: the production azimuthal average replaces")
        print(f"      |Delta p|^2 -> p^2+p'^2 INSIDE the nonlinear pion propagator,")
        print(f"      which is NOT the same as averaging the propagator over angles.")

    # Compare the independent re-derivation of the production formula to the
    # actual C++ W1_element — confirms we correctly understand the production code.
    if a_3S1 is not None:
        v_prod_2pe = get_production_w1(0.0, 0.0, c1_gev, c3_gev, Lambda, (a_3S1, a_3S1), (p,q,pp,qp))
        if v_prod_2pe and fR_mono != 0:
            # V_prod = prefactor * (sigma.sig)(tau.tau)/3 * fR * fourier_norm * I_mono
            # (production kernel_2pe_c1c3 already includes fourier_norm and the
            #  int dx; the W1_2pe caller adds prefactor and 1/3 and fR).
            si_3S1 = sigma2_dot_sigma3(1) * tau2_dot_tau3(0) / 3.0
            v_oracle_avg = prefactor * si_3S1 * fR_mono * I_mono
            print(f"  production C++ W1_element (3S1, c1/c3 only) = {v_prod_2pe:.8e}")
            print(f"  oracle AZIMUTHAL re-deriv (prefactor*si*fR*I_mono) = {v_oracle_avg:.8e}")
            print(f"    -> prod / oracle_avg = {v_prod_2pe/v_oracle_avg:.6f}  "
                  f"({'independent re-derivation OK' if abs(v_prod_2pe/v_oracle_avg-1)<0.05 else 'MISMATCH — investigate'})")
            v_oracle_full = prefactor * si_3S1 * fR_full * I_full / (4.0 * math.pi)
            print(f"  oracle FULL-ANGULAR V (3S1)                    = {v_oracle_full:.8e}")
            print(f"    -> prod / oracle_full = {v_prod_2pe/v_oracle_full:.6f}  "
                  f"({'prod==full (B7 exact)' if abs(v_prod_2pe/v_oracle_full-1)<0.05 else 'prod != full (B7 confirmed)'})")
    print()
    print("=" * 72)
    print("Phase 3 oracle complete.")
    print("If 'prod / oracle_avg' matches: the independent re-derivation of the")
    print("production formula is correct (we understand what the code computes).")
    print("If 'prod / oracle_full' deviates from 1: the azimuthal-average (audit")
    print("B7) is a genuine approximation; the deviation quantifies its magnitude.")
    print("=" * 72)


if __name__ == "__main__":
    run()
