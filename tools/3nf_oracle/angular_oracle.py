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
2. c_D 1PE-contact, spectator S wave: the pair vertex is a contact (L_2N=0),
   and for l=l'=0 the only angular integral is over the spectator direction
   (x = cos q.q'). The enabled production sector does this exactly. Higher-l
   and rank-2 pieces remain fail-closed.
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
MPI_MEV = 138.039             # pion mass [MeV]  (include/constants.h)
FPI_MEV = 92.2                # f_pi [MeV]       (include/constants.h)
GA = 1.289                    # g_A              (include/constants.h)
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
def oracle_cD_rank0_swave(cD, Lambda, p, q, pp, qp):
    """Independent c_D 1S0->3S1 spectator-S-wave matrix element.

    Explicit Pauli-basis states give
      <S23=0|sigma1.sigma3|S23=1> = sqrt(3),
      <T23=1|tau1.tau3|T23=0>     = sqrt(3),
    hence the rank-0 recoupling is +1 after the factor 1/3. The spectator
    integral is evaluated in closed form, independently of production GL.
    """
    if q <= 0.0 or qp <= 0.0:
        raise ValueError("closed-form cD oracle requires q,q' > 0")
    numerator = (q + qp) ** 2 + MPI ** 2
    denominator = (q - qp) ** 2 + MPI ** 2
    integ = math.log(numerator / denominator) / (2.0 * q * qp)
    coeff = -GA * cD / (8.0 * FPI**4 * LAMBDA_CHI) * 2.0
    fR = regulator_gauss(p, q, Lambda) * regulator_gauss(pp, qp, Lambda)
    recoupling = (math.sqrt(3.0) * math.sqrt(3.0)) / 3.0
    return coeff * fR * recoupling * FOURIER_NORM * integ


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
    Lambda_fm = Lambda / HBARC
    cE = -0.02914
    c1_gev, c3_gev = -0.81, -3.2
    c1_fm = c1_gev * HBARC / 1000.0
    c3_fm = c3_gev * HBARC / 1000.0
    a_3S1 = resolve_channel(cE, 0, 0, 0, Lambda, 0,1,1,0, 0,1, 1,1)
    a_1S0 = resolve_channel(cE, 0, 0, 0, Lambda, 0,0,0,1, 0,1, 1,1)
    if a_3S1 is None or a_1S0 is None:
        print("[oracle] 1S0/3S1 channel not found; aborting"); return

    # =========================================================================
    # CALIBRATION 1: c_E (closed-form, no angular integral).
    # =========================================================================
    print("--- calibration 1: c_E contact (closed-form) ---")
    mom = (0.5, 0.4, 0.6, 0.7)
    v_prod = get_production_w1(cE, 0, 0, 0, Lambda, (a_3S1, a_3S1), mom)
    v_oracle = oracle_cE(cE, Lambda_fm, *mom, S_pair=1, T_pair=0)
    r = v_oracle/v_prod if v_prod else float('nan')
    print(f"  c_E 3S1: prod={v_prod:.6e} oracle={v_oracle:.6e}  oracle/prod={r:.4f}")
    print(f"  -> {'PASS (within constant precision)' if abs(r-1)<0.02 else 'FAIL'}")
    print()

    # =========================================================================
    # CALIBRATION 2: c_D 1PE-contact. Pair vertex is contact (L_2N=0), so the
    # ONLY angular dependence is the spectator x = cos(q.q'). The production
    # code integrates over x exactly. The full angular integral (= spectator
    # x-integral with NO pair-angle dependence) MUST match production.
    # This is the decisive normalization calibration.
    # =========================================================================
    print("--- calibration 2: c_D rank-0 1S0->3S1 (Pauli + closed form) ---")
    cD = -0.2
    for (p, q, pp, qp) in [(0.5,0.4,0.6,0.7), (1.0,1.0,1.0,1.0), (0.3,0.8,0.5,0.3)]:
        v_prod_cD = get_production_w1(0.0, cD, 0.0, 0.0, Lambda,
                                     (a_1S0, a_3S1), (p,q,pp,qp))
        v_oracle_cD = oracle_cD_rank0_swave(cD, Lambda_fm, p,q,pp,qp)
        ratio = v_prod_cD / v_oracle_cD if v_prod_cD is not None else float('nan')
        verdict = "PASS" if abs(ratio - 1.0) < 2e-6 else "FAIL"
        print(f"  p={p} q={q} p'={pp} q'={qp}: prod={v_prod_cD:.6e} "
              f"oracle={v_oracle_cD:.6e} prod/oracle={ratio:.8f} {verdict}")
    print()

    # =========================================================================
    # KEY CHECK: c_1/c_3 2PE rank-0, multi-point scan.
    # For each momentum point: compute the full pair-angular integral and
    # compare to the monopole-approximated integral. The ratio
    #   R = V_full / V_mono = [int dOmega/(4pi) int dx kernel_full]
    #                          / [int dx kernel_mono]
    # is 1 if the monopole approximation is exact. R > 1 means the production
    # code UNDERESTIMATES (it misses aligned-momentum configs where the pion
    # propagator blows up).
    # =========================================================================
    print("--- key check: c_1/c_3 2PE rank-0 monopole-approximation scan ---")
    print(f"  c1={c1_gev} GeV^-1, c3={c3_gev} GeV^-1, Lambda={Lambda} MeV")
    print(f"  {'p':>5} {'q':>5} {'p\'':>5} {'q\'':>5} | "
          f"{'V_prod':>12} {'V_mono_oracle':>14} {'V_full_oracle':>14} | "
          f"{'prod/mono':>9} {'prod/full':>9} {'B7 err%':>8}")
    print("  " + "-" * 95)
    points = [
        (0.5, 0.4, 0.6, 0.7),
        (0.5, 0.5, 0.5, 0.5),   # diagonal
        (1.0, 0.5, 0.8, 0.6),
        (0.3, 0.3, 0.7, 0.7),   # large transfer
        (1.0, 1.0, 1.5, 1.5),   # high momentum
        (0.2, 0.2, 0.2, 0.2),   # low momentum (near threshold)
    ]
    si_3S1 = sigma2_dot_sigma3(1) * tau2_dot_tau3(0) / 3.0
    prefactor = (GA / (2.0 * FPI)) ** 2 * FOURIER_NORM
    for (p, q, pp, qp) in points:
        I_full, fR_f = oracle_2pe_rank0_full(c1_fm, c3_fm, Lambda_fm, p, q, pp, qp)
        I_mono, fR_m = oracle_2pe_rank0_production_formula(c1_fm, c3_fm, Lambda_fm, p, q, pp, qp)
        v_mono = prefactor * si_3S1 * fR_m * I_mono
        v_full = prefactor * si_3S1 * fR_f * I_full / (4.0 * math.pi)
        v_prod = get_production_w1(0.0, 0.0, c1_gev, c3_gev, Lambda, (a_3S1, a_3S1), (p,q,pp,qp))
        if v_prod and abs(v_mono) > 1e-30 and abs(v_full) > 1e-30:
            pm = v_prod / v_mono
            pf = v_prod / v_full
            b7err = 100.0 * (1.0 - pf)   # positive = production underestimates
            print(f"  {p:5.2f} {q:5.2f} {pp:5.2f} {qp:5.2f} | "
                  f"{v_prod:12.4e} {v_mono:14.4e} {v_full:14.4e} | "
                  f"{pm:9.4f} {pf:9.4f} {b7err:+7.1f}%")
        else:
            print(f"  {p:5.2f} {q:5.2f} {pp:5.2f} {qp:5.2f} | (zero or None)")
    print()
    print("  Interpretation:")
    print("    prod/mono ~ 1.00  : independent re-derivation of the production formula OK")
    print("    prod/full < 1     : production UNDERESTIMATES the full angular integral")
    print("                         (B7 monopole approximation: |Delta p|^2 -> p^2+p'^2")
    print("                          inside the nonlinear pion propagator)")
    print("    prod/full > 1     : production OVERESTIMATES")
    print()
    print("=" * 72)
    print("Phase 3 oracle complete.")
    print("=" * 72)


if __name__ == "__main__":
    run()
