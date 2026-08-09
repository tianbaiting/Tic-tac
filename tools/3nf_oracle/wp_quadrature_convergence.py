#!/usr/bin/env python3
"""
tools/3nf_oracle/wp_quadrature_convergence.py — Phase 4 WP-cell quadrature
convergence report.

Tests whether the legacy 1-point midpoint rule (Np_per_WP_W1 = Nq_per_WP_W1 = 1)
is converged for the W^(1) bin matrix element, by
comparing it to higher-order Gauss-Legendre cell quadrature (N=2,4,8 points per
bin per dimension).

For each (p_bin, q_bin, p_bin', q_bin') the bin-averaged W^(1)_WP is
    W1_WP = [sum over N^4 quad points] w_p w_q w_p' w_q' * p*q*p'*q' * W1(p,q,p',q')
            * 1/sqrt(dp*dq*dp'*dq') / hbarc^5
mirroring src/interactions/w1_pw_cache.cpp:228-260. The W1 values come from the
production C++ driver print_w1_element (so we test the REAL W1_element, not a
re-transcription).

The c_E contact kernel is momentum-independent (only the regulator varies), so
its bin average should converge instantly (midpoint ~ exact). The c_D and c1/c3
kernels have pion propagators that vary within the bin, so the midpoint may not
be converged — this quantifies the discretization error.
"""
import math
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DRIVER = REPO / "build" / "tools" / "3nf_oracle" / "print_w1_element"
HBARC = 197.327

try:
    import numpy as np
    def gl_nodes(n, a, b):
        x, w = np.polynomial.legendre.leggauss(n)
        return 0.5*(b-a)*x + 0.5*(a+b), 0.5*(b-a)*w
except ImportError:
    sys.exit("requires numpy")


def batch_w1(cE, cD, c1, c3, Lambda, points):
    """Call print_w1_element in batch mode; return list of W1 values."""
    a_3S1 = 4  # 3S1 channel index in the test pw_statespace (L0 S1 J1 T0 l0 2j1 2J3=1 2T3=1)
    inp = "".join(f"{a_3S1} {a_3S1} {p:.10f} {q:.10f} {pp:.10f} {qp:.10f}\n"
                  for (p, q, pp, qp) in points)
    r = subprocess.run([str(DRIVER), str(cE), str(cD), str(c1), str(c3), str(Lambda)],
                       input=inp, capture_output=True, text=True, timeout=120)
    vals = []
    for line in r.stdout.splitlines():
        if line.startswith("W1 "):
            vals.append(float(line.split()[-1]))
    return vals


def bin_average(cE, cD, c1, c3, Lambda, bin_bounds, N):
    """Compute the bin-averaged W1_WP for N-point Gauss quadrature per dimension.

    bin_bounds = (p_lo, p_hi, q_lo, q_hi, pp_lo, pp_hi, qp_lo, qp_hi) in fm^-1.
    """
    p_lo, p_hi, q_lo, q_hi, pp_lo, pp_hi, qp_lo, qp_hi = bin_bounds
    px, pw = gl_nodes(N, p_lo, p_hi)
    qx, qw = gl_nodes(N, q_lo, q_hi)
    ppx, ppw = gl_nodes(N, pp_lo, pp_hi)
    qpx, qpw = gl_nodes(N, qp_lo, qp_hi)
    points = []
    weights = []
    for pi, pwi in zip(px, pw):
        for qi, qwi in zip(qx, qw):
            for ppi, ppiw in zip(ppx, ppw):
                for qpi, qpiw in zip(qpx, qpw):
                    points.append((pi, qi, ppi, qpi))
                    weights.append((pwi, qwi, ppiw, qpiw))
    vals = batch_w1(cE, cD, c1, c3, Lambda, points)
    inv_hbarc5 = (1.0/HBARC)**5
    dp = p_hi - p_lo; dq = q_hi - q_lo; dpp = pp_hi - pp_lo; dqp = qp_hi - qp_lo
    bin_norm = 1.0 / math.sqrt(dp * dq * dpp * dqp)
    accum = 0.0
    for v, (p, q, pp, qp), (pwi, qwi, ppiw, qpiw) in zip(vals, points, weights):
        accum += (p * q * pp * qp) * (pwi * qwi * ppiw * qpiw) * v
    return accum * bin_norm * inv_hbarc5


def convergence_table(cE, cD, c1, c3, Lambda, bin_bounds, label):
    """Print convergence of the bin average for N=1,2,4,8."""
    print(f"\n  {label}")
    print(f"    bin: p=[{bin_bounds[0]:.3f},{bin_bounds[1]:.3f}] "
          f"q=[{bin_bounds[2]:.3f},{bin_bounds[3]:.3f}] "
          f"p'=[{bin_bounds[4]:.3f},{bin_bounds[5]:.3f}] "
          f"q'=[{bin_bounds[6]:.3f},{bin_bounds[7]:.3f}] fm^-1")
    vals = {}
    for N in [1, 2, 4, 8]:
        vals[N] = bin_average(cE, cD, c1, c3, Lambda, bin_bounds, N)
    ref = vals[8]
    print(f"    {'N':>3} {'W1_WP [MeV]':>16} {'rel. diff vs N=8':>18}")
    for N in [1, 2, 4, 8]:
        rdiff = abs(vals[N] - ref) / abs(ref) if abs(ref) > 0 else float('nan')
        tag = " (legacy midpoint diagnostic)" if N == 1 else ""
        print(f"    {N:3d} {vals[N]:16.8e} {rdiff:18.2e}{tag}")
    return ref


def run():
    if not DRIVER.exists():
        sys.exit(f"driver not found: {DRIVER} (run cmake --build build --target print_w1_element)")
    print("=" * 72)
    print("Phase 4 — WP-cell quadrature convergence report")
    print("=" * 72)
    Lambda = 500.0
    cE, cD, c1, c3 = -0.02914, 0.0, -0.81, -3.2

    # Representative WP bins (fm^-1). Use small bins (typical of a Chebyshev
    # grid at moderate Np) so the midpoint-vs-converged difference is visible.
    print("\nc_E contact (momentum-independent kernel — should converge instantly):")
    convergence_table(cE, 0, 0, 0, Lambda,
                      (0.40, 0.60, 0.30, 0.50, 0.50, 0.70, 0.40, 0.60),
                      "c_E only")

    print("\nc_1/c_3 2PE (pion propagator varies within the bin):")
    convergence_table(0.0, 0.0, c1, c3, Lambda,
                      (0.40, 0.60, 0.30, 0.50, 0.50, 0.70, 0.40, 0.60),
                      "c1/c3 2PE rank-0")

    print("\nc_1/c_3 2PE, wider bins (more variation within bin):")
    convergence_table(0.0, 0.0, c1, c3, Lambda,
                      (0.20, 0.80, 0.20, 0.80, 0.20, 0.80, 0.20, 0.80),
                      "c1/c3 2PE, wide bins")

    print("\nc_E + c_D + c_1/c_3 combined (production-like):")
    convergence_table(cE, -1.0, c1, c3, Lambda,
                      (0.40, 0.60, 0.30, 0.50, 0.50, 0.70, 0.40, 0.60),
                      "c_E + c_D + c_1/c_3")

    print("\n" + "=" * 72)
    print("Interpretation:")
    print("  rel. diff vs N=8 < 1e-3 : midpoint (N=1) is converged for this cell")
    print("  rel. diff vs N=8 > 1e-2  : midpoint is NOT converged; the legacy")
    print("                              Np_per_WP_W1=Nq_per_WP_W1=1 setting introduces")
    print("                              a discretization error that grows with bin width.")
    print("  Default N=2 is a safer baseline, not a universal convergence certificate;")
    print("  wide bins in this report require N=4 for 1e-3-level convergence.")
    print("=" * 72)


if __name__ == "__main__":
    run()
