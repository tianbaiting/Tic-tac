#!/usr/bin/env python3
"""
tools/3nf_oracle/wp_quadrature_convergence.py — Phase 4 WP-cell quadrature
convergence report.

Tests whether the legacy 1-point midpoint rule (Np_per_WP_W1 = Nq_per_WP_W1 = 1)
is converged for the W^(1) bin matrix element, by comparing it to higher-order
Gauss-Legendre cell quadrature.  Both the legacy approximation and the complete
Hebeler-factorized projector are supported.

For each (p_bin, q_bin, p_bin', q_bin') the bin-averaged W^(1)_WP is
    W1_WP = [sum over N^4 quad points] w_p w_q w_p' w_q' * p*q*p'*q' * W1(p,q,p',q')
            * 1/sqrt(dp*dq*dp'*dq') * hbarc,
where this script's momenta and weights are in fm^-1.  This is exactly the
production MeV-variable formula, whose final factor is hbarc^-5, after changing
all four integration variables and the WP normalization.  The W1 values come
from the production C++ driver print_w1_element (so we test the REAL W1_element,
not a re-transcription).

The c_E contact operator is momentum-independent apart from the regulator, so
its bin average generally converges rapidly. The c_D and c1/c3
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


def batch_w1(cE, cD, c1, c3, Lambda, points, *, c4=0.0,
             factorized=False, transfer_order=6, alpha_r=4, alpha_c=4):
    """Call print_w1_element in batch mode; return list of W1 values."""
    inp = "".join(f"{alpha_r} {alpha_c} {p:.10f} {q:.10f} {pp:.10f} {qp:.10f}\n"
                  for (p, q, pp, qp) in points)
    command = [str(DRIVER)]
    if factorized:
        command += ["--factorized", str(cE), str(cD), str(c1), str(c3),
                    str(c4), str(Lambda), str(transfer_order)]
    else:
        command += [str(cE), str(cD), str(c1), str(c3), str(Lambda)]
    r = subprocess.run(command,
                       input=inp, capture_output=True, text=True, timeout=120)
    vals = []
    for line in r.stdout.splitlines():
        if line.startswith("W1 "):
            vals.append(float(line.split()[-1]))
    return vals


def bin_average(cE, cD, c1, c3, Lambda, bin_bounds, N, **driver_options):
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
    vals = batch_w1(cE, cD, c1, c3, Lambda, points, **driver_options)
    if len(vals) != len(points):
        raise RuntimeError(f"driver returned {len(vals)} values for {len(points)} points")
    dp = p_hi - p_lo; dq = q_hi - q_lo; dpp = pp_hi - pp_lo; dqp = qp_hi - qp_lo
    bin_norm = 1.0 / math.sqrt(dp * dq * dpp * dqp)
    accum = 0.0
    for v, (p, q, pp, qp), (pwi, qwi, ppiw, qpiw) in zip(vals, points, weights):
        accum += (p * q * pp * qp) * (pwi * qwi * ppiw * qpiw) * v
    # The integration variables above are in fm^-1.  Converting the production
    # cache formula from MeV variables gives hbarc^(8-2-5) = hbarc: four radial
    # measures contribute hbarc^8, the WP normalization hbarc^-2, and the
    # W1(fm^5)->MeV conversion hbarc^-5.
    return accum * bin_norm * HBARC


def convergence_table(cE, cD, c1, c3, Lambda, bin_bounds, label,
                      orders=(1, 2, 4, 8), **driver_options):
    """Print convergence of the bin average for successive quadrature orders."""
    print(f"\n  {label}")
    print(f"    bin: p=[{bin_bounds[0]:.3f},{bin_bounds[1]:.3f}] "
          f"q=[{bin_bounds[2]:.3f},{bin_bounds[3]:.3f}] "
          f"p'=[{bin_bounds[4]:.3f},{bin_bounds[5]:.3f}] "
          f"q'=[{bin_bounds[6]:.3f},{bin_bounds[7]:.3f}] fm^-1")
    vals = {}
    for N in orders:
        vals[N] = bin_average(
            cE, cD, c1, c3, Lambda, bin_bounds, N, **driver_options
        )
    reference_order = orders[-1]
    ref = vals[reference_order]
    print(f"    {'N':>3} {'W1_WP [MeV]':>16} "
          f"{'rel. diff vs N='+str(reference_order):>18}")
    for N in orders:
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

    print("\nComplete factorized c1/c3/cE, ordinary cell:")
    convergence_table(-0.205, -0.2, -0.81, -3.2, Lambda,
                      (0.40, 0.60, 0.30, 0.50, 0.50, 0.70, 0.40, 0.60),
                      "full factorized diagonal", orders=(1, 2, 4, 6),
                      c4=5.4, factorized=True, transfer_order=6,
                      alpha_r=4, alpha_c=4)

    print("\nComplete factorized c4/cD transition, ordinary cell:")
    convergence_table(-0.205, -0.2, -0.81, -3.2, Lambda,
                      (0.40, 0.60, 0.30, 0.50, 0.50, 0.70, 0.40, 0.60),
                      "full factorized 3S1-to-1S0", orders=(1, 2, 4),
                      c4=5.4, factorized=True, transfer_order=6,
                      alpha_r=4, alpha_c=0)

    print("\n" + "=" * 72)
    print("Interpretation:")
    print("  Each table uses its last listed order as the numerical reference.")
    print("  rel. diff < 1e-3 : that quadrature order is converged for this cell")
    print("  rel. diff > 1e-2 : that quadrature order is NOT converged; the legacy")
    print("                              Np_per_WP_W1=Nq_per_WP_W1=1 setting introduces")
    print("                              a discretization error that grows with bin width.")
    print("  Default N=2 is a safer baseline, not a universal convergence certificate;")
    print("  wide bins in this report require N=4 for 1e-3-level convergence.")
    print("=" * 72)


if __name__ == "__main__":
    run()
