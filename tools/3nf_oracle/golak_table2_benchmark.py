#!/usr/bin/env python3
"""Reproduce the simplest Golak et al. (2010) Table 2 five-angle integrals.

This is an independent published-number normalization gate for the full-vector
oracle.  It transcribes Eqs. (18), (19), (23), (25), and (26) of
Eur. Phys. J. A 43, 241 (2010), not any Tic-tac production PWD formula.

The two matrix elements implemented here isolate complementary physics:

* G(1,1): c1/c3 and tau2.tau3;
* G(2,1): c4 and the imaginary tau1.(tau2 x tau3) recoupling.

Both are five-dimensional Gauss-Legendre angular integrals.  The published
standard-PWD targets are 443.618 and 1200.219 fm^5, respectively.
"""

from __future__ import annotations

import argparse
import math
from typing import Dict, Tuple

import numpy as np


HBARC = 197.327
GA = 1.29
FPI = 92.4 / HBARC
MPI = 138.0 / HBARC
C1 = -0.81 * HBARC / 1000.0
C3 = -3.4 * HBARC / 1000.0
C4 = +3.4 * HBARC / 1000.0

P = 1.0
Q = 2.0
PP = 3.0
QP = 4.0

PUBLISHED = {"G(1,1)": 443.618, "G(2,1)": 1200.219}


def _nodes_weights(n: int, lo: float, hi: float) -> Tuple[np.ndarray, np.ndarray]:
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * ((hi - lo) * x + hi + lo), 0.5 * (hi - lo) * w


def integrate(n: int) -> Dict[str, float]:
    """Evaluate the two Eq. (25) integrals with N points in every dimension."""
    cosines, cosine_weights = _nodes_weights(n, -1.0, 1.0)
    phis, phi_weights = _nodes_weights(n, 0.0, 2.0 * math.pi)

    sin_cos = np.sqrt(np.maximum(0.0, 1.0 - cosines * cosines))
    cos_phi = np.cos(phis)
    sin_phi = np.sin(phis)

    # Primed vectors on a phi_p' x phi_q' mesh for each polar-angle pair.
    total_11 = 0.0
    total_21 = 0.0j
    phi_weight_grid = phi_weights[:, None] * phi_weights[None, :]

    for iq, (cq, sq, wq) in enumerate(zip(cosines, sin_cos, cosine_weights)):
        del iq
        p_in = np.array([0.0, 0.0, P])
        q_in = Q * np.array([sq, 0.0, cq])
        k2_in = p_in - 0.5 * q_in
        k3_in = -p_in - 0.5 * q_in

        for ipp, (cpp, spp, wpp) in enumerate(zip(cosines, sin_cos, cosine_weights)):
            del ipp
            p_out = PP * np.stack(
                (spp * cos_phi, spp * sin_phi, np.full(n, cpp)), axis=-1
            )
            for iqp, (cqp, sqp, wqp) in enumerate(zip(cosines, sin_cos, cosine_weights)):
                del iqp
                q_out = QP * np.stack(
                    (sqp * cos_phi, sqp * sin_phi, np.full(n, cqp)), axis=-1
                )

                q2 = p_out[:, None, :] - 0.5 * q_out[None, :, :] - k2_in
                q3 = -p_out[:, None, :] - 0.5 * q_out[None, :, :] - k3_in
                q2_sq = np.einsum("...i,...i->...", q2, q2)
                q3_sq = np.einsum("...i,...i->...", q3, q3)
                q2_dot_q3 = np.einsum("...i,...i->...", q2, q3)

                denominator = (q2_sq + MPI * MPI) * (q3_sq + MPI * MPI)
                axial = (GA / (2.0 * FPI)) ** 2 / denominator
                f1 = axial * (
                    -4.0 * C1 * MPI * MPI / (FPI * FPI)
                    + 2.0 * C3 * q2_dot_q3 / (FPI * FPI)
                )
                f2 = axial * C4 / (FPI * FPI)

                # Eq. (22): I1(t'=1,t=1)=+1.
                g11 = -(f1 * q2_dot_q3) / (16.0 * math.pi**2)

                # Eq. (23): I2(t'=0,t=1)=i 2 sqrt(3)(-1)^(t'+1)
                # = -i 2 sqrt(3).  The plain-text PDF extraction loses the 2
                # immediately before the radical; direct 8-state Pauli
                # enumeration independently gives -i*3.464101615... .  Keep
                # the complex factors explicit so a phase typo cannot hide in
                # a real shortcut.
                i2_01 = -2.0j * math.sqrt(3.0)
                gram = q2_dot_q3 * q2_dot_q3 - q2_sq * q3_sq
                g21 = (-1.0j / (16.0 * math.pi**2 * math.sqrt(3.0))) * f2 * i2_01 * gram

                angular_weight = wq * wpp * wqp * phi_weight_grid
                total_11 += float(np.sum(angular_weight * g11))
                total_21 += complex(np.sum(angular_weight * g21))

    # Fixing p-hat=z and phi_q=0 removes the common global SO(3) orientation
    # volume 4pi*2pi=8pi^2 from the original eight-angle integral (Sec. 4).
    rotational_volume = 8.0 * math.pi**2
    total_11 *= rotational_volume
    total_21 *= rotational_volume
    if abs(total_21.imag) > 5e-10:
        raise ArithmeticError(f"G(2,1) should be real after isospin recoupling: {total_21}")
    return {"G(1,1)": total_11, "G(2,1)": total_21.real}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--order", type=int, default=24)
    parser.add_argument("--relative-tolerance", type=float, default=8e-4)
    args = parser.parse_args()
    values = integrate(args.order)
    failed = False
    print(f"Golak 2010 Table 2 five-angle benchmark, N={args.order}")
    for name, target in PUBLISHED.items():
        value = values[name]
        relative = abs(value - target) / abs(target)
        print(f"  {name:>6}: oracle={value:12.6f}  published={target:12.6f}  rel={relative:.3e}")
        failed = failed or relative > args.relative_tolerance
    if failed:
        raise SystemExit("published Table 2 tolerance was not met")


if __name__ == "__main__":
    main()
