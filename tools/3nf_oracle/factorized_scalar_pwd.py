#!/usr/bin/env python3
"""Hebeler three-integral projector for a spin-independent local 3NF kernel.

This module is the first, deliberately small building block of the scalable
N2LO projector.  It transcribes Hebeler et al., Phys. Rev. C 91, 044001
(2015), Eq. (6), in the pair-23/spectator-1 Jacobi convention used by
Tic-tac.  The formula performs the five angular integrations analytically and
leaves integrations over the two transfer magnitudes and their relative
angle.

The returned values are *raw angular partial-wave matrix elements*.  The
Tic-tac ``(2*pi)^-6`` Jacobi Fourier normalization is intentionally not
applied here, so the constant-kernel diagonal limit is ``(4*pi)^2``.  Keeping
that boundary explicit prevents the factorized projector from hiding a
normalization convention inside its numerical quadrature.
"""

from __future__ import annotations

from functools import lru_cache
import math
from typing import Callable

import numpy as np

from full_vector_five_angle_pwd import LSChannel, _cg, bipolar_harmonic


ScalarKernel = Callable[[float, float, float], complex]


def _nodes_weights(order: int, lo: float, hi: float) -> tuple[np.ndarray, np.ndarray]:
    if order < 1:
        raise ValueError("quadrature order must be positive")
    nodes, weights = np.polynomial.legendre.leggauss(order)
    return 0.5 * ((hi - lo) * nodes + hi + lo), 0.5 * (hi - lo) * weights


@lru_cache(maxsize=None)
def _orbital_cg(
    left_l: int,
    left_m: int,
    right_l: int,
    right_m: int,
    total_l: int,
    total_m: int,
) -> float:
    return _cg(
        2 * left_l,
        2 * right_l,
        2 * total_l,
        2 * left_m,
        2 * right_m,
        2 * total_m,
    )


def uncoupled_orbital_kernel(
    bra_l_pair: int,
    bra_m_pair: int,
    bra_lambda: int,
    bra_m_lambda: int,
    ket_l_pair: int,
    ket_m_pair: int,
    ket_lambda: int,
    ket_m_lambda: int,
    momenta: tuple[float, float, float, float],
    scalar_kernel: ScalarKernel,
    order: int,
) -> complex:
    """Evaluate Hebeler Eq. (6) for uncoupled orbital projections.

    ``momenta`` are ``(p, q, p_prime, q_prime)``.  ``scalar_kernel`` receives
    ``(delta_p, delta_q, cos_theta_delta_p_delta_q)`` in that order.
    """

    p, q, pp, qp = map(float, momenta)
    if min(p, q, pp, qp) <= 0.0:
        raise ValueError("Eq. (6) requires strictly positive external momenta")

    # Hebeler's magnetic-projection selection rule, Eq. (6), first line.
    if ket_m_pair - bra_m_pair != bra_m_lambda - ket_m_lambda:
        return 0.0j

    lbar_min = max(
        abs(bra_l_pair - ket_l_pair),
        abs(bra_lambda - ket_lambda),
    )
    lbar_max = min(
        bra_l_pair + ket_l_pair,
        bra_lambda + ket_lambda,
    )
    if lbar_min > lbar_max:
        return 0.0j

    delta_p_nodes, delta_p_weights = _nodes_weights(order, abs(pp - p), pp + p)
    delta_q_nodes, delta_q_weights = _nodes_weights(order, abs(qp - q), qp + q)
    x_nodes, x_weights = _nodes_weights(order, -1.0, 1.0)

    phase = -1.0 if (ket_m_pair + bra_m_lambda) % 2 else 1.0
    prefactor = phase * 2.0 * (2.0 * math.pi) ** 4 / (p * pp * q * qp)
    total = 0.0j

    for lbar in range(lbar_min, lbar_max + 1):
        pair_cg = _orbital_cg(
            bra_l_pair,
            -bra_m_pair,
            ket_l_pair,
            ket_m_pair,
            lbar,
            -bra_m_pair + ket_m_pair,
        )
        spectator_cg = _orbital_cg(
            bra_lambda,
            -bra_m_lambda,
            ket_lambda,
            ket_m_lambda,
            lbar,
            -bra_m_lambda + ket_m_lambda,
        )
        angular_coefficient = pair_cg * spectator_cg / (2 * lbar + 1)
        if angular_coefficient == 0.0:
            continue

        lbar_integral = 0.0j
        for delta_p, weight_p in zip(delta_p_nodes, delta_p_weights):
            cosine_p = (pp * pp - p * p - delta_p * delta_p) / (2.0 * delta_p * p)
            cosine_p = float(np.clip(cosine_p, -1.0, 1.0))
            sine_p = math.sqrt(max(0.0, 1.0 - cosine_p * cosine_p))
            p_hat = np.array([sine_p, 0.0, cosine_p])
            pp_vector = p * p_hat + np.array([0.0, 0.0, delta_p])
            pp_hat = pp_vector / np.linalg.norm(pp_vector)
            pair_bipolar = bipolar_harmonic(
                bra_l_pair,
                ket_l_pair,
                lbar,
                0,
                pp_hat,
                p_hat,
            )

            for delta_q, weight_q in zip(delta_q_nodes, delta_q_weights):
                cosine_q = (qp * qp - q * q - delta_q * delta_q) / (2.0 * delta_q * q)
                cosine_q = float(np.clip(cosine_q, -1.0, 1.0))
                sine_q = math.sqrt(max(0.0, 1.0 - cosine_q * cosine_q))
                q_hat = np.array([sine_q, 0.0, cosine_q])
                qp_vector = q * q_hat + np.array([0.0, 0.0, delta_q])
                qp_hat = qp_vector / np.linalg.norm(qp_vector)
                spectator_bipolar = bipolar_harmonic(
                    bra_lambda,
                    ket_lambda,
                    lbar,
                    0,
                    qp_hat,
                    q_hat,
                )

                relative_integral = 0.0j
                for x, weight_x in zip(x_nodes, x_weights):
                    relative_integral += (
                        weight_x
                        * np.polynomial.legendre.legval(x, [0.0] * lbar + [1.0])
                        * scalar_kernel(float(delta_p), float(delta_q), float(x))
                    )

                lbar_integral += (
                    weight_p
                    * delta_p
                    * weight_q
                    * delta_q
                    * pair_bipolar
                    * spectator_bipolar
                    * relative_integral
                )

        total += angular_coefficient * lbar_integral

    return prefactor * total


def project_scalar_ls(
    bra: LSChannel,
    ket: LSChannel,
    momenta: tuple[float, float, float, float],
    scalar_kernel: ScalarKernel,
    order: int,
) -> complex:
    """Project a spin-independent scalar local kernel between LS channels."""

    if bra.two_total_J != ket.two_total_J:
        return 0.0j
    if bra.s_pair != ket.s_pair or bra.two_total_S != ket.two_total_S:
        return 0.0j

    two_j = bra.two_total_J
    two_s = bra.two_total_S
    total = 0.0j

    for two_m_j in range(-two_j, two_j + 1, 2):
        for two_m_s in range(-two_s, two_s + 1, 2):
            difference = two_m_j - two_m_s
            if difference % 2:
                continue
            m_total_l = difference // 2
            if abs(m_total_l) > bra.total_L or abs(m_total_l) > ket.total_L:
                continue
            bra_outer = _cg(
                2 * bra.total_L,
                two_s,
                two_j,
                2 * m_total_l,
                two_m_s,
                two_m_j,
            )
            ket_outer = _cg(
                2 * ket.total_L,
                two_s,
                two_j,
                2 * m_total_l,
                two_m_s,
                two_m_j,
            )
            if bra_outer == 0.0 or ket_outer == 0.0:
                continue

            for bra_m_pair in range(-bra.l_pair, bra.l_pair + 1):
                bra_m_lambda = m_total_l - bra_m_pair
                if abs(bra_m_lambda) > bra.lambda_spectator:
                    continue
                bra_orbital = _orbital_cg(
                    bra.l_pair,
                    bra_m_pair,
                    bra.lambda_spectator,
                    bra_m_lambda,
                    bra.total_L,
                    m_total_l,
                )
                if bra_orbital == 0.0:
                    continue

                for ket_m_pair in range(-ket.l_pair, ket.l_pair + 1):
                    ket_m_lambda = m_total_l - ket_m_pair
                    if abs(ket_m_lambda) > ket.lambda_spectator:
                        continue
                    ket_orbital = _orbital_cg(
                        ket.l_pair,
                        ket_m_pair,
                        ket.lambda_spectator,
                        ket_m_lambda,
                        ket.total_L,
                        m_total_l,
                    )
                    if ket_orbital == 0.0:
                        continue

                    total += (
                        bra_outer
                        * ket_outer
                        * bra_orbital
                        * ket_orbital
                        * uncoupled_orbital_kernel(
                            bra.l_pair,
                            bra_m_pair,
                            bra.lambda_spectator,
                            bra_m_lambda,
                            ket.l_pair,
                            ket_m_pair,
                            ket.lambda_spectator,
                            ket_m_lambda,
                            momenta,
                            scalar_kernel,
                            order,
                        )
                    )

    return total / (two_j + 1.0)

