#!/usr/bin/env python3
"""Finite-rank spin-momentum extensions of the Hebeler scalar projector.

The implementation expands Cartesian Jacobi vectors multiplying a spherical
harmonic into the finite ``l-1,l+1`` harmonic basis.  Spin matrices and all
magnetic-projection sums are then momentum independent; the remaining radial
objects are the scalar three-integral kernels from ``factorized_scalar_pwd``.

The module exposes all momentum-dependent spectator-1 N2LO structures.  It is
kept in the independent Python oracle stack until the factorized results have
been validated term by term and the same algorithm is ported to production
C++.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import itertools
import math
from typing import Iterable

import numpy as np

from factorized_scalar_pwd import ScalarKernel, project_scalar_ls, uncoupled_orbital_kernel
from full_vector_five_angle_pwd import LSChannel, _cg, _isospin_state, _OP


@dataclass(frozen=True)
class _ChannelTerm:
    m_pair: int
    m_spectator: int
    spin_index: int
    coefficient: complex


@dataclass(frozen=True)
class _CartesianTerm:
    coefficient: complex
    coordinate_axes: tuple[tuple[str, int], ...]
    spin_matrix: np.ndarray


_CARTESIAN_Y1 = {
    0: {-1: math.sqrt(2.0 * math.pi / 3.0), 1: -math.sqrt(2.0 * math.pi / 3.0)},
    1: {-1: 1.0j * math.sqrt(2.0 * math.pi / 3.0), 1: 1.0j * math.sqrt(2.0 * math.pi / 3.0)},
    2: {0: math.sqrt(4.0 * math.pi / 3.0)},
}


@lru_cache(maxsize=None)
def _multiply_once(l_value: int, m_value: int, axis: int) -> tuple[tuple[int, int, complex], ...]:
    """Expand ``x_hat[axis] Y_lm`` in normalized spherical harmonics."""
    result: dict[tuple[int, int], complex] = {}
    for mu, cartesian_coefficient in _CARTESIAN_Y1[axis].items():
        new_m = m_value + mu
        for new_l in range(abs(l_value - 1), l_value + 2):
            if abs(new_m) > new_l:
                continue
            product_coefficient = (
                math.sqrt(3.0 / (4.0 * math.pi) * (2 * l_value + 1) / (2 * new_l + 1))
                * _cg(2 * l_value, 2, 2 * new_l, 0, 0, 0)
                * _cg(2 * l_value, 2, 2 * new_l, 2 * m_value, 2 * mu, 2 * new_m)
            )
            value = cartesian_coefficient * product_coefficient
            if abs(value) > 1.0e-15:
                result[(new_l, new_m)] = result.get((new_l, new_m), 0.0j) + value
    return tuple((l_new, m_new, value) for (l_new, m_new), value in result.items())


@lru_cache(maxsize=None)
def _multiply_axes(
    l_value: int,
    m_value: int,
    axes: tuple[int, ...],
) -> tuple[tuple[int, int, complex], ...]:
    states: dict[tuple[int, int], complex] = {(l_value, m_value): 1.0 + 0.0j}
    for axis in axes:
        updated: dict[tuple[int, int], complex] = {}
        for (current_l, current_m), current_coefficient in states.items():
            for new_l, new_m, factor in _multiply_once(current_l, current_m, axis):
                key = (new_l, new_m)
                updated[key] = updated.get(key, 0.0j) + current_coefficient * factor
        states = updated
    return tuple((l_new, m_new, value) for (l_new, m_new), value in states.items())


@lru_cache(maxsize=None)
def _channel_terms(channel: LSChannel, two_m_j: int) -> tuple[_ChannelTerm, ...]:
    accumulated: dict[tuple[int, int, int], complex] = {}
    for m_pair in range(-channel.l_pair, channel.l_pair + 1):
        for m_spectator in range(
            -channel.lambda_spectator, channel.lambda_spectator + 1
        ):
            m_total_l = m_pair + m_spectator
            orbital_cg = _cg(
                2 * channel.l_pair,
                2 * channel.lambda_spectator,
                2 * channel.total_L,
                2 * m_pair,
                2 * m_spectator,
                2 * m_total_l,
            )
            if orbital_cg == 0.0:
                continue
            for two_m_pair_spin in range(-2 * channel.s_pair, 2 * channel.s_pair + 1, 2):
                for two_m_one in (-1, 1):
                    two_m_s = two_m_pair_spin + two_m_one
                    spin_cg = _cg(
                        2 * channel.s_pair,
                        1,
                        channel.two_total_S,
                        two_m_pair_spin,
                        two_m_one,
                        two_m_s,
                    )
                    total_cg = _cg(
                        2 * channel.total_L,
                        channel.two_total_S,
                        channel.two_total_J,
                        2 * m_total_l,
                        two_m_s,
                        two_m_j,
                    )
                    if spin_cg == 0.0 or total_cg == 0.0:
                        continue
                    for two_m_two in (-1, 1):
                        for two_m_three in (-1, 1):
                            pair_cg = _cg(
                                1,
                                1,
                                2 * channel.s_pair,
                                two_m_two,
                                two_m_three,
                                two_m_pair_spin,
                            )
                            if pair_cg == 0.0:
                                continue

                            def bit(two_m: int) -> int:
                                return 0 if two_m == 1 else 1

                            spin_index = (
                                4 * bit(two_m_one)
                                + 2 * bit(two_m_two)
                                + bit(two_m_three)
                            )
                            key = (m_pair, m_spectator, spin_index)
                            accumulated[key] = accumulated.get(key, 0.0j) + (
                                orbital_cg * spin_cg * total_cg * pair_cg
                            )
    return tuple(
        _ChannelTerm(m_pair, m_spectator, spin_index, coefficient)
        for (m_pair, m_spectator, spin_index), coefficient in accumulated.items()
        if abs(coefficient) > 1.0e-15
    )


def _coordinate_slot(coordinate: str) -> tuple[str, str]:
    return {
        "p_bra": ("bra", "pair"),
        "q_bra": ("bra", "spectator"),
        "p_ket": ("ket", "pair"),
        "q_ket": ("ket", "spectator"),
    }[coordinate]


def _coordinate_magnitude(
    coordinate: str,
    momenta: tuple[float, float, float, float],
) -> float:
    p, q, pp, qp = momenta
    return {"p_bra": pp, "q_bra": qp, "p_ket": p, "q_ket": q}[coordinate]


def _transformed_slots(
    channel: LSChannel,
    term: _ChannelTerm,
    axes_by_slot: dict[tuple[str, str], tuple[int, ...]],
    side: str,
) -> tuple[
    tuple[tuple[int, int, complex], ...],
    tuple[tuple[int, int, complex], ...],
]:
    pair = _multiply_axes(
        channel.l_pair,
        term.m_pair,
        axes_by_slot.get((side, "pair"), ()),
    )
    spectator = _multiply_axes(
        channel.lambda_spectator,
        term.m_spectator,
        axes_by_slot.get((side, "spectator"), ()),
    )
    if side == "bra":
        pair = tuple((l_value, m_value, np.conj(value)) for l_value, m_value, value in pair)
        spectator = tuple(
            (l_value, m_value, np.conj(value))
            for l_value, m_value, value in spectator
        )
    return pair, spectator


def _project_cartesian_terms(
    bra: LSChannel,
    ket: LSChannel,
    momenta: tuple[float, float, float, float],
    scalar_kernel: ScalarKernel,
    cartesian_terms: Iterable[_CartesianTerm],
    order: int,
) -> complex:
    if bra.two_total_J != ket.two_total_J:
        return 0.0j

    orbital_weights: dict[tuple[int, ...], complex] = {}
    m_count = bra.two_total_J + 1
    for two_m_j in range(-bra.two_total_J, bra.two_total_J + 1, 2):
        bra_terms = _channel_terms(bra, two_m_j)
        ket_terms = _channel_terms(ket, two_m_j)
        for cartesian in cartesian_terms:
            axes_by_slot: dict[tuple[str, str], tuple[int, ...]] = {}
            radial_coefficient = cartesian.coefficient
            for coordinate, axis in cartesian.coordinate_axes:
                slot = _coordinate_slot(coordinate)
                axes_by_slot[slot] = axes_by_slot.get(slot, ()) + (axis,)
                radial_coefficient *= _coordinate_magnitude(coordinate, momenta)

            for bra_term in bra_terms:
                bra_pair, bra_spectator = _transformed_slots(
                    bra, bra_term, axes_by_slot, "bra"
                )
                for ket_term in ket_terms:
                    spin_element = cartesian.spin_matrix[
                        bra_term.spin_index, ket_term.spin_index
                    ]
                    if abs(spin_element) < 1.0e-15:
                        continue
                    ket_pair, ket_spectator = _transformed_slots(
                        ket, ket_term, axes_by_slot, "ket"
                    )
                    state_coefficient = (
                        radial_coefficient
                        * np.conj(bra_term.coefficient)
                        * ket_term.coefficient
                        * spin_element
                        / m_count
                    )
                    for (
                        (bra_l_pair, bra_m_pair, bra_pair_coefficient),
                        (bra_lambda, bra_m_lambda, bra_spectator_coefficient),
                        (ket_l_pair, ket_m_pair, ket_pair_coefficient),
                        (ket_lambda, ket_m_lambda, ket_spectator_coefficient),
                    ) in itertools.product(
                        bra_pair, bra_spectator, ket_pair, ket_spectator
                    ):
                        key = (
                            bra_l_pair,
                            bra_m_pair,
                            bra_lambda,
                            bra_m_lambda,
                            ket_l_pair,
                            ket_m_pair,
                            ket_lambda,
                            ket_m_lambda,
                        )
                        orbital_weights[key] = orbital_weights.get(key, 0.0j) + (
                            state_coefficient
                            * bra_pair_coefficient
                            * bra_spectator_coefficient
                            * ket_pair_coefficient
                            * ket_spectator_coefficient
                        )

    result = 0.0j
    orbital_cache: dict[tuple[int, ...], complex] = {}
    for key, coefficient in orbital_weights.items():
        if abs(coefficient) < 2.0e-14:
            continue
        if key not in orbital_cache:
            orbital_cache[key] = uncoupled_orbital_kernel(
                *key,
                momenta,
                scalar_kernel,
                order,
            )
        result += coefficient * orbital_cache[key]
    return result


@lru_cache(maxsize=None)
def _q_transfer_dot_terms(second_particle: int) -> tuple[_CartesianTerm, ...]:
    """Expand ``(sigma1.Deltaq)(sigma_i.Deltaq)`` in q/q' vectors."""
    sources = (("q_bra", 1.0), ("q_ket", -1.0))
    terms = []
    for (first_coordinate, first_sign), (second_coordinate, second_sign) in itertools.product(
        sources, sources
    ):
        for first_axis in range(3):
            for second_axis in range(3):
                terms.append(_CartesianTerm(
                    first_sign * second_sign,
                    (
                        (first_coordinate, first_axis),
                        (second_coordinate, second_axis),
                    ),
                    _OP._SIGMA[0][first_axis]
                    @ _OP._SIGMA[second_particle - 1][second_axis],
                ))
    return tuple(terms)


@lru_cache(maxsize=1)
def _q2_q3_dot_terms() -> tuple[_CartesianTerm, ...]:
    """Expand ``(sigma2.Q2)(sigma3.Q3)`` in p,p',q,q' vectors."""
    q2_sources = (
        ("p_bra", +1.0),
        ("p_ket", -1.0),
        ("q_bra", -0.5),
        ("q_ket", +0.5),
    )
    q3_sources = (
        ("p_bra", -1.0),
        ("p_ket", +1.0),
        ("q_bra", -0.5),
        ("q_ket", +0.5),
    )
    terms = []
    for (q2_coordinate, q2_sign), (q3_coordinate, q3_sign) in itertools.product(
        q2_sources, q3_sources
    ):
        for q2_axis in range(3):
            for q3_axis in range(3):
                terms.append(_CartesianTerm(
                    q2_sign * q3_sign,
                    (
                        (q2_coordinate, q2_axis),
                        (q3_coordinate, q3_axis),
                    ),
                    _OP._SIGMA[1][q2_axis] @ _OP._SIGMA[2][q3_axis],
                ))
    return tuple(terms)


def _levi_civita_terms() -> tuple[tuple[int, int, int, int], ...]:
    return (
        (0, 1, 2, +1),
        (1, 2, 0, +1),
        (2, 0, 1, +1),
        (1, 0, 2, -1),
        (2, 1, 0, -1),
        (0, 2, 1, -1),
    )


@lru_cache(maxsize=1)
def _c4_cartesian_terms() -> tuple[_CartesianTerm, ...]:
    """Expand ``(sigma2.Q2)(sigma3.Q3)sigma1.(Q2 x Q3)``.

    In spectator-1 Jacobi variables ``Q2 x Q3 = -Delta p x Delta q``.
    """
    q2_sources = (
        ("p_bra", +1.0),
        ("p_ket", -1.0),
        ("q_bra", -0.5),
        ("q_ket", +0.5),
    )
    q3_sources = (
        ("p_bra", -1.0),
        ("p_ket", +1.0),
        ("q_bra", -0.5),
        ("q_ket", +0.5),
    )
    delta_p_sources = (("p_bra", +1.0), ("p_ket", -1.0))
    delta_q_sources = (("q_bra", +1.0), ("q_ket", -1.0))
    terms = []
    for (
        (q2_coordinate, q2_sign),
        (q3_coordinate, q3_sign),
        (delta_p_coordinate, delta_p_sign),
        (delta_q_coordinate, delta_q_sign),
    ) in itertools.product(
        q2_sources,
        q3_sources,
        delta_p_sources,
        delta_q_sources,
    ):
        for q2_axis in range(3):
            for q3_axis in range(3):
                for delta_p_axis, delta_q_axis, spin_axis, epsilon in _levi_civita_terms():
                    terms.append(_CartesianTerm(
                        -q2_sign
                        * q3_sign
                        * delta_p_sign
                        * delta_q_sign
                        * epsilon,
                        (
                            (q2_coordinate, q2_axis),
                            (q3_coordinate, q3_axis),
                            (delta_p_coordinate, delta_p_axis),
                            (delta_q_coordinate, delta_q_axis),
                        ),
                        _OP._SIGMA[1][q2_axis]
                        @ _OP._SIGMA[2][q3_axis]
                        @ _OP._SIGMA[0][spin_axis],
                    ))
    return tuple(terms)


def project_c1_c3_ls(
    bra: LSChannel,
    ket: LSChannel,
    momenta: tuple[float, float, float, float],
    constants: _OP.N2LOConstants,
    c1_gev_inverse: float,
    c3_gev_inverse: float,
    order: int,
    two_m_t: int = 1,
) -> dict[str, complex]:
    """Raw complete spectator-1 ``c1`` and ``c3`` LS matrix elements."""
    iso_bra = _isospin_state(bra, two_m_t)
    iso_ket = _isospin_state(ket, two_m_t)
    iso23 = np.vdot(iso_bra, _OP._TAU23 @ iso_ket)
    m_pi = constants.m_pi

    def propagators(delta_p: float, delta_q: float, cosine: float) -> float:
        common = delta_p**2 + 0.25 * delta_q**2 + m_pi**2
        d2 = common - delta_p * delta_q * cosine
        d3 = common + delta_p * delta_q * cosine
        return 1.0 / (d2 * d3)

    def c3_kernel(delta_p: float, delta_q: float, cosine: float) -> float:
        q2_dot_q3 = -delta_p**2 + 0.25 * delta_q**2
        return q2_dot_q3 * propagators(delta_p, delta_q, cosine)

    spin_c1 = _project_cartesian_terms(
        bra, ket, momenta, propagators, _q2_q3_dot_terms(), order
    )
    spin_c3 = _project_cartesian_terms(
        bra, ket, momenta, c3_kernel, _q2_q3_dot_terms(), order
    )
    common = constants.g_a**2 / (4.0 * constants.f_pi**4)
    c1_fm = constants.gev_inverse_to_fm(c1_gev_inverse)
    c3_fm = constants.gev_inverse_to_fm(c3_gev_inverse)
    return {
        "c1": common * (-4.0 * c1_fm * m_pi**2) * iso23 * spin_c1,
        "c3": common * (2.0 * c3_fm) * iso23 * spin_c3,
    }


def project_c4_ls(
    bra: LSChannel,
    ket: LSChannel,
    momenta: tuple[float, float, float, float],
    constants: _OP.N2LOConstants,
    c4_gev_inverse: float,
    order: int,
    two_m_t: int = 1,
) -> complex:
    """Raw complete spectator-1 ``c4`` LS matrix element."""
    iso_bra = _isospin_state(bra, two_m_t)
    iso_ket = _isospin_state(ket, two_m_t)
    iso_cross = np.vdot(
        iso_bra, _OP._TAU1_DOT_TAU2_CROSS_TAU3 @ iso_ket
    )
    m_pi = constants.m_pi

    def propagators(delta_p: float, delta_q: float, cosine: float) -> float:
        common = delta_p**2 + 0.25 * delta_q**2 + m_pi**2
        d2 = common - delta_p * delta_q * cosine
        d3 = common + delta_p * delta_q * cosine
        return 1.0 / (d2 * d3)

    spin_c4 = _project_cartesian_terms(
        bra, ket, momenta, propagators, _c4_cartesian_terms(), order
    )
    c4_fm = constants.gev_inverse_to_fm(c4_gev_inverse)
    return (
        constants.g_a**2
        / (4.0 * constants.f_pi**4)
        * c4_fm
        * iso_cross
        * spin_c4
    )


def project_cD_ls(
    bra: LSChannel,
    ket: LSChannel,
    momenta: tuple[float, float, float, float],
    constants: _OP.N2LOConstants,
    c_d: float,
    order: int,
    two_m_t: int = 1,
) -> complex:
    """Raw complete spectator-1 ``cD`` matrix element in the LS basis."""
    iso_bra = _isospin_state(bra, two_m_t)
    iso_ket = _isospin_state(ket, two_m_t)
    iso12 = np.vdot(iso_bra, _OP._TAU12 @ iso_ket)
    iso13 = np.vdot(iso_bra, _OP._TAU13 @ iso_ket)
    m_pi = constants.m_pi

    def pion_kernel(_delta_p: float, delta_q: float, _cosine: float) -> float:
        return 1.0 / (delta_q * delta_q + m_pi * m_pi)

    spin12 = _project_cartesian_terms(
        bra,
        ket,
        momenta,
        pion_kernel,
        _q_transfer_dot_terms(2),
        order,
    )
    spin13 = _project_cartesian_terms(
        bra,
        ket,
        momenta,
        pion_kernel,
        _q_transfer_dot_terms(3),
        order,
    )
    d_lec = c_d / (constants.f_pi**2 * constants.lambda_chi)
    return (
        -constants.g_a
        * d_lec
        / (8.0 * constants.f_pi**2)
        * (iso12 * spin12 + iso13 * spin13)
    )


def project_cE_ls(
    bra: LSChannel,
    ket: LSChannel,
    momenta: tuple[float, float, float, float],
    constants: _OP.N2LOConstants,
    c_e: float,
    order: int,
    two_m_t: int = 1,
) -> complex:
    """Raw spectator-1 contact matrix element before ``(2*pi)^-6``."""
    iso_bra = _isospin_state(bra, two_m_t)
    iso_ket = _isospin_state(ket, two_m_t)
    iso23 = np.vdot(iso_bra, _OP._TAU23 @ iso_ket)
    scalar = project_scalar_ls(
        bra,
        ket,
        momenta,
        lambda _delta_p, _delta_q, _cosine: 1.0,
        order,
    )
    e_lec = c_e / (constants.f_pi**4 * constants.lambda_chi)
    return e_lec * iso23 * scalar
