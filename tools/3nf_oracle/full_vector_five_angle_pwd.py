#!/usr/bin/env python3
"""Generic five-angle LS-coupled projection of the full-vector N2LO oracle.

The implementation follows Golak et al., Eur. Phys. J. A 43, 241 (2010),
Eqs. (12)--(16): ``p_hat`` is fixed to z, ``phi_q`` to zero, and the removed
global orientation contributes ``8*pi**2``.  Spin and isospin states are built
directly from Clebsch--Gordan coefficients in an uncoupled three-particle
m-scheme basis.  No Tic-tac production recoupling routine is used.

This is deliberately a transparent reference projector, not a fast production
algorithm.  It supports every c1/c3/c4/cD/cE spin-momentum structure and all LS
channels allowed by the supplied quantum numbers.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import importlib.util
import math
from pathlib import Path
import sys
from typing import Dict, Iterable, Tuple

import numpy as np
from scipy.special import sph_harm
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan, wigner_9j


_THIS_DIR = Path(__file__).resolve().parent
_OPERATOR_PATH = _THIS_DIR / "full_vector_n2lo_oracle.py"
_SPEC = importlib.util.spec_from_file_location("full_vector_n2lo_operator", _OPERATOR_PATH)
_OP = importlib.util.module_from_spec(_SPEC)
assert _SPEC.loader is not None
sys.modules[_SPEC.name] = _OP
_SPEC.loader.exec_module(_OP)


# Tic-tac uses one (2*pi)^-3 Fourier factor for each independent relative
# coordinate.  A three-body Jacobi state has the two coordinates p and q, so a
# raw angular PWD such as Golak's G(beta',beta) must be multiplied by this
# factor before it is consumed by Tic-tac's radial/WP measure.
THREE_BODY_FOURIER_NORMALIZATION = (2.0 * math.pi) ** -6


@lru_cache(maxsize=None)
def _cg(
    two_j1: int,
    two_j2: int,
    two_j: int,
    two_m1: int,
    two_m2: int,
    two_m: int,
) -> float:
    if two_m1 + two_m2 != two_m:
        return 0.0
    return float(
        clebsch_gordan(
            Rational(two_j1, 2),
            Rational(two_j2, 2),
            Rational(two_j, 2),
            Rational(two_m1, 2),
            Rational(two_m2, 2),
            Rational(two_m, 2),
        )
    )


def _m_values(two_j: int) -> range:
    return range(-two_j, two_j + 1, 2)


@lru_cache(maxsize=None)
def _w9j(two_args: Tuple[int, ...]) -> float:
    if len(two_args) != 9:
        raise ValueError("a Wigner 9j symbol needs nine arguments")
    return float(wigner_9j(*(Rational(value, 2) for value in two_args)))


def _product_index(two_m1: int, two_m2: int, two_m3: int) -> int:
    def bit(two_m: int) -> int:
        if two_m == 1:
            return 0
        if two_m == -1:
            return 1
        raise ValueError("one-particle spin/isospin projection must be +/-1/2")
    return 4 * bit(two_m1) + 2 * bit(two_m2) + bit(two_m3)


@lru_cache(maxsize=None)
def coupled_three_half_state(pair_angular_momentum: int, two_total: int, two_m_total: int) -> np.ndarray:
    """Return |((1/2,1/2) pair,1/2) total M> in particle order (1,2,3)."""
    two_pair = 2 * pair_angular_momentum
    result = np.zeros(8, dtype=complex)
    for two_m_pair in _m_values(two_pair):
        for two_m1 in (-1, 1):
            outer = _cg(two_pair, 1, two_total, two_m_pair, two_m1, two_m_total)
            if outer == 0.0:
                continue
            for two_m2 in (-1, 1):
                for two_m3 in (-1, 1):
                    inner = _cg(1, 1, two_pair, two_m2, two_m3, two_m_pair)
                    if inner != 0.0:
                        result[_product_index(two_m1, two_m2, two_m3)] += outer * inner
    norm = np.linalg.norm(result)
    if not np.isclose(norm, 1.0, atol=2e-14):
        raise ArithmeticError(f"coupled three-half state has norm {norm}")
    return result


@dataclass(frozen=True)
class LSChannel:
    """Golak Eq. (12) channel (l, s, lambda, L, S, J, t, T)."""

    l_pair: int
    s_pair: int
    lambda_spectator: int
    total_L: int
    two_total_S: int
    two_total_J: int
    t_pair: int
    two_total_T: int = 1

    def __post_init__(self) -> None:
        if (self.l_pair + self.s_pair + self.t_pair) % 2 != 1:
            raise ValueError("channel violates pair-23 antisymmetry: (-1)^(l+s+t) must be -1")
        if abs(2 * self.total_L - self.two_total_S) > self.two_total_J:
            raise ValueError("L and S cannot couple to the requested J")
        if 2 * self.total_L + self.two_total_S < self.two_total_J:
            raise ValueError("L and S cannot couple to the requested J")


@dataclass(frozen=True)
class JjChannel:
    """Tic-tac/Miller Jj channel (l,s,j,lambda,I,J,t,T)."""

    l_pair: int
    s_pair: int
    j_pair: int
    lambda_spectator: int
    two_j_spectator: int
    two_total_J: int
    t_pair: int
    two_total_T: int = 1

    def __post_init__(self) -> None:
        if (self.l_pair + self.s_pair + self.t_pair) % 2 != 1:
            raise ValueError("channel violates pair-23 antisymmetry: (-1)^(l+s+t) must be -1")
        if not abs(self.l_pair - self.s_pair) <= self.j_pair <= self.l_pair + self.s_pair:
            raise ValueError("pair l and s cannot couple to j_pair")
        if self.two_j_spectator not in (
            abs(2 * self.lambda_spectator - 1),
            2 * self.lambda_spectator + 1,
        ):
            raise ValueError("spectator lambda and spin 1/2 cannot couple to j_spectator")
        if abs(2 * self.j_pair - self.two_j_spectator) > self.two_total_J:
            raise ValueError("j_pair and j_spectator cannot couple to total J")
        if 2 * self.j_pair + self.two_j_spectator < self.two_total_J:
            raise ValueError("j_pair and j_spectator cannot couple to total J")

    def ls_expansion(self) -> Tuple[Tuple[LSChannel, float], ...]:
        """Return |Jj> = sum_(L,S) coefficient |LS> via the unitary 9j."""
        result = []
        for total_l in range(
            abs(self.l_pair - self.lambda_spectator),
            self.l_pair + self.lambda_spectator + 1,
        ):
            for two_total_s in range(abs(2 * self.s_pair - 1), 2 * self.s_pair + 2, 2):
                if abs(2 * total_l - two_total_s) > self.two_total_J:
                    continue
                if 2 * total_l + two_total_s < self.two_total_J:
                    continue
                coefficient = math.sqrt(
                    (2 * self.j_pair + 1)
                    * (self.two_j_spectator + 1)
                    * (2 * total_l + 1)
                    * (two_total_s + 1)
                ) * _w9j((
                    2 * self.l_pair, 2 * self.s_pair, 2 * self.j_pair,
                    2 * self.lambda_spectator, 1, self.two_j_spectator,
                    2 * total_l, two_total_s, self.two_total_J,
                ))
                if abs(coefficient) > 2e-15:
                    result.append((LSChannel(
                        self.l_pair,
                        self.s_pair,
                        self.lambda_spectator,
                        total_l,
                        two_total_s,
                        self.two_total_J,
                        self.t_pair,
                        self.two_total_T,
                    ), coefficient))
        norm = sum(coefficient * coefficient for _, coefficient in result)
        if not math.isclose(norm, 1.0, abs_tol=2e-13):
            raise ArithmeticError(f"Jj-to-LS expansion has norm {norm}")
        return tuple(result)


def _angles(vector: np.ndarray) -> Tuple[float, float]:
    radius = float(np.linalg.norm(vector))
    if radius == 0.0:
        raise ValueError("direction vector must be nonzero")
    theta = math.acos(float(np.clip(vector[2] / radius, -1.0, 1.0)))
    phi = math.atan2(float(vector[1]), float(vector[0]))
    return theta, phi


def bipolar_harmonic(
    l_pair: int,
    lambda_spectator: int,
    total_L: int,
    m_L: int,
    p_direction: np.ndarray,
    q_direction: np.ndarray,
) -> complex:
    theta_p, phi_p = _angles(p_direction)
    theta_q, phi_q = _angles(q_direction)
    value = 0.0j
    for m_l in range(-l_pair, l_pair + 1):
        m_lambda = m_L - m_l
        if abs(m_lambda) > lambda_spectator:
            continue
        coefficient = _cg(
            2 * l_pair,
            2 * lambda_spectator,
            2 * total_L,
            2 * m_l,
            2 * m_lambda,
            2 * m_L,
        )
        value += coefficient * sph_harm(m_l, l_pair, phi_p, theta_p) * sph_harm(
            m_lambda, lambda_spectator, phi_q, theta_q
        )
    return complex(value)


def angular_spin_state(
    channel: LSChannel,
    p_direction: np.ndarray,
    q_direction: np.ndarray,
    two_m_j: int,
) -> np.ndarray:
    """Spin vector including the two coupled spherical harmonics for fixed M_J."""
    result = np.zeros(8, dtype=complex)
    for m_L in range(-channel.total_L, channel.total_L + 1):
        two_m_s = two_m_j - 2 * m_L
        if abs(two_m_s) > channel.two_total_S or (two_m_s - channel.two_total_S) % 2:
            continue
        outer = _cg(
            2 * channel.total_L,
            channel.two_total_S,
            channel.two_total_J,
            2 * m_L,
            two_m_s,
            two_m_j,
        )
        if outer == 0.0:
            continue
        orbital = bipolar_harmonic(
            channel.l_pair,
            channel.lambda_spectator,
            channel.total_L,
            m_L,
            p_direction,
            q_direction,
        )
        spin = coupled_three_half_state(channel.s_pair, channel.two_total_S, two_m_s)
        result += outer * orbital * spin
    return result


def _pair_spectator_spin_product(
    s_pair: int,
    two_m_pair: int,
    two_m_spectator: int,
) -> np.ndarray:
    """Return |m1> tensor |(1/2,1/2)s_pair,m_pair> in particle order 1,2,3."""
    result = np.zeros(8, dtype=complex)
    for two_m2 in (-1, 1):
        for two_m3 in (-1, 1):
            inner = _cg(1, 1, 2 * s_pair, two_m2, two_m3, two_m_pair)
            if inner != 0.0:
                result[_product_index(two_m_spectator, two_m2, two_m3)] += inner
    return result


def angular_spin_state_jj(
    channel: JjChannel,
    p_direction: np.ndarray,
    q_direction: np.ndarray,
    two_m_j: int,
) -> np.ndarray:
    """Direct Miller/Tic-tac Jj-coupled angular-spin state for fixed M_J."""
    theta_p, phi_p = _angles(p_direction)
    theta_q, phi_q = _angles(q_direction)
    result = np.zeros(8, dtype=complex)
    for m_l in range(-channel.l_pair, channel.l_pair + 1):
        y_l = sph_harm(m_l, channel.l_pair, phi_p, theta_p)
        for two_m_pair in _m_values(2 * channel.s_pair):
            two_m_j_pair = 2 * m_l + two_m_pair
            if abs(two_m_j_pair) > 2 * channel.j_pair:
                continue
            pair_cg = _cg(
                2 * channel.l_pair,
                2 * channel.s_pair,
                2 * channel.j_pair,
                2 * m_l,
                two_m_pair,
                two_m_j_pair,
            )
            if pair_cg == 0.0:
                continue
            for m_lambda in range(-channel.lambda_spectator, channel.lambda_spectator + 1):
                y_lambda = sph_harm(
                    m_lambda,
                    channel.lambda_spectator,
                    phi_q,
                    theta_q,
                )
                for two_m_spectator in (-1, 1):
                    two_m_j_spectator = 2 * m_lambda + two_m_spectator
                    spectator_cg = _cg(
                        2 * channel.lambda_spectator,
                        1,
                        channel.two_j_spectator,
                        2 * m_lambda,
                        two_m_spectator,
                        two_m_j_spectator,
                    )
                    total_cg = _cg(
                        2 * channel.j_pair,
                        channel.two_j_spectator,
                        channel.two_total_J,
                        two_m_j_pair,
                        two_m_j_spectator,
                        two_m_j,
                    )
                    coefficient = pair_cg * spectator_cg * total_cg
                    if coefficient != 0.0:
                        result += (
                            coefficient
                            * y_l
                            * y_lambda
                            * _pair_spectator_spin_product(
                                channel.s_pair,
                                two_m_pair,
                                two_m_spectator,
                            )
                        )
    return result


def _isospin_state(channel: LSChannel, two_m_t: int) -> np.ndarray:
    return coupled_three_half_state(channel.t_pair, channel.two_total_T, two_m_t)


def _nodes_weights(order: int, lo: float, hi: float) -> Tuple[np.ndarray, np.ndarray]:
    x, w = np.polynomial.legendre.leggauss(order)
    return 0.5 * ((hi - lo) * x + hi + lo), 0.5 * (hi - lo) * w


class FiveAngleProjector:
    def __init__(self, constants: _OP.N2LOConstants | None = None) -> None:
        self.constants = constants or _OP.N2LOConstants.tictac()

    def project(
        self,
        bra: LSChannel,
        ket: LSChannel,
        momenta: Tuple[float, float, float, float],
        lecs: _OP.N2LOLECs,
        order: int,
        two_m_t: int = 1,
    ) -> Dict[str, complex]:
        return self._project(
            bra, ket, momenta, lecs, order, two_m_t, angular_spin_state
        )

    def project_jj_direct(
        self,
        bra: JjChannel,
        ket: JjChannel,
        momenta: Tuple[float, float, float, float],
        lecs: _OP.N2LOLECs,
        order: int,
        two_m_t: int = 1,
    ) -> Dict[str, complex]:
        """Project directly in the Tic-tac Jj basis using explicit m sums."""
        return self._project(
            bra, ket, momenta, lecs, order, two_m_t, angular_spin_state_jj
        )

    def project_jj_recoupled(
        self,
        bra: JjChannel,
        ket: JjChannel,
        momenta: Tuple[float, float, float, float],
        lecs: _OP.N2LOLECs,
        order: int,
        two_m_t: int = 1,
    ) -> Dict[str, complex]:
        """Project in LS channels and transform with the unitary 9j map."""
        totals = {name: 0.0j for name in ("c1", "c3", "c4", "cD", "cE")}
        for ls_bra, bra_coefficient in bra.ls_expansion():
            for ls_ket, ket_coefficient in ket.ls_expansion():
                block = self.project(ls_bra, ls_ket, momenta, lecs, order, two_m_t)
                factor = bra_coefficient * ket_coefficient
                for name, value in block.items():
                    totals[name] += factor * value
        return totals

    def _project(
        self,
        bra: LSChannel | JjChannel,
        ket: LSChannel | JjChannel,
        momenta: Tuple[float, float, float, float],
        lecs: _OP.N2LOLECs,
        order: int,
        two_m_t: int,
        state_builder,
    ) -> Dict[str, complex]:
        """Project all five components; momenta are (p,q,p',q') in fm^-1."""
        if bra.two_total_J != ket.two_total_J:
            return {name: 0.0j for name in ("c1", "c3", "c4", "cD", "cE")}
        if bra.two_total_T != ket.two_total_T:
            return {name: 0.0j for name in ("c1", "c3", "c4", "cD", "cE")}
        if order < 1:
            raise ValueError("order must be positive")

        p, q, pp, qp = map(float, momenta)
        cosines, cosine_weights = _nodes_weights(order, -1.0, 1.0)
        phis, phi_weights = _nodes_weights(order, 0.0, 2.0 * math.pi)
        sines = np.sqrt(np.maximum(0.0, 1.0 - cosines * cosines))

        p_in = np.array([0.0, 0.0, p])
        iso_bra = _isospin_state(bra, two_m_t)
        iso_ket = _isospin_state(ket, two_m_t)
        iso23 = np.vdot(iso_bra, _OP._TAU23 @ iso_ket)
        iso_cross = np.vdot(iso_bra, _OP._TAU1_DOT_TAU2_CROSS_TAU3 @ iso_ket)
        iso12 = np.vdot(iso_bra, _OP._TAU12 @ iso_ket)
        iso13 = np.vdot(iso_bra, _OP._TAU13 @ iso_ket)

        c = self.constants
        f_pi = c.f_pi
        m_pi = c.m_pi
        lambda_chi = c.lambda_chi
        c1 = c.gev_inverse_to_fm(lecs.c1_gev_inverse)
        c3 = c.gev_inverse_to_fm(lecs.c3_gev_inverse)
        c4 = c.gev_inverse_to_fm(lecs.c4_gev_inverse)
        d_lec = lecs.c_d / (f_pi * f_pi * lambda_chi)
        e_lec = lecs.c_e / (f_pi**4 * lambda_chi)

        totals = {name: 0.0j for name in ("c1", "c3", "c4", "cD", "cE")}
        two_j = bra.two_total_J
        m_j_values = tuple(_m_values(two_j))

        for cq, sq, wq in zip(cosines, sines, cosine_weights):
            q_in = q * np.array([sq, 0.0, cq])
            ket_spin = {
                two_m_j: state_builder(ket, p_in, q_in, two_m_j)
                for two_m_j in m_j_values
            }
            for cpp, spp, wpp in zip(cosines, sines, cosine_weights):
                for phip, wphip in zip(phis, phi_weights):
                    p_out = pp * np.array([spp * math.cos(phip), spp * math.sin(phip), cpp])
                    for cqp, sqp, wqp in zip(cosines, sines, cosine_weights):
                        for phiqp, wphiqp in zip(phis, phi_weights):
                            q_out = qp * np.array(
                                [sqp * math.cos(phiqp), sqp * math.sin(phiqp), cqp]
                            )
                            bra_spin = {
                                two_m_j: state_builder(bra, p_out, q_out, two_m_j)
                                for two_m_j in m_j_values
                            }
                            q1, q2, q3 = _OP.SpectatorOneN2LOOracle.transfers(
                                p_out, q_out, p_in, q_in
                            )
                            d2 = float(np.dot(q2, q2) + m_pi * m_pi)
                            d3 = float(np.dot(q3, q3) + m_pi * m_pi)
                            d1 = float(np.dot(q1, q1) + m_pi * m_pi)
                            common = c.g_a * c.g_a / (4.0 * f_pi**4 * d2 * d3)
                            q2q3 = float(np.dot(q2, q3))

                            spin23 = _OP._pauli_dot(2, q2) @ _OP._pauli_dot(3, q3)
                            spin4 = spin23 @ _OP._pauli_dot(1, np.cross(q2, q3))
                            spin_d2 = _OP._pauli_dot(1, q1) @ _OP._pauli_dot(2, q1)
                            spin_d3 = _OP._pauli_dot(1, q1) @ _OP._pauli_dot(3, q1)
                            spin_matrices = {
                                "c1": spin23,
                                "c3": spin23,
                                "c4": spin4,
                                "cD2": spin_d2,
                                "cD3": spin_d3,
                                "cE": np.eye(8, dtype=complex),
                            }
                            spin_me = {}
                            for name, matrix in spin_matrices.items():
                                spin_me[name] = sum(
                                    np.vdot(bra_spin[mj], matrix @ ket_spin[mj])
                                    for mj in m_j_values
                                ) / (two_j + 1.0)

                            weight = wq * wpp * wphip * wqp * wphiqp
                            totals["c1"] += weight * common * (-4.0 * c1 * m_pi**2) * iso23 * spin_me["c1"]
                            totals["c3"] += weight * common * (2.0 * c3 * q2q3) * iso23 * spin_me["c3"]
                            totals["c4"] += weight * common * c4 * iso_cross * spin_me["c4"]
                            totals["cD"] += weight * (
                                -c.g_a * d_lec / (8.0 * f_pi * f_pi * d1)
                            ) * (iso12 * spin_me["cD2"] + iso13 * spin_me["cD3"])
                            totals["cE"] += weight * e_lec * iso23 * spin_me["cE"]

        rotational_volume = 8.0 * math.pi**2
        return {name: rotational_volume * value for name, value in totals.items()}

    @staticmethod
    def to_tictac_normalization(raw_components: Dict[str, complex]) -> Dict[str, complex]:
        """Convert raw Golak/Epelbaum PWD values to Tic-tac's Jacobi basis."""
        return {
            name: THREE_BODY_FOURIER_NORMALIZATION * value
            for name, value in raw_components.items()
        }


GOLAK_BETA = {
    1: LSChannel(0, 0, 0, 0, 1, 1, 1),
    2: LSChannel(0, 1, 0, 0, 1, 1, 0),
    3: LSChannel(0, 1, 2, 2, 3, 1, 0),
    4: LSChannel(1, 0, 1, 0, 1, 1, 0),
    5: LSChannel(1, 0, 1, 1, 1, 1, 0),
    6: LSChannel(1, 1, 1, 0, 1, 1, 1),
    7: LSChannel(1, 1, 1, 1, 1, 1, 1),
    8: LSChannel(1, 1, 1, 1, 3, 1, 1),
    9: LSChannel(1, 1, 1, 2, 3, 1, 1),
    10: LSChannel(2, 1, 0, 2, 3, 1, 0),
    11: LSChannel(2, 0, 2, 0, 1, 1, 1),
    12: LSChannel(2, 0, 2, 1, 1, 1, 1),
    13: LSChannel(2, 1, 2, 0, 1, 1, 0),
    14: LSChannel(2, 1, 2, 1, 1, 1, 0),
    15: LSChannel(2, 1, 2, 1, 3, 1, 0),
    16: LSChannel(2, 1, 2, 2, 3, 1, 0),
}


if __name__ == "__main__":
    projector = FiveAngleProjector(_OP.N2LOConstants(197.327, 1.29, 92.4, 138.0))
    values = projector.project(
        GOLAK_BETA[1],
        GOLAK_BETA[1],
        (1.0, 2.0, 3.0, 4.0),
        _OP.N2LOLECs(c1_gev_inverse=-0.81, c3_gev_inverse=-3.4),
        order=6,
    )
    print(f"generic five-angle G(1,1), N=6: {(values['c1'] + values['c3']).real:.9f} fm^5")
