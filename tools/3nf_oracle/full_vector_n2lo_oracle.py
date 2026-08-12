#!/usr/bin/env python3
"""Independent full-vector m-scheme oracle for the local chiral N2LO 3NF.

This module evaluates the spectator-1 component ``u1 = V^(1)`` before any
partial-wave projection.  It deliberately uses explicit 8x8 Pauli matrices in
the three-particle spin and isospin product spaces; no production recoupling or
partial-wave routine is imported.

Conventions
-----------
For equal masses in the three-body centre-of-mass frame,

    k1 = q,  k2 = p - q/2,  k3 = -p - q/2,

and Q_i = k_i' - k_i.  Thus Q1+Q2+Q3=0.  The formulas are the spectator-1
pieces obtained from Epelbaum et al., Phys. Rev. C 66, 064001 (2002),
Eqs. (2.2), (2.10), and (2.12): the c1/c3/c4 piece has the subleading pi-pi-N
vertex on particle 1, while the cD pion is carried by particle 1.  The cE
component is E tau2.tau3.  The factor 1/2 in the full ordered-pair sum is
cancelled by its two (2,3)/(3,2) occurrences.

All momenta and mass scales passed to the algebra are in fm^-1.  c1/c3/c4 are
accepted in GeV^-1 and converted to fm.  The returned momentum-space operator
has units fm^5 *before* any plane-wave/Fourier state-normalisation factor.  The
caller may supply such a factor explicitly; its default is one so that a
Fourier convention cannot be smuggled into the physics oracle.

The 64-dimensional basis is ``|spin1 spin2 spin3> tensor
|isospin1 isospin2 isospin3>``, with each one-particle up state preceding down.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterable, Mapping, Tuple

import numpy as np


Vector = np.ndarray
Matrix = np.ndarray


_I2 = np.eye(2, dtype=complex)
_PAULI = (
    np.array([[0.0, 1.0], [1.0, 0.0]], dtype=complex),
    np.array([[0.0, -1.0j], [1.0j, 0.0]], dtype=complex),
    np.array([[1.0, 0.0], [0.0, -1.0]], dtype=complex),
)
_I8 = np.eye(8, dtype=complex)


def _embed_one_body(particle: int, one_body: Matrix) -> Matrix:
    """Embed a 2x2 matrix on particle 1, 2, or 3 in an 8-state space."""
    if particle not in (1, 2, 3):
        raise ValueError(f"particle must be 1, 2, or 3, got {particle}")
    factors = [_I2, _I2, _I2]
    factors[particle - 1] = one_body
    return np.kron(np.kron(factors[0], factors[1]), factors[2])


def _three_particle_pauli() -> Tuple[Tuple[Matrix, Matrix, Matrix], ...]:
    return tuple(
        tuple(_embed_one_body(particle, _PAULI[axis]) for axis in range(3))
        for particle in (1, 2, 3)
    )


_SIGMA = _three_particle_pauli()
_TAU = _three_particle_pauli()


def _as_vector(value: Iterable[float], name: str) -> Vector:
    out = np.asarray(tuple(value), dtype=float)
    if out.shape != (3,):
        raise ValueError(f"{name} must contain three Cartesian components")
    if not np.all(np.isfinite(out)):
        raise ValueError(f"{name} contains a non-finite component")
    return out


def _pauli_dot(particle: int, vector: Vector, operators=_SIGMA) -> Matrix:
    return sum(vector[a] * operators[particle - 1][a] for a in range(3))


def _pair_dot(first: int, second: int, operators=_TAU) -> Matrix:
    return sum(
        operators[first - 1][a] @ operators[second - 1][a]
        for a in range(3)
    )


def _triple_cross_123(operators=_TAU) -> Matrix:
    """Return tau1.(tau2 x tau3), with epsilon_xyz=+1."""
    result = np.zeros((8, 8), dtype=complex)
    for a, b, c, sign in (
        (0, 1, 2, +1), (1, 2, 0, +1), (2, 0, 1, +1),
        (1, 0, 2, -1), (2, 1, 0, -1), (0, 2, 1, -1),
    ):
        result += sign * (
            operators[0][a] @ operators[1][b] @ operators[2][c]
        )
    return result


_TAU12 = _pair_dot(1, 2)
_TAU13 = _pair_dot(1, 3)
_TAU23 = _pair_dot(2, 3)
_TAU1_DOT_TAU2_CROSS_TAU3 = _triple_cross_123()


@dataclass(frozen=True)
class N2LOConstants:
    """Dimensionful constants kept explicit rather than imported from production."""

    hbar_c_mev_fm: float
    g_a: float
    f_pi_mev: float
    m_pi_mev: float
    lambda_chi_mev: float = 700.0

    @classmethod
    def tictac(cls) -> "N2LOConstants":
        return cls(197.327, 1.289, 92.2, 138.039, 700.0)

    @classmethod
    def epelbaum_2002(cls) -> "N2LOConstants":
        # The paper explicitly states g_A=1.276, f_pi=92.4 MeV and
        # Lambda_chi=700 MeV.  It denotes the pion mass by M_pi without fixing
        # one numerical charge average in Eq. (2.2), so 138.039 MeV is recorded
        # here as a transparent external choice.
        return cls(197.327, 1.276, 92.4, 138.039, 700.0)

    @property
    def f_pi(self) -> float:
        return self.f_pi_mev / self.hbar_c_mev_fm

    @property
    def m_pi(self) -> float:
        return self.m_pi_mev / self.hbar_c_mev_fm

    @property
    def lambda_chi(self) -> float:
        return self.lambda_chi_mev / self.hbar_c_mev_fm

    def gev_inverse_to_fm(self, value: float) -> float:
        return value * self.hbar_c_mev_fm / 1000.0


@dataclass(frozen=True)
class N2LOLECs:
    c1_gev_inverse: float = 0.0
    c3_gev_inverse: float = 0.0
    c4_gev_inverse: float = 0.0
    c_d: float = 0.0
    c_e: float = 0.0


class SpectatorOneN2LOOracle:
    """Evaluate full Cartesian spectator-1 N2LO operators in the m scheme."""

    def __init__(
        self,
        constants: N2LOConstants | None = None,
        cutoff_mev: float | None = None,
        state_normalization: float = 1.0,
    ) -> None:
        self.constants = constants or N2LOConstants.tictac()
        self.cutoff = (
            None
            if cutoff_mev is None
            else float(cutoff_mev) / self.constants.hbar_c_mev_fm
        )
        self.state_normalization = float(state_normalization)
        if self.cutoff is not None and self.cutoff <= 0.0:
            raise ValueError("cutoff_mev must be positive or None")
        if not np.isfinite(self.state_normalization):
            raise ValueError("state_normalization must be finite")

    @staticmethod
    def single_particle_momenta(p: Iterable[float], q: Iterable[float]) -> Tuple[Vector, ...]:
        p_vec = _as_vector(p, "p")
        q_vec = _as_vector(q, "q")
        return q_vec, p_vec - 0.5 * q_vec, -p_vec - 0.5 * q_vec

    @classmethod
    def transfers(
        cls,
        p_bra: Iterable[float],
        q_bra: Iterable[float],
        p_ket: Iterable[float],
        q_ket: Iterable[float],
    ) -> Tuple[Vector, Vector, Vector]:
        final = cls.single_particle_momenta(p_bra, q_bra)
        initial = cls.single_particle_momenta(p_ket, q_ket)
        return tuple(final[i] - initial[i] for i in range(3))

    def regulator(self, p: Vector, q: Vector) -> float:
        if self.cutoff is None:
            return 1.0
        invariant = float(np.dot(p, p) + 0.75 * np.dot(q, q))
        return float(np.exp(-((invariant / (self.cutoff * self.cutoff)) ** 2)))

    def component_operators(
        self,
        p_bra: Iterable[float],
        q_bra: Iterable[float],
        p_ket: Iterable[float],
        q_ket: Iterable[float],
        lecs: N2LOLECs,
    ) -> Dict[str, Matrix]:
        """Return separate c1, c3, c4, cD, and cE 64x64 matrices in fm^5."""
        p_out = _as_vector(p_bra, "p_bra")
        q_out = _as_vector(q_bra, "q_bra")
        p_in = _as_vector(p_ket, "p_ket")
        q_in = _as_vector(q_ket, "q_ket")
        q1, q2, q3 = self.transfers(p_out, q_out, p_in, q_in)
        if not np.allclose(q1 + q2 + q3, 0.0, atol=2e-14, rtol=0.0):
            raise ArithmeticError("Jacobi transfer construction violated momentum conservation")

        c = self.constants
        f_pi = c.f_pi
        m_pi = c.m_pi
        lambda_chi = c.lambda_chi
        c1 = c.gev_inverse_to_fm(lecs.c1_gev_inverse)
        c3 = c.gev_inverse_to_fm(lecs.c3_gev_inverse)
        c4 = c.gev_inverse_to_fm(lecs.c4_gev_inverse)

        regulator = (
            self.regulator(p_out, q_out)
            * self.regulator(p_in, q_in)
            * self.state_normalization
        )
        zero = np.zeros((64, 64), dtype=complex)

        d2 = float(np.dot(q2, q2) + m_pi * m_pi)
        d3 = float(np.dot(q3, q3) + m_pi * m_pi)
        spin_2pi = _pauli_dot(2, q2) @ _pauli_dot(3, q3)
        tpe_common = c.g_a * c.g_a / (4.0 * f_pi**4 * d2 * d3)

        components: Dict[str, Matrix] = {}
        components["c1"] = (
            regulator
            * tpe_common
            * (-4.0 * c1 * m_pi * m_pi)
            * np.kron(spin_2pi, _TAU23)
            if c1 != 0.0 else zero.copy()
        )
        components["c3"] = (
            regulator
            * tpe_common
            * (2.0 * c3 * float(np.dot(q2, q3)))
            * np.kron(spin_2pi, _TAU23)
            if c3 != 0.0 else zero.copy()
        )

        cross_23 = np.cross(q2, q3)
        spin_c4 = spin_2pi @ _pauli_dot(1, cross_23)
        components["c4"] = (
            regulator
            * tpe_common
            * c4
            * np.kron(spin_c4, _TAU1_DOT_TAU2_CROSS_TAU3)
            if c4 != 0.0 else zero.copy()
        )

        d1 = float(np.dot(q1, q1) + m_pi * m_pi)
        d_lec = lecs.c_d / (f_pi * f_pi * lambda_chi)
        spin_1 = _pauli_dot(1, q1)
        spin_d2 = spin_1 @ _pauli_dot(2, q1)
        spin_d3 = spin_1 @ _pauli_dot(3, q1)
        d_operator = np.kron(spin_d2, _TAU12) + np.kron(spin_d3, _TAU13)
        components["cD"] = (
            regulator
            * (-c.g_a * d_lec / (8.0 * f_pi * f_pi * d1))
            * d_operator
            if lecs.c_d != 0.0 else zero.copy()
        )

        e_lec = lecs.c_e / (f_pi**4 * lambda_chi)
        components["cE"] = (
            regulator * e_lec * np.kron(_I8, _TAU23)
            if lecs.c_e != 0.0 else zero.copy()
        )
        return components

    def operator(self, *args, **kwargs) -> Matrix:
        return sum(self.component_operators(*args, **kwargs).values())

    @staticmethod
    def product_basis_state(bits: Tuple[int, int, int]) -> Vector:
        """Return |b1 b2 b3>, where 0 is up and 1 is down."""
        if len(bits) != 3 or any(bit not in (0, 1) for bit in bits):
            raise ValueError("bits must be a three-entry tuple containing only 0 or 1")
        out = np.zeros(8, dtype=complex)
        out[4 * bits[0] + 2 * bits[1] + bits[2]] = 1.0
        return out

    @staticmethod
    def matrix_element(
        operator: Matrix,
        spin_bra: Vector,
        isospin_bra: Vector,
        spin_ket: Vector,
        isospin_ket: Vector,
    ) -> complex:
        bra = np.kron(np.asarray(spin_bra), np.asarray(isospin_bra))
        ket = np.kron(np.asarray(spin_ket), np.asarray(isospin_ket))
        return complex(np.vdot(bra, operator @ ket))


def _self_check() -> None:
    oracle = SpectatorOneN2LOOracle(cutoff_mev=500.0)
    lecs = N2LOLECs(-0.81, -3.2, 5.4, -0.2, -0.205)
    p = np.array([0.31, -0.22, 0.48])
    q = np.array([-0.19, 0.41, 0.27])
    pp = np.array([-0.37, 0.16, 0.52])
    qp = np.array([0.24, -0.33, 0.29])
    terms = oracle.component_operators(pp, qp, p, q, lecs)
    reverse = oracle.component_operators(p, q, pp, qp, lecs)
    for name in terms:
        np.testing.assert_allclose(terms[name], reverse[name].conj().T, atol=2e-12)
    q_transfer = oracle.transfers(pp, qp, p, q)
    np.testing.assert_allclose(sum(q_transfer), np.zeros(3), atol=2e-14)
    print("full-vector N2LO spectator-1 oracle self-check: PASS")
    for name, matrix in terms.items():
        print(f"  {name:>2}: ||u1||_F = {np.linalg.norm(matrix):.12e} fm^5")


if __name__ == "__main__":
    _self_check()
