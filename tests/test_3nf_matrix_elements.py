#!/usr/bin/env python3
"""L2: Independent Python computation of W^(1) matrix elements as oracle."""

import math
import unittest
from sympy import S, sqrt, Rational
from sympy.physics.wigner import clebsch_gordan, wigner_6j, wigner_9j


# Physical constants (must match include/constants.h)
GAxial = 1.289
FPI_MEV = 92.2
MPI_MEV = 138.039
HBARC = 197.327
LAMBDA_CHI_MEV = 700.0

FPI_FM = FPI_MEV / HBARC
MPI_FM = MPI_MEV / HBARC
LAMBDA_CHI_FM = LAMBDA_CHI_MEV / HBARC
FPI4_INV = FPI_FM ** 4


def gaussian_regulator(p, q, Lambda_fm):
    return math.exp(-(p**2 + 0.75 * q**2) / Lambda_fm**2)


def w1_contact(T_2N, c_E, p_r, q_r, p_c, q_c, Lambda_fm):
    """Contact term W^(1)_CT for diagonal alpha."""
    tau23 = 2.0 * T_2N * (T_2N + 1) - 3.0
    coeff = c_E / (FPI4_INV * LAMBDA_CHI_FM)
    f_bra = gaussian_regulator(p_r, q_r, Lambda_fm)
    f_ket = gaussian_regulator(p_c, q_c, Lambda_fm)
    return coeff * tau23 * f_bra * f_ket


def w1_2pe_scalar(S_2N, T_2N, c1, c3, p_r, q_r, p_c, q_c, Lambda_fm):
    """2PE scalar (rank-0) term W^(1)_2PE for diagonal alpha.

    Implements the monopole/scalar approximation:
      (sigma2.q2)(sigma3.q3) -> (1/3)(sigma2.sigma3)(q2.q3)
    with angle-averaged pion propagators.
    """
    if c1 == 0.0 and c3 == 0.0:
        return 0.0

    sig23 = 2.0 * S_2N * (S_2N + 1) - 3.0
    tau23 = 2.0 * T_2N * (T_2N + 1) - 3.0

    A_sq = p_r**2 + p_c**2
    B_sq = 0.25 * (q_r**2 + q_c**2)
    Q_sq = A_sq + B_sq
    q2_dot_q3 = B_sq - A_sq

    denom = Q_sq + MPI_FM**2
    prop = 1.0 / (denom**2)

    c1_part = -4.0 * c1 * MPI_FM**2 * q2_dot_q3
    c3_part = 2.0 * c3 * q2_dot_q3**2

    scalar = (1.0 / 3.0) * sig23 * tau23 * (c1_part + c3_part) * prop

    coeff = GAxial**2 / (4.0 * FPI_FM**4)

    f_bra = gaussian_regulator(p_r, q_r, Lambda_fm)
    f_ket = gaussian_regulator(p_c, q_c, Lambda_fm)

    return coeff * scalar * f_bra * f_ket


class TestTwoPionExchangeTerm(unittest.TestCase):
    def test_c1_sign_T0S0(self):
        """For T_2N=0, S_2N=0: sig23=-3, tau23=-3, sig*tau=+9.
        With c1<0 and p>q/2 (so q2.q3<0):
          c1_part = -4*c1*mpi^2*q2.q3 = -4*(-0.81)*mpi^2*(-3.28) < 0
        scalar = (1/3)*9*c1_part*prop < 0. coeff>0, so result < 0."""
        val = w1_2pe_scalar(0, 0, -0.81, 0.0, 1.5, 0.8, 1.2, 1.0, 500.0 / HBARC)
        self.assertLess(val, 0.0)

    def test_c1_sign_T1S1(self):
        """For T_2N=1, S_2N=1: sig23=+1, tau23=+1, sig*tau=+1.
        Same momentum config: c1_part<0, so result < 0."""
        val = w1_2pe_scalar(1, 1, -0.81, 0.0, 1.5, 0.8, 1.2, 1.0, 500.0 / HBARC)
        self.assertLess(val, 0.0)

    def test_zero_when_c1c3_zero(self):
        val = w1_2pe_scalar(1, 0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 500.0 / HBARC)
        self.assertEqual(val, 0.0)

    def test_regulator_suppresses_high_momentum(self):
        val_low = abs(w1_2pe_scalar(1, 1, -0.81, -3.2, 0.5, 0.5, 0.5, 0.5, 500.0 / HBARC))
        val_high = abs(w1_2pe_scalar(1, 1, -0.81, -3.2, 5.0, 5.0, 5.0, 5.0, 500.0 / HBARC))
        self.assertGreater(val_low / (val_high + 1e-300), 100.0)

    def test_c3_dominates_at_large_momentum_transfer(self):
        """c3 term goes as (q2.q3)^2 while c1 goes as (q2.q3), so c3 dominates
        at large momentum transfers."""
        # Use momenta where A_sq >> B_sq so |q2.q3| is large
        p, q = 3.0, 0.5
        val_c1_only = abs(w1_2pe_scalar(1, 1, -0.81, 0.0, p, q, p, q, 500.0 / HBARC))
        val_c3_only = abs(w1_2pe_scalar(1, 1, 0.0, -3.2, p, q, p, q, 500.0 / HBARC))
        self.assertGreater(val_c3_only, val_c1_only)


class TestContactTerm(unittest.TestCase):
    def test_T0_negative(self):
        """T_2N=0 gives tau23=-3, so contact term is negative for positive c_E."""
        val = w1_contact(0, 1.0, 1.0, 1.0, 1.0, 1.0, 500.0 / HBARC)
        self.assertLess(val, 0.0)

    def test_T1_positive(self):
        """T_2N=1 gives tau23=+1, so contact term is positive for positive c_E."""
        val = w1_contact(1, 1.0, 1.0, 1.0, 1.0, 1.0, 500.0 / HBARC)
        self.assertGreater(val, 0.0)

    def test_regulator_suppresses_high_momentum(self):
        """High momenta should be strongly suppressed."""
        val_low = abs(w1_contact(1, 1.0, 0.5, 0.5, 0.5, 0.5, 500.0 / HBARC))
        val_high = abs(w1_contact(1, 1.0, 5.0, 5.0, 5.0, 5.0, 500.0 / HBARC))
        self.assertGreater(val_low / (val_high + 1e-300), 100.0)

    def test_zero_cE(self):
        self.assertEqual(w1_contact(1, 0.0, 1.0, 1.0, 1.0, 1.0, 500.0 / HBARC), 0.0)


if __name__ == "__main__":
    unittest.main()
