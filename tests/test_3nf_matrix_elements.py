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
