"""Validation of the finite-rank Hebeler N2LO factorization."""

from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest


REPO = Path(__file__).resolve().parents[1]
ORACLE_DIR = REPO / "tools" / "3nf_oracle"
sys.path.insert(0, str(ORACLE_DIR))

import factorized_n2lo_pwd as FACTORIZED  # noqa: E402
import full_vector_five_angle_pwd as PWD  # noqa: E402


class FactorizedN2LOPWDTests(unittest.TestCase):
    def setUp(self):
        self.constants = PWD._OP.N2LOConstants.tictac()
        self.momenta = (0.43, 0.61, 0.78, 0.69)
        self.direct = PWD.FiveAngleProjector(self.constants)

    def test_factorized_two_pion_terms_reproduce_golak_table2(self):
        constants = PWD._OP.N2LOConstants(197.327, 1.29, 92.4, 138.0)
        momenta = (1.0, 2.0, 3.0, 4.0)
        c13 = FACTORIZED.project_c1_c3_ls(
            PWD.GOLAK_BETA[1],
            PWD.GOLAK_BETA[1],
            momenta,
            constants,
            c1_gev_inverse=-0.81,
            c3_gev_inverse=-3.4,
            order=12,
        )
        c4 = FACTORIZED.project_c4_ls(
            PWD.GOLAK_BETA[2],
            PWD.GOLAK_BETA[1],
            momenta,
            constants,
            c4_gev_inverse=3.4,
            order=12,
        )
        values = {
            "G(1,1)": (c13["c1"] + c13["c3"]).real,
            "G(2,1)": c4.real,
        }
        published = {"G(1,1)": 443.618, "G(2,1)": 1200.219}
        for name, target in published.items():
            with self.subTest(matrix_element=name):
                self.assertLess(abs(values[name] - target) / abs(target), 3.0e-4)
        self.assertLess(abs(c13["c1"].imag + c13["c3"].imag), 1.0e-11)
        self.assertLess(abs(c4.imag), 1.0e-11)

    def test_complete_cD_matches_five_angle_s_and_spectator_d_transitions(self):
        cases = (
            (PWD.GOLAK_BETA[2], PWD.GOLAK_BETA[1]),
            (PWD.GOLAK_BETA[1], PWD.GOLAK_BETA[3]),
            (PWD.GOLAK_BETA[3], PWD.GOLAK_BETA[1]),
        )
        lecs = PWD._OP.N2LOLECs(c_d=-0.2)
        for bra, ket in cases:
            with self.subTest(bra=bra, ket=ket):
                factorized = FACTORIZED.project_cD_ls(
                    bra,
                    ket,
                    self.momenta,
                    self.constants,
                    c_d=-0.2,
                    order=8,
                )
                direct = self.direct.project(
                    bra, ket, self.momenta, lecs, order=8
                )["cD"]
                tolerance = 1.0e-4 * max(1.0, abs(factorized.real))
                self.assertAlmostEqual(
                    factorized.real, direct.real, delta=tolerance
                )
                self.assertLess(abs(factorized.imag), 1.0e-10)
                self.assertLess(abs(direct.imag), 1.0e-10)

    def test_factorized_contact_matches_five_angle_normalization(self):
        bra = PWD.GOLAK_BETA[1]
        ket = PWD.GOLAK_BETA[1]
        factorized = FACTORIZED.project_cE_ls(
            bra,
            ket,
            self.momenta,
            self.constants,
            c_e=-0.205,
            order=6,
        )
        direct = self.direct.project(
            bra,
            ket,
            self.momenta,
            PWD._OP.N2LOLECs(c_e=-0.205),
            order=2,
        )["cE"]
        self.assertAlmostEqual(factorized.real, direct.real, places=10)
        self.assertLess(abs(factorized.imag), 1.0e-11)
        expected_scalar_factor = (4.0 * math.pi) ** 2
        e_lec = -0.205 / (self.constants.f_pi**4 * self.constants.lambda_chi)
        self.assertAlmostEqual(
            factorized.real, expected_scalar_factor * e_lec, places=10
        )


if __name__ == "__main__":
    unittest.main()

