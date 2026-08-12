"""Structural tests for the independent full-vector N2LO 3NF oracle."""

import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np


REPO = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO / "tools" / "3nf_oracle" / "full_vector_n2lo_oracle.py"
SPEC = importlib.util.spec_from_file_location("full_vector_n2lo_oracle", MODULE_PATH)
ORACLE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = ORACLE
SPEC.loader.exec_module(ORACLE)

GOLAK_PATH = REPO / "tools" / "3nf_oracle" / "golak_table2_benchmark.py"
GOLAK_SPEC = importlib.util.spec_from_file_location("golak_table2_benchmark", GOLAK_PATH)
GOLAK = importlib.util.module_from_spec(GOLAK_SPEC)
assert GOLAK_SPEC.loader is not None
sys.modules[GOLAK_SPEC.name] = GOLAK
GOLAK_SPEC.loader.exec_module(GOLAK)

PWD_PATH = REPO / "tools" / "3nf_oracle" / "full_vector_five_angle_pwd.py"
PWD_SPEC = importlib.util.spec_from_file_location("full_vector_five_angle_pwd", PWD_PATH)
PWD = importlib.util.module_from_spec(PWD_SPEC)
assert PWD_SPEC.loader is not None
sys.modules[PWD_SPEC.name] = PWD
PWD_SPEC.loader.exec_module(PWD)


class FullVectorN2LOOracleTests(unittest.TestCase):
    def setUp(self):
        self.oracle = ORACLE.SpectatorOneN2LOOracle(cutoff_mev=None)
        self.p = np.array([0.31, -0.22, 0.48])
        self.q = np.array([-0.19, 0.41, 0.27])
        self.pp = np.array([-0.37, 0.16, 0.52])
        self.qp = np.array([0.24, -0.33, 0.29])

    def test_jacobi_transfers_conserve_momentum(self):
        transfers = self.oracle.transfers(self.pp, self.qp, self.p, self.q)
        np.testing.assert_allclose(sum(transfers), np.zeros(3), atol=2e-14, rtol=0.0)

    def test_each_lec_is_independently_switchable_and_nonzero(self):
        cases = {
            "c1": ORACLE.N2LOLECs(c1_gev_inverse=1.0),
            "c3": ORACLE.N2LOLECs(c3_gev_inverse=1.0),
            "c4": ORACLE.N2LOLECs(c4_gev_inverse=1.0),
            "cD": ORACLE.N2LOLECs(c_d=1.0),
            "cE": ORACLE.N2LOLECs(c_e=1.0),
        }
        for selected, lecs in cases.items():
            with self.subTest(selected=selected):
                terms = self.oracle.component_operators(self.pp, self.qp, self.p, self.q, lecs)
                self.assertGreater(np.linalg.norm(terms[selected]), 1e-12)
                for name, value in terms.items():
                    if name != selected:
                        self.assertEqual(np.count_nonzero(value), 0)

    def test_cE_ordered_pair_counting_and_pair_isospin_eigenvalues(self):
        terms = self.oracle.component_operators(
            self.pp, self.qp, self.p, self.q, ORACLE.N2LOLECs(c_e=1.0)
        )
        op = terms["cE"]
        up = self.oracle.product_basis_state((0, 0, 0))
        pair_singlet = (
            self.oracle.product_basis_state((0, 0, 1))
            - self.oracle.product_basis_state((0, 1, 0))
        ) / np.sqrt(2.0)
        spin = up
        e_lec = 1.0 / (
            self.oracle.constants.f_pi**4 * self.oracle.constants.lambda_chi
        )
        triplet_value = self.oracle.matrix_element(op, spin, up, spin, up)
        singlet_value = self.oracle.matrix_element(
            op, spin, pair_singlet, spin, pair_singlet
        )
        self.assertAlmostEqual(triplet_value.real, +e_lec, places=13)
        self.assertAlmostEqual(singlet_value.real, -3.0 * e_lec, places=13)
        self.assertAlmostEqual(triplet_value.imag, 0.0, places=14)
        self.assertAlmostEqual(singlet_value.imag, 0.0, places=14)

    def test_reverse_kernel_is_hermitian_for_every_component(self):
        lecs = ORACLE.N2LOLECs(-0.81, -3.2, 5.4, -0.2, -0.205)
        forward = self.oracle.component_operators(self.pp, self.qp, self.p, self.q, lecs)
        reverse = self.oracle.component_operators(self.p, self.q, self.pp, self.qp, lecs)
        for name in forward:
            with self.subTest(component=name):
                np.testing.assert_allclose(
                    forward[name], reverse[name].conj().T, atol=2e-12, rtol=2e-13
                )

    @staticmethod
    def _permutation_23() -> np.ndarray:
        result = np.zeros((8, 8), dtype=complex)
        for b1 in (0, 1):
            for b2 in (0, 1):
                for b3 in (0, 1):
                    old = 4 * b1 + 2 * b2 + b3
                    new = 4 * b1 + 2 * b3 + b2
                    result[new, old] = 1.0
        return result

    def test_spectator_component_is_symmetric_under_pair_exchange(self):
        lecs = ORACLE.N2LOLECs(-0.81, -3.2, 5.4, -0.2, -0.205)
        original = self.oracle.component_operators(self.pp, self.qp, self.p, self.q, lecs)
        exchanged = self.oracle.component_operators(-self.pp, self.qp, -self.p, self.q, lecs)
        p23 = self._permutation_23()
        full_p23 = np.kron(p23, p23)
        for name in original:
            with self.subTest(component=name):
                np.testing.assert_allclose(
                    exchanged[name],
                    full_p23 @ original[name] @ full_p23,
                    atol=2e-12,
                    rtol=2e-13,
                )

    def test_c4_cross_product_zero_and_nonzero_discriminator(self):
        lecs = ORACLE.N2LOLECs(c4_gev_inverse=1.0)
        generic = self.oracle.component_operators(self.pp, self.qp, self.p, self.q, lecs)["c4"]
        self.assertGreater(np.linalg.norm(generic), 1e-10)

        # Pure spectator transfer with Delta p=0 gives Q2=Q3=-Q1/2, so
        # Q2 x Q3 vanishes exactly and isolates the c4 cross product.
        collinear = self.oracle.component_operators(
            self.p, self.qp, self.p, self.q, lecs
        )["c4"]
        self.assertLess(np.linalg.norm(collinear), 1e-14)
        self.assertGreater(np.max(np.abs(generic.imag)), 1e-10)

    def test_c4_isospin_recoupling_has_the_full_two_sqrt_three_factor(self):
        basis = self.oracle.product_basis_state
        t23_zero = (basis((0, 0, 1)) - basis((0, 1, 0))) / np.sqrt(2.0)
        t23_one = (
            np.sqrt(2.0 / 3.0) * basis((1, 0, 0))
            - np.sqrt(1.0 / 6.0) * (basis((0, 0, 1)) + basis((0, 1, 0)))
        )
        value = np.vdot(
            t23_zero,
            ORACLE._TAU1_DOT_TAU2_CROSS_TAU3 @ t23_one,
        )
        self.assertAlmostEqual(value.real, 0.0, places=14)
        self.assertAlmostEqual(value.imag, -2.0 * np.sqrt(3.0), places=13)

    def test_published_golak_table2_c1c3_and_c4_normalization(self):
        values = GOLAK.integrate(12)
        for name, target in GOLAK.PUBLISHED.items():
            with self.subTest(matrix_element=name):
                relative = abs(values[name] - target) / abs(target)
                self.assertLess(relative, 3.0e-4)

    def test_generic_five_angle_projector_matches_golak_closed_integrands(self):
        constants = PWD._OP.N2LOConstants(197.327, 1.29, 92.4, 138.0)
        projector = PWD.FiveAngleProjector(constants)
        order = 4
        closed = GOLAK.integrate(order)
        scalar = projector.project(
            PWD.GOLAK_BETA[1],
            PWD.GOLAK_BETA[1],
            (1.0, 2.0, 3.0, 4.0),
            PWD._OP.N2LOLECs(c1_gev_inverse=-0.81, c3_gev_inverse=-3.4),
            order,
        )
        c4 = projector.project(
            PWD.GOLAK_BETA[2],
            PWD.GOLAK_BETA[1],
            (1.0, 2.0, 3.0, 4.0),
            PWD._OP.N2LOLECs(c4_gev_inverse=3.4),
            order,
        )
        self.assertAlmostEqual((scalar["c1"] + scalar["c3"]).real, closed["G(1,1)"], places=10)
        self.assertAlmostEqual(c4["c4"].real, closed["G(2,1)"], places=10)
        self.assertLess(abs(scalar["c1"].imag + scalar["c3"].imag), 1e-12)
        self.assertLess(abs(c4["c4"].imag), 1e-12)

    def test_generic_contact_projection_has_four_pi_squared_angular_factor(self):
        projector = PWD.FiveAngleProjector(PWD._OP.N2LOConstants.tictac())
        projected = projector.project(
            PWD.GOLAK_BETA[1],
            PWD.GOLAK_BETA[1],
            (0.4, 0.5, 0.6, 0.7),
            PWD._OP.N2LOLECs(c_e=1.0),
            order=2,
        )["cE"]
        e_lec = 1.0 / (projector.constants.f_pi**4 * projector.constants.lambda_chi)
        expected = (4.0 * np.pi) ** 2 * e_lec  # t_pair=1 => tau23=+1
        self.assertAlmostEqual(projected.real, expected, places=10)
        self.assertLess(abs(projected.imag), 1e-12)

    def test_tictac_contact_normalization_is_one_over_four_pi_to_fourth(self):
        projector = PWD.FiveAngleProjector(PWD._OP.N2LOConstants.tictac())
        raw = projector.project(
            PWD.GOLAK_BETA[1],
            PWD.GOLAK_BETA[1],
            (0.4, 0.5, 0.6, 0.7),
            PWD._OP.N2LOLECs(c_e=1.0),
            order=2,
        )
        normalized = projector.to_tictac_normalization(raw)["cE"]
        e_lec = 1.0 / (projector.constants.f_pi**4 * projector.constants.lambda_chi)
        expected = e_lec / (4.0 * np.pi**4)  # t_pair=1 => tau23=+1
        self.assertAlmostEqual(normalized.real, expected, places=12)
        self.assertLess(abs(normalized.imag), 1e-12)

    def test_jj_to_ls_unitary_recoupling_matches_direct_angular_state(self):
        channel = PWD.JjChannel(
            l_pair=1,
            s_pair=1,
            j_pair=1,
            lambda_spectator=1,
            two_j_spectator=1,
            two_total_J=1,
            t_pair=1,
        )
        p_direction = np.array([0.3, 0.4, 0.5])
        q_direction = np.array([-0.2, 0.7, 0.1])
        expansion = channel.ls_expansion()
        self.assertGreater(len(expansion), 1)
        self.assertAlmostEqual(sum(coefficient**2 for _, coefficient in expansion), 1.0, places=13)
        for two_m_j in (-1, 1):
            direct = PWD.angular_spin_state_jj(
                channel, p_direction, q_direction, two_m_j
            )
            recoupled = sum(
                coefficient
                * PWD.angular_spin_state(
                    ls_channel, p_direction, q_direction, two_m_j
                )
                for ls_channel, coefficient in expansion
            )
            np.testing.assert_allclose(direct, recoupled, atol=3e-14, rtol=2e-13)

    def test_direct_jj_projection_matches_unitary_ls_transform(self):
        projector = PWD.FiveAngleProjector(PWD._OP.N2LOConstants.tictac())
        momenta = (0.4, 0.7, 0.9, 1.1)
        lecs = PWD._OP.N2LOLECs(-0.81, -3.2, 5.4, -0.2, -0.205)
        cases = (
            (
                PWD.JjChannel(1, 1, 1, 1, 1, 1, 1),
                PWD.JjChannel(1, 1, 1, 1, 1, 1, 1),
            ),
            (
                PWD.JjChannel(0, 1, 1, 0, 1, 1, 0),
                PWD.JjChannel(0, 0, 0, 0, 1, 1, 1),
            ),
        )
        for bra, ket in cases:
            with self.subTest(bra=bra, ket=ket):
                direct = projector.project_jj_direct(
                    bra, ket, momenta, lecs, order=2
                )
                recoupled = projector.project_jj_recoupled(
                    bra, ket, momenta, lecs, order=2
                )
                for name in direct:
                    np.testing.assert_allclose(
                        direct[name], recoupled[name], atol=8e-13, rtol=3e-13
                    )
        self.assertGreater(
            abs(projector.project_jj_direct(*cases[1], momenta, lecs, order=2)["c4"]),
            1.0,
        )


if __name__ == "__main__":
    unittest.main()
