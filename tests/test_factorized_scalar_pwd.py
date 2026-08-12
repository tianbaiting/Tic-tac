"""Independent checks of the Hebeler three-integral scalar PWD building block."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path
import sys
import unittest

import numpy as np


REPO = Path(__file__).resolve().parents[1]
ORACLE_DIR = REPO / "tools" / "3nf_oracle"

PWD_PATH = ORACLE_DIR / "full_vector_five_angle_pwd.py"
PWD_SPEC = importlib.util.spec_from_file_location("full_vector_five_angle_pwd", PWD_PATH)
PWD = importlib.util.module_from_spec(PWD_SPEC)
assert PWD_SPEC.loader is not None
sys.modules[PWD_SPEC.name] = PWD
PWD_SPEC.loader.exec_module(PWD)

FACTORIZED_PATH = ORACLE_DIR / "factorized_scalar_pwd.py"
FACTORIZED_SPEC = importlib.util.spec_from_file_location(
    "factorized_scalar_pwd", FACTORIZED_PATH
)
FACTORIZED = importlib.util.module_from_spec(FACTORIZED_SPEC)
assert FACTORIZED_SPEC.loader is not None
sys.modules[FACTORIZED_SPEC.name] = FACTORIZED
FACTORIZED_SPEC.loader.exec_module(FACTORIZED)


def direct_five_angle_scalar_ls(
    bra: PWD.LSChannel,
    ket: PWD.LSChannel,
    momenta: tuple[float, float, float, float],
    scalar_kernel,
    order: int,
) -> complex:
    """Structurally independent Golak five-angle projection of a scalar kernel."""

    if bra.two_total_J != ket.two_total_J:
        return 0.0j
    p, q, pp, qp = momenta
    cosines, cosine_weights = PWD._nodes_weights(order, -1.0, 1.0)
    phis, phi_weights = PWD._nodes_weights(order, 0.0, 2.0 * math.pi)
    sines = np.sqrt(np.maximum(0.0, 1.0 - cosines * cosines))
    p_in = np.array([0.0, 0.0, p])
    m_j_values = tuple(PWD._m_values(bra.two_total_J))
    total = 0.0j

    for cosine_q, sine_q, weight_q in zip(cosines, sines, cosine_weights):
        q_in = q * np.array([sine_q, 0.0, cosine_q])
        ket_states = {
            m_j: PWD.angular_spin_state(ket, p_in, q_in, m_j)
            for m_j in m_j_values
        }
        for cosine_pp, sine_pp, weight_pp in zip(cosines, sines, cosine_weights):
            for phi_pp, weight_phi_pp in zip(phis, phi_weights):
                p_out = pp * np.array(
                    [sine_pp * math.cos(phi_pp), sine_pp * math.sin(phi_pp), cosine_pp]
                )
                for cosine_qp, sine_qp, weight_qp in zip(
                    cosines, sines, cosine_weights
                ):
                    for phi_qp, weight_phi_qp in zip(phis, phi_weights):
                        q_out = qp * np.array(
                            [
                                sine_qp * math.cos(phi_qp),
                                sine_qp * math.sin(phi_qp),
                                cosine_qp,
                            ]
                        )
                        bra_states = {
                            m_j: PWD.angular_spin_state(bra, p_out, q_out, m_j)
                            for m_j in m_j_values
                        }
                        delta_p = p_out - p_in
                        delta_q = q_out - q_in
                        magnitude_p = float(np.linalg.norm(delta_p))
                        magnitude_q = float(np.linalg.norm(delta_q))
                        relative_cosine = float(
                            np.dot(delta_p, delta_q) / (magnitude_p * magnitude_q)
                        )
                        spin_orbital_overlap = sum(
                            np.vdot(bra_states[m_j], ket_states[m_j])
                            for m_j in m_j_values
                        ) / (bra.two_total_J + 1.0)
                        total += (
                            weight_q
                            * weight_pp
                            * weight_phi_pp
                            * weight_qp
                            * weight_phi_qp
                            * scalar_kernel(magnitude_p, magnitude_q, relative_cosine)
                            * spin_orbital_overlap
                        )

    return 8.0 * math.pi**2 * total


def polynomial_local_kernel(delta_p: float, delta_q: float, cosine: float) -> float:
    """Finite-rank scalar kernel with monopole, dipole, and quadrupole pieces."""
    return (
        1.0
        + 0.11 * delta_p**2
        + 0.07 * delta_q**2
        + 0.05 * delta_p * delta_q * cosine
        + 0.03
        * (delta_p * delta_q) ** 2
        * 0.5
        * (3.0 * cosine**2 - 1.0)
    )


class FactorizedScalarPWDTests(unittest.TestCase):
    def test_constant_contact_normalization_and_s_wave_selection(self):
        momenta = (0.37, 0.52, 0.81, 0.66)
        s_wave = FACTORIZED.project_scalar_ls(
            PWD.GOLAK_BETA[1],
            PWD.GOLAK_BETA[1],
            momenta,
            lambda _dp, _dq, _x: 1.0,
            order=8,
        )
        p_wave = FACTORIZED.project_scalar_ls(
            PWD.GOLAK_BETA[4],
            PWD.GOLAK_BETA[4],
            momenta,
            lambda _dp, _dq, _x: 1.0,
            order=8,
        )
        self.assertAlmostEqual(s_wave.real, (4.0 * math.pi) ** 2, places=10)
        self.assertLess(abs(s_wave.imag), 1.0e-12)
        self.assertLess(abs(p_wave), 1.0e-12)

    def test_finite_rank_kernel_matches_independent_five_angle_projection(self):
        momenta = (0.43, 0.61, 0.78, 0.69)
        cases = (
            (PWD.GOLAK_BETA[4], PWD.GOLAK_BETA[4]),
            (
                PWD.GOLAK_BETA[3],
                PWD.LSChannel(2, 1, 0, 2, 3, 1, 0),
            ),
        )
        for bra, ket in cases:
            with self.subTest(bra=bra, ket=ket):
                factorized = FACTORIZED.project_scalar_ls(
                    bra, ket, momenta, polynomial_local_kernel, order=8
                )
                direct = direct_five_angle_scalar_ls(
                    bra, ket, momenta, polynomial_local_kernel, order=6
                )
                self.assertLess(abs(factorized.imag), 2.0e-11)
                self.assertLess(abs(direct.imag), 2.0e-11)
                self.assertAlmostEqual(factorized.real, direct.real, delta=2.0e-6)


if __name__ == "__main__":
    unittest.main()
