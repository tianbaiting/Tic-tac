#!/usr/bin/env python3
"""Regression guard for complete-N2LO W1 wave-packet cell integration."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO / "tools" / "3nf_oracle" / "wp_quadrature_convergence.py"
SPEC = importlib.util.spec_from_file_location("wp_quadrature_convergence", MODULE_PATH)
WP = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = WP
SPEC.loader.exec_module(WP)


@unittest.skipUnless(WP.DRIVER.is_file(), "build print_w1_element before running this test")
class TestFactorizedWpQuadrature(unittest.TestCase):
    def test_wide_c4_cd_transition_requires_more_than_two_points(self):
        """The difficult complete-force cell converges at N=4, not N=1 or N=2."""
        bounds = (0.2, 0.8, 0.2, 0.8, 0.2, 0.8, 0.2, 0.8)
        options = dict(
            c4=5.4,
            factorized=True,
            transfer_order=6,
            alpha_r=4,
            alpha_c=0,
        )
        values = {
            order: WP.bin_average(
                -0.205, -0.2, -0.81, -3.2, 500.0,
                bounds, order, **options,
            )
            for order in (1, 2, 4, 6)
        }
        reference = values[6]
        relative_error = {
            order: abs(values[order] - reference) / abs(reference)
            for order in (1, 2, 4)
        }

        self.assertGreater(relative_error[1], 0.15)
        self.assertGreater(relative_error[2], 3.0e-3)
        self.assertLess(relative_error[4], 1.0e-5)
        self.assertAlmostEqual(reference / 3.184203956289e-2, 1.0, places=10)


if __name__ == "__main__":
    unittest.main()
