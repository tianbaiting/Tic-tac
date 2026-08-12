#!/usr/bin/env python3
"""Tests for the complete-3NF reduced-solver comparison report."""

from __future__ import annotations

import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
MODULE_PATH = REPO / "tools" / "3nf_oracle" / "compare_reduced_solver_outputs.py"
SPEC = importlib.util.spec_from_file_location("compare_reduced_solver_outputs", MODULE_PATH)
CROSSCHECK = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = CROSSCHECK
SPEC.loader.exec_module(CROSSCHECK)


class TestReducedSolverCrosscheck(unittest.TestCase):
    def test_agreement_and_honesty_are_reported_separately(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            pade = root / "pade"
            dense = root / "dense"
            pade.mkdir()
            dense.mkdir()

            common = (
                "Running program for:\n"
                "Np_WP: 4\n"
                "Energy input file: /checkout/CPP/Input/lab_energies.txt\n"
                f"Cache root: {root / 'cache'}\n"
                f"P123-matrix read/write folder: {root / 'cache' / 'p123'}\n"
            )
            (pade / "run_parameters.txt").write_text(
                common + "Solve Faddeev with LAPACK: false\n", encoding="utf-8"
            )
            (dense / "run_parameters.txt").write_text(
                common + "Solve Faddeev with LAPACK: true\n", encoding="utf-8"
            )

            filename = "U_PW_elements_Np_4_Nq_3_JP_1_1_Jmax_1_PSI_0.txt"
            header = "# Deuteron BE: -0.011400617 MeV\n"
            (pade / filename).write_text(
                header + "1.0 0.5 0 +1.0000000000e+00+2.0000000000e+00j\n",
                encoding="utf-8",
            )
            (dense / filename).write_text(
                header + "1.0 0.5 0 +1.0000000010e+00+2.0000000000e+00j\n",
                encoding="utf-8",
            )
            convergence = pade / filename.replace(
                "U_PW_elements_", "U_PW_convergence_"
            ).replace("_PSI_0", "")
            convergence.write_text("0 0 0 2 8\n", encoding="utf-8")

            report = CROSSCHECK.build_report(pade, dense)
            self.assertAlmostEqual(
                report["aggregate"]["max_absolute_difference_mev"], 1.0e-9,
                delta=2.0e-16,
            )
            self.assertFalse(report["aggregate"]["all_pade_elements_truly_converged"])
            self.assertEqual(report["files"][0]["pade_status"]["num_max_order_truncated"], 1)
            self.assertEqual(report["parameters"]["Energy input file"], "CPP/Input/lab_energies.txt")
            self.assertNotIn("Cache root", report["parameters"])


if __name__ == "__main__":
    unittest.main()
