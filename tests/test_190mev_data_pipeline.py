#!/usr/bin/env python3

import math
import tempfile
import unittest
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = ROOT / "examples"
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from compare_Ay_experiment import parse_u_file, read_experimental_iT11, read_experimental_dsigma  # noqa: E402


class Test190MeVSolverValidationHelpers(unittest.TestCase):
    def test_parse_solver_u_file_and_proxy(self):
        content = (
            "# header\n"
            "1.0000000000000000e+02  7.0000000000000000e+01        14   "
            "+1.0000000000000000e+00+0.0000000000000000e+00j   "
            "+0.0000000000000000e+00+1.0000000000000000e+00j   "
            "+0.0000000000000000e+00+1.0000000000000000e+00j   "
            "+1.0000000000000000e+00+0.0000000000000000e+00j\n"
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "U_PW_elements_Np_20_Nq_20_JP_1_1_Jmax_2_PSI_0.txt"
            fpath.write_text(content, encoding="utf-8")

            points = parse_u_file(fpath)
            self.assertEqual(len(points), 1)
            self.assertEqual(points[0].parity, "+")

            # f_no_flip = 2+0j, f_flip = 0+2j => ay_proxy = Im(conj(f0)*f1)/(|f0|^2+|f1|^2)=4/8=0.5
            self.assertTrue(math.isclose(points[0].ay_proxy, 0.5, rel_tol=0.0, abs_tol=1e-12))
            self.assertTrue(math.isclose(points[0].dsigma_proxy, 8.0, rel_tol=0.0, abs_tol=1e-12))

    def test_experiment_data_parsers(self):
        it11 = read_experimental_iT11(ROOT / "data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt")
        dsig = read_experimental_dsigma(ROOT / "data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt")

        self.assertGreater(len(it11["values"]), 10)
        self.assertGreater(len(dsig["values"]), 10)


if __name__ == "__main__":
    unittest.main()
