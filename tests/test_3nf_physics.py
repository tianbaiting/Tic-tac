#!/usr/bin/env python3
"""L5: Physics validation — 3NF should improve Ay agreement with experiment."""

import subprocess
import tempfile
import unittest
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = ROOT / "examples"
DATA_DIR = ROOT / "data" / "DataOfCrosssectionAndPol"

if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))


class TestPhysicsValidation(unittest.TestCase):
    """
    Run 190 MeV benchmark with and without 3NF, compare Ay/iT11 to experiment.
    This test is slow (minutes) and should be run explicitly, not in CI.
    """

    @unittest.skipUnless(
        (ROOT / "build" / "bin" / "Tic-tac").exists(),
        "Solver not built"
    )
    def test_3nf_improves_ay(self):
        """With 3NF enabled, Ay discrepancy should be reduced vs 2NF-only."""
        # This test runs the full pipeline via examples/deuteron_proton_Ay.py
        # with and without 3NF, and compares RMSE against experiment.
        #
        # Implementation depends on the observable extraction code in
        # examples/compare_Ay_experiment.py.
        #
        # For now, this is a placeholder that documents the intended test.
        # The full implementation will be added once the 1PE-CT and 2PE terms
        # are verified to produce correct matrix elements.
        self.skipTest("Full physics validation requires completed 3NF implementation")


if __name__ == "__main__":
    unittest.main()
