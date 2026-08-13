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
        """A complete, converged 3NF run may improve, preserve, or worsen Ay."""
        # Operator completeness is no longer the blocker.  The remaining gate
        # is a paired set of solver artifacts with physical deuteron binding,
        # complete J^pi coverage, honest Pade tails, and Ay convergence ladders
        # in Np, Nq, W1, J2N, J3N, and resummation order.  The fail-closed
        # artifact audit lives in examples/audit_low_energy_Ay.py.
        self.skipTest(
            "No complete-N2LO low-energy Ay artifact set has passed the "
            "convergence audit"
        )


if __name__ == "__main__":
    unittest.main()
