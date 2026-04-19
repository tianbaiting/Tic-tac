#!/usr/bin/env python3
"""L4: Regression tests — 3NF disabled must reproduce 2NF-only baseline."""

import subprocess
import tempfile
import unittest
from pathlib import Path
import hashlib

ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT / "build" / "bin" / "Tic-tac"
BASELINE = ROOT / "tests" / "data" / "baseline_2nf_only"


def run_solver(output_dir, three_nucleon_force="none", c_D=0.0, c_E=0.0):
    cmd = [
        str(SOLVER),
        "Np_WP=5", "Nq_WP=5", "J_2N_max=1", "two_J_3N_max=1",
        "potential_model=N2LOopt",
        f"three_nucleon_force={three_nucleon_force}",
        f"c_D={c_D}", f"c_E={c_E}",
        f"output_folder={output_dir}",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    if result.returncode != 0:
        raise RuntimeError(f"Solver failed: {result.stderr}\n{result.stdout[-500:]}")


def hash_u_files(directory):
    """Hash all U_PW_elements files for bit-for-bit comparison."""
    hashes = {}
    for f in sorted(Path(directory).glob("U_PW_elements_*.txt")):
        content = f.read_bytes()
        hashes[f.name] = hashlib.sha256(content).hexdigest()
    return hashes


class TestRegression(unittest.TestCase):
    def test_none_reproduces_baseline(self):
        """three_nucleon_force=none must produce identical output."""
        with tempfile.TemporaryDirectory() as tmpdir:
            run_solver(tmpdir, "none")
            baseline_hashes = hash_u_files(BASELINE)
            test_hashes = hash_u_files(tmpdir)
            self.assertTrue(len(baseline_hashes) > 0, "No baseline U_PW files found")
            self.assertEqual(baseline_hashes, test_hashes,
                             "three_nucleon_force=none output differs from baseline")


if __name__ == "__main__":
    unittest.main()
