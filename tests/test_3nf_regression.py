#!/usr/bin/env python3
"""L4: Regression tests — 3NF disabled must reproduce 2NF-only baseline.

The production solver uses OpenMP reductions whose floating-point sum order is
non-deterministic across runs / thread counts, so the U_PW_elements output files
are NOT bit-for-bit reproducible (the run-to-run max|delta| is ~1e-10 at the
Np=Nq=5 grid used here — the OpenMP non-determinism floor). This test therefore
compares the U-matrix elements NUMERICALLY (max abs diff < tol) rather than by
SHA-256 hash. The baseline was generated with three_nucleon_force=none on the
master branch.

Also verifies the Phase 7 parser fix: `energy_input_file=<path>.txt` must be
parseable as an inline key=value override (previously misparsed as an input-list
file because the arg ended in .txt).
"""

import math
import os
import re
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT / "build" / "bin" / "Tic-tac"
BASELINE = ROOT / "tests" / "data" / "baseline_2nf_only"
ENERGY_FILE = ROOT / "CPP" / "Input" / "lab_energies.txt"

U_TOL = 1e-7

# matches a complex token like +7.625e-04-1.903e-03j or -1.21e-01+4.04e-01j
_CPLX_RE = re.compile(r"([+-]?\d+\.\d+[eE][+-]?\d+)([+-]\d+\.\d+[eE][+-]?\d+)j")


def run_solver(output_dir, three_nucleon_force="none"):
    cmd = [
        str(SOLVER),
        "Np_WP=5", "Nq_WP=5", "J_2N_max=1", "two_J_3N_max=1",
        "potential_model=N2LOopt",
        f"three_nucleon_force={three_nucleon_force}",
        "c_D=0.0", "c_E=0.0",
        f"energy_input_file={ENERGY_FILE}",
        f"output_folder={output_dir}",
    ]
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = "2"
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=120, env=env)
    if result.returncode != 0:
        raise RuntimeError(f"Solver failed: {result.stderr}\n{result.stdout[-500:]}")


def parse_u_file(path):
    """Parse a U_PW_elements_*.txt data rows into (Tlab, Ecm, q_idx, [complex...])."""
    rows = []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) < 4:
            continue
        try:
            tlab = float(parts[0]); ecm = float(parts[1]); qidx = int(parts[2])
        except ValueError:
            continue
        us = []
        for tok in parts[3:]:
            m = _CPLX_RE.match(tok)
            if m:
                us.append(complex(float(m.group(1)), float(m.group(2))))
            else:
                try:
                    us.append(float(tok))
                except ValueError:
                    break
        rows.append((tlab, ecm, qidx, us))
    return rows


class TestRegression(unittest.TestCase):
    def _diff_against_baseline(self, test_dir):
        max_diff = 0.0
        n_compared = 0
        for bf in sorted(BASELINE.glob("U_PW_elements_*.txt")):
            tf = Path(test_dir) / bf.name
            self.assertTrue(tf.exists(), f"missing output file {tf.name}")
            base_rows = parse_u_file(bf)
            test_rows = parse_u_file(tf)
            self.assertEqual(len(base_rows), len(test_rows),
                             f"{bf.name}: row count mismatch {len(base_rows)} vs {len(test_rows)}")
            for (bt, be, bq, bu), (tt, te, tq, tu) in zip(base_rows, test_rows):
                self.assertEqual(bq, tq, f"{bf.name}: q_idx mismatch")
                self.assertEqual(len(bu), len(tu), f"{bf.name}: U count mismatch")
                for bv, tv in zip(bu, tu):
                    d = abs(bv - tv)
                    if d > max_diff:
                        max_diff = d
                    n_compared += 1
        self.assertTrue(n_compared > 0, "no U-matrix elements compared")
        return max_diff, n_compared

    def test_none_reproduces_baseline_within_tolerance(self):
        """three_nucleon_force=none must reproduce the 2NF-only baseline
        to within the OpenMP non-determinism floor (NOT bit-for-bit)."""
        with tempfile.TemporaryDirectory() as tmpdir:
            run_solver(tmpdir, "none")
            max_diff, n = self._diff_against_baseline(tmpdir)
            print(f"  [regression] {n} U-elements compared; max|Δ| vs baseline = {max_diff:.3e}")
            self.assertLess(max_diff, U_TOL,
                f"three_nucleon_force=none differs from baseline by {max_diff:.3e} "
                f"(tol {U_TOL:.0e}). This indicates a REAL regression in the 2NF-only "
                f"path (not OpenMP noise, which is ~1e-10).")

    def test_parser_energy_input_file_inline(self):
        """Phase 7: energy_input_file=<path>.txt must be parseable inline."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # if the parser bug is present, this raises RuntimeError because
            # the solver exits with 'Unable to open file lab_energies.txt'.
            run_solver(tmpdir, "none")
            # reaching here means the parser accepted the inline override
            self.assertTrue((Path(tmpdir) / "run_parameters.txt").exists())


if __name__ == "__main__":
    unittest.main()
