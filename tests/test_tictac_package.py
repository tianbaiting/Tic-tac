#!/usr/bin/env python3
"""Tests for the tictac observable package and the examples/ compatibility shims.

Locks in the Phase 9 migration: the partial-wave amplitude / observable physics
now lives in the importable package ``python/tictac/`` (with ``examples/`` as
thin re-export shims), and the numerical content is unchanged.
"""

import importlib
import math
import sys
import unittest
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PY_PKG = REPO / "python"
EXAMPLES = REPO / "examples"


def _with_path(path):
    sys.path.insert(0, str(path))


class TestTictacPackage(unittest.TestCase):
    def test_public_api_is_importable(self):
        _with_path(PY_PKG)
        import tictac
        for name in ["JPiBlock", "WPGridBin", "assemble_m_matrix",
                     "observables_from_M", "clebsch_gordan", "wigner_6j",
                     "parse_u_file_meta", "select_latest_u_file_family"]:
            self.assertTrue(hasattr(tictac, name), f"missing public symbol {name}")

    def test_observables_module_reexports(self):
        _with_path(PY_PKG)
        from tictac import observables
        self.assertIs(observables.observables_from_M,
                      importlib.import_module("tictac.amplitudes").observables_from_M)

    def test_clebsch_gordan_known_values(self):
        _with_path(PY_PKG)
        from tictac import clebsch_gordan
        # <1/2 +1/2, 1/2 -1/2 | 1 0> = 1/sqrt(2)
        self.assertAlmostEqual(clebsch_gordan(0.5, 0.5, 0.5, -0.5, 1.0, 0.0),
                               1.0 / math.sqrt(2.0), places=12)
        # <1/2 +1/2, 1/2 +1/2 | 1 1> = 1
        self.assertAlmostEqual(clebsch_gordan(0.5, 0.5, 0.5, 0.5, 1.0, 1.0),
                               1.0, places=12)
        # triangle-forbidden -> 0
        self.assertAlmostEqual(clebsch_gordan(0.5, 0.5, 0.5, 0.5, 0.0, 0.0),
                               0.0, places=12)


class TestExamplesShims(unittest.TestCase):
    """The example scripts import `pw_amplitudes` / `solver_u_file_utils` by name;
    the shims must re-export the SAME objects the package defines."""

    def test_pw_amplitudes_shim_matches_package(self):
        _with_path(PY_PKG)
        _with_path(EXAMPLES)
        import pw_amplitudes as shim
        import tictac.amplitudes as pkg
        self.assertIs(shim.assemble_m_matrix, pkg.assemble_m_matrix)
        self.assertIs(shim.observables_from_M, pkg.observables_from_M)
        self.assertIs(shim.JPiBlock, pkg.JPiBlock)
        self.assertIs(shim.parse_u_file_meta,
                      importlib.import_module("tictac.io").parse_u_file_meta)

    def test_solver_u_file_utils_shim_matches_package(self):
        _with_path(PY_PKG)
        _with_path(EXAMPLES)
        import solver_u_file_utils as shim
        import tictac.io as pkg
        self.assertIs(shim.select_latest_u_file_family, pkg.select_latest_u_file_family)
        self.assertIs(shim.parse_u_file_meta, pkg.parse_u_file_meta)


if __name__ == "__main__":
    unittest.main()
