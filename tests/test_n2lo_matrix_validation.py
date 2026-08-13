"""Regression checks for the checked-in broad N2LO 3NF matrix table."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
import sys
import unittest


REPO = Path(__file__).resolve().parents[1]
GENERATOR = REPO / "tools" / "3nf_oracle" / "generate_n2lo_matrix_table.py"
TABLE = REPO / "output" / "validation" / "n2lo_3nf_matrix_elements.json"


def load_generator():
    spec = importlib.util.spec_from_file_location("n2lo_matrix_table_generator", GENERATOR)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class N2LOMatrixValidationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.generator = load_generator()
        cls.table = json.loads(TABLE.read_text())

    def test_case_manifest_has_broad_bidirectional_coverage(self):
        cases = self.generator.matrix_cases()
        self.assertEqual(len(cases), 26)
        self.assertEqual(len({case.name for case in cases}), len(cases))
        self.assertEqual(len({case.hermitian_pair for case in cases}), 13)
        self.assertEqual({case.direction for case in cases}, {"forward", "reverse"})
        self.assertEqual(
            {case.space for case in cases},
            {"J12_T12", "J12_J32_T12_T32"},
        )
        coverage = {tag for case in cases for tag in case.coverage}
        self.assertTrue({
            "diagonal",
            "off-diagonal",
            "S-wave",
            "P-wave",
            "pair tensor S-D",
            "spectator D-wave",
            "negative parity",
            "J=3/2",
            "T=3/2",
        }.issubset(coverage))

    def test_checked_in_table_passes_all_declared_checks(self):
        report = self.table
        self.assertEqual(report["schema"], "tictac.n2lo_3nf_matrix_validation.v1")
        summary = report["summary"]
        self.assertTrue(summary["passed"])
        self.assertEqual(summary["matrix_elements"], 26)
        self.assertEqual(summary["component_comparisons"], 130)
        self.assertEqual(summary["passed_matrix_elements"], 26)
        self.assertEqual(summary["hermitian_pairs"], 13)
        self.assertEqual(summary["passed_hermitian_pairs"], 13)
        self.assertTrue(all(row["passed"] for row in report["matrix_elements"]))
        self.assertTrue(all(pair["passed"] for pair in report["hermiticity"]))

    def test_every_component_has_signed_nonzero_and_selection_rule_cases(self):
        report = self.table
        for component in self.generator.COMPONENTS:
            with self.subTest(component=component):
                records = [
                    row["components"][component]
                    for row in report["matrix_elements"]
                ]
                nonzero = [
                    item for item in records if item["classification"] == "nonzero"
                ]
                zeros = [
                    item
                    for item in records
                    if item["classification"] == "selection-rule zero"
                ]
                self.assertGreaterEqual(len(nonzero), 8)
                self.assertGreaterEqual(len(zeros), 8)
                signed = [item["production_high_order"]["real"] for item in nonzero]
                self.assertTrue(any(value > 0.0 for value in signed))
                self.assertTrue(any(value < 0.0 for value in signed))

    def test_production_and_separate_factorized_oracle_agree_tightly(self):
        errors = [
            row["components"][component]["absolute_errors"][
                "production_vs_factorized"
            ]
            for row in self.table["matrix_elements"]
            for component in self.generator.COMPONENTS
        ]
        self.assertLess(max(errors), 5.0e-11)


if __name__ == "__main__":
    unittest.main()
