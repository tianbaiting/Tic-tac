"""Tests for the fail-closed low-energy nd Ay acceptance audit."""

from __future__ import annotations

from dataclasses import replace
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np


REPO = Path(__file__).resolve().parents[1]
AUDITOR_PATH = REPO / "examples" / "audit_low_energy_Ay.py"
PREFLIGHT_PATH = (
    REPO / "output" / "validation" / "n2lo_3nf_low_energy_Ay_preflight.json"
)


def load_auditor():
    spec = importlib.util.spec_from_file_location("low_energy_ay_auditor", AUDITOR_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


class LowEnergyAyAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.audit = load_auditor()

    def make_artifact_dir(self, root: Path, model: str) -> Path:
        root.mkdir(parents=True)
        (root / "run_parameters.txt").write_text(
            "\n".join([
                "Running program for:",
                "two_J_3N_max:                    3",
                "Np_WP:                           24",
                "Nq_WP:                           24",
                "J_2N_max:                        3",
                "Nphi:                            48",
                "Nx:                              48",
                "chebyshev sparseness:            1.000",
                "chebyshev scale:                 100.000",
                "Np_per_WP:                       8",
                "Nq_per_WP:                       8",
                "Np_per_WP_W1:                    4",
                "Nq_per_WP_W1:                    4",
                "Nangle_3NF:                      10",
                "Padé maximum diagonal order:     24",
                "Tensor-force on:                 true",
                "Isospin-breaking in 1S0:         true",
                "Mid-point approximation:         false",
                "Calculate breakup amplitudes:    false",
                "Potential model:                 Idaho_N3LO",
                f"Three-nucleon force:             {model}",
                "3NF LEC c_D:                     -0.200",
                "3NF LEC c_E:                     -0.205",
                "3NF cutoff Lambda_3NF [MeV]:     500.000",
                "p-momentum grid type:            chebyshev",
                "p-momentum grid input file:      ",
                "q-momentum grid type:            chebyshev",
                "q-momentum grid input file:      ",
                "Energy input file:               energies.txt",
                "",
            ]),
            encoding="utf-8",
        )
        for two_j in (1, 3):
            for parity in (-1, 1):
                suffix = f"_Np_24_Nq_24_JP_{two_j}_{parity}_Jmax_3"
                (root / f"U_PW_elements{suffix}_PSI_0.txt").write_text(
                    "# Deuteron BE: -2.22457000 MeV\n", encoding="utf-8"
                )
                (root / f"U_PW_convergence{suffix}.txt").write_text(
                    "# row col q Conv order\n0 0 0 1 14\n", encoding="utf-8"
                )
        return root

    def test_inspect_run_parses_binding_blocks_and_pade(self):
        with tempfile.TemporaryDirectory() as temporary:
            run_dir = self.make_artifact_dir(Path(temporary) / "run", "none")
            result = self.audit.inspect_run(run_dir)
        self.assertEqual(result.parameters["np_wp"], 24)
        self.assertEqual(result.parameters["pade_max_order"], 24)
        self.assertAlmostEqual(result.binding_energy_mev, -2.22457)
        self.assertTrue(result.binding_headers_consistent)
        self.assertEqual(result.missing_blocks, [])
        self.assertEqual(result.pade_counts, {
            "converged": 4, "truncated": 0, "unknown": 0,
        })

    def test_core_gates_accept_only_matched_complete_model_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            two_nf = self.audit.inspect_run(
                self.make_artifact_dir(root / "two", "none")
            )
            three_nf = self.audit.inspect_run(
                self.make_artifact_dir(
                    root / "three", self.audit.COMPLETE_3NF_MODEL
                )
            )
            gates = self.audit.core_gates(two_nf, three_nf, 0.02)
            self.assertTrue(all(item["passed"] for item in gates.values()))

            approximate = replace(
                three_nf,
                parameters=three_nf.parameters | {
                    "three_nucleon_force": "chiral_N2LO_c1c3cDcE_approx"
                },
            )
            self.assertFalse(
                self.audit.core_gates(two_nf, approximate, 0.02)[
                    "model_identity"
                ]["passed"]
            )

    def test_curve_metrics_include_peak_rmse_and_chi2(self):
        metrics = self.audit.curve_metrics(
            np.asarray([0.0, 90.0, 180.0]),
            np.asarray([0.0, 0.20, 0.0]),
            np.asarray([90.0]),
            np.asarray([0.25]),
            np.asarray([0.05]),
        )
        self.assertAlmostEqual(metrics["ay_max"], 0.20)
        self.assertAlmostEqual(metrics["theta_at_ay_max_deg"], 90.0)
        self.assertAlmostEqual(metrics["peak_deficit_model_minus_experiment"], -0.05)
        self.assertAlmostEqual(metrics["rmse"], 0.05)
        self.assertAlmostEqual(metrics["chi2"], 1.0)
        self.assertAlmostEqual(metrics["maximum_absolute_residual"], 0.05)

    def test_convergence_reference_must_change_only_requested_dimension(self):
        with tempfile.TemporaryDirectory() as temporary:
            run = self.audit.inspect_run(
                self.make_artifact_dir(Path(temporary) / "main", "none")
            )
        np_reference = replace(
            run, parameters=run.parameters | {"np_wp": 20}
        )
        np_changes = self.audit.ladder_parameter_changes(run, np_reference)
        self.assertTrue(self.audit.valid_ladder_change(np_changes, "Np"))
        self.assertFalse(self.audit.valid_ladder_change(np_changes, "Nq"))
        self.assertFalse(self.audit.valid_ladder_change({}, "Np"))

        contaminated = replace(
            run,
            parameters=run.parameters | {"np_wp": 20, "j_2n_max": 2},
        )
        contaminated_changes = self.audit.ladder_parameter_changes(
            run, contaminated
        )
        self.assertFalse(
            self.audit.valid_ladder_change(contaminated_changes, "Np")
        )

    def test_checked_preflight_withholds_unconverged_interpretation(self):
        report = json.loads(PREFLIGHT_PATH.read_text(encoding="utf-8"))
        self.assertFalse(report["accepted_physics_benchmark"])
        self.assertFalse(report["interpretation_allowed"])
        failed = {
            name for name, item in report["gates"].items() if not item["passed"]
        }
        self.assertEqual(failed, {
            "model_identity",
            "deuteron_binding",
            "pade_honesty",
            "Ay_convergence_ladders",
            "effect_resolved_above_numerical_uncertainty",
        })
        self.assertTrue(report["gates"]["paired_numerical_setup"]["passed"])
        self.assertTrue(
            report["gates"]["partial_wave_block_completeness"]["passed"]
        )


if __name__ == "__main__":
    unittest.main()
