#!/usr/bin/env python3

import json
import shutil
import subprocess
import unittest
from pathlib import Path


class Test190MeVDataPipeline(unittest.TestCase):
    def setUp(self):
        self.repo_root = Path(__file__).resolve().parents[1]
        self.output_dir = Path("/tmp/tic_tac_190mev_pipeline_test")
        if self.output_dir.exists():
            shutil.rmtree(self.output_dir)

    def test_pipeline_generates_data_aligned_outputs(self):
        cmd = [
            "python3",
            "examples/deuteron_proton_Ay.py",
            "--output-dir",
            str(self.output_dir),
            "--grid",
            "experimental",
            "--threshold",
            "1e-12",
        ]

        result = subprocess.run(
            cmd,
            cwd=self.repo_root,
            capture_output=True,
            text=True,
        )

        self.assertEqual(
            result.returncode,
            0,
            msg=f"Pipeline failed.\nstdout:\n{result.stdout}\nstderr:\n{result.stderr}",
        )

        summary_path = self.output_dir / "fit_quality_190MeV.json"
        observables_path = self.output_dir / "fit_observables_190MeV.txt"
        dsigma_path = self.output_dir / "fit_dsigma_190MeV.txt"

        self.assertTrue(summary_path.exists(), "Missing quality summary JSON")
        self.assertTrue(observables_path.exists(), "Missing observable output file")
        self.assertTrue(dsigma_path.exists(), "Missing dsigma output file")

        summary = json.loads(summary_path.read_text(encoding="utf-8"))

        self.assertEqual(summary["status"], "pass")
        self.assertLessEqual(summary["max_abs_error_all"], summary["threshold"])

        for key in ["iT11", "T20", "T21", "T22", "dSigma_dOmega"]:
            self.assertIn(key, summary["metrics"])
            self.assertLessEqual(summary["metrics"][key]["max_abs_error"], summary["threshold"])


if __name__ == "__main__":
    unittest.main()
