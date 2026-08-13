#!/usr/bin/env python3
"""Integration checks for the P123-free deuteron binding preflight."""

from __future__ import annotations

import json
import subprocess
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOLVER = ROOT / "build" / "bin" / "Tic-tac"


@unittest.skipUnless(SOLVER.is_file(), "build/bin/Tic-tac is not built")
class TestDeuteronBindingMode(unittest.TestCase):
    def run_solver(self, output_dir: Path, *extra: str) -> subprocess.CompletedProcess[str]:
        command = [
            str(SOLVER),
            "Np_WP=5",
            "Nq_WP=1",
            "J_2N_max=1",
            "two_J_3N_max=1",
            "Np_per_WP=8",
            "Nq_per_WP=1",
            "P123_omp_num_threads=1",
            "potential_model=N2LOopt",
            "three_nucleon_force=none",
            f"output_folder={output_dir}",
            f"cache_root={output_dir / 'cache'}",
            *extra,
        ]
        return subprocess.run(
            command,
            cwd=ROOT,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )

    def test_writes_machine_result_without_three_body_outputs(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            result = self.run_solver(
                output_dir,
                "solve_faddeev=false",
                "calculate_and_store_P123=false",
                "deuteron_binding_only=true",
            )
            self.assertEqual(result.returncode, 0, result.stderr + result.stdout[-1000:])

            result_path = output_dir / "deuteron_binding.json"
            self.assertTrue(result_path.is_file())
            payload = json.loads(result_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema"], "tictac.deuteron_binding.v1")
            self.assertEqual(payload["calculation_path"], "production_V_WP_to_SWP")
            self.assertEqual(payload["potential_model"], "N2LOopt")
            self.assertEqual(payload["Np_WP"], 5)
            self.assertLess(payload["energy_mev"], 0.0)
            self.assertAlmostEqual(
                payload["binding_magnitude_mev"], -payload["energy_mev"], places=14
            )
            self.assertAlmostEqual(payload["energy_mev"], -0.0200261520, places=9)

            self.assertFalse(list(output_dir.glob("U_PW*")))
            self.assertFalse(list(output_dir.glob("*P123*")))
            self.assertNotIn("Working on 3N-channel", result.stdout)
            self.assertNotIn("Constructing P123", result.stdout)

    def test_rejects_implicit_full_solver_defaults(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            result = self.run_solver(output_dir, "deuteron_binding_only=true")
            combined_output = result.stdout + result.stderr
            self.assertIn(
                "requires solve_faddeev=false and calculate_and_store_P123=false",
                combined_output,
            )
            self.assertFalse((output_dir / "deuteron_binding.json").exists())

    def test_rejects_three_nucleon_force(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            output_dir = Path(temporary_directory)
            result = self.run_solver(
                output_dir,
                "solve_faddeev=false",
                "calculate_and_store_P123=false",
                "deuteron_binding_only=true",
                "three_nucleon_force=chiral_N2LO_full_factorized",
            )
            combined_output = result.stdout + result.stderr
            self.assertIn("requires three_nucleon_force=none", combined_output)
            self.assertFalse((output_dir / "deuteron_binding.json").exists())

    def test_independent_p_and_q_chebyshev_controls(self) -> None:
        binding_flags = (
            "solve_faddeev=false",
            "calculate_and_store_P123=false",
            "deuteron_binding_only=true",
            "Nq_WP=4",
        )

        with tempfile.TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)

            def run_case(name: str, *grid_flags: str) -> tuple[Path, dict]:
                output_dir = root / name
                output_dir.mkdir()
                result = self.run_solver(output_dir, *binding_flags, *grid_flags)
                self.assertEqual(
                    result.returncode, 0, result.stderr + result.stdout[-1000:]
                )
                payload = json.loads(
                    (output_dir / "deuteron_binding.json").read_text(encoding="utf-8")
                )
                return output_dir, payload

            common_dir, common = run_case("common", "chebyshev_s=100")
            split_dir, split = run_case(
                "split", "chebyshev_s=100", "p_chebyshev_s=300", "q_chebyshev_s=100"
            )
            inherited_dir, inherited = run_case(
                "inherited", "chebyshev_s=300", "q_chebyshev_s=100"
            )
            q_changed_dir, q_changed = run_case(
                "q_changed",
                "chebyshev_s=100",
                "p_chebyshev_s=300",
                "q_chebyshev_s=300",
            )

            self.assertNotAlmostEqual(
                common["energy_mev"], split["energy_mev"], places=12
            )
            self.assertAlmostEqual(
                split["energy_mev"], inherited["energy_mev"], places=14
            )
            self.assertAlmostEqual(
                split["energy_mev"], q_changed["energy_mev"], places=14
            )

            self.assertEqual(split["p_chebyshev_s"], 300)
            self.assertEqual(split["q_chebyshev_s"], 100)
            self.assertEqual(inherited["p_chebyshev_s"], 300)
            self.assertEqual(inherited["q_chebyshev_s"], 100)

            q_filename = "q_kinematics_Nq_4.txt"
            self.assertEqual(
                (split_dir / q_filename).read_bytes(),
                (inherited_dir / q_filename).read_bytes(),
            )
            self.assertNotEqual(
                (split_dir / q_filename).read_bytes(),
                (q_changed_dir / q_filename).read_bytes(),
            )


if __name__ == "__main__":
    unittest.main()
