#!/usr/bin/env python3

import json
import math
import subprocess
import tempfile
import unittest
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[1]
EXAMPLES_DIR = ROOT / "examples"
if str(EXAMPLES_DIR) not in sys.path:
    sys.path.insert(0, str(EXAMPLES_DIR))

from compare_Ay_experiment import read_experimental_iT11, read_experimental_dsigma  # noqa: E402
from deuteron_proton_Ay import SolverRunConfig as SingleEnergyRunConfig, write_solver_inputs as write_single_inputs  # noqa: E402
from run_dpol_p_observables import SolverRunConfig as MultiEnergyRunConfig, write_solver_inputs as write_multi_inputs  # noqa: E402
from pw_amplitudes import parse_u_file as pw_parse_u_file  # noqa: E402
from tlab_utils import format_tlab_dir_name  # noqa: E402


_MINIMAL_U_CONTENT = (
    "# Elastic Nd-scattering U-matrix elements\n"
    "#      Name   row-idx   col-idx        l'      2*j'         l       2*j\n"
    "        U00         0         0         0         1         0         1\n"
    "        U01         0         1         0         1         2         3\n"
    "        U10         1         0         2         3         0         1\n"
    "        U11         1         1         2         3         2         3\n"
    "# header\n"
    "1.0000000000000000e+02  7.0000000000000000e+01        14   "
    "+1.0000000000000000e+00+0.0000000000000000e+00j   "
    "+0.0000000000000000e+00+1.0000000000000000e+00j   "
    "+0.0000000000000000e+00+1.0000000000000000e+00j   "
    "+1.0000000000000000e+00+0.0000000000000000e+00j\n"
)

_MINIMAL_Q_KINEMATICS = (
    "# BOUNDARIES:\n"
    "# idx  q[MeV]  Ecm[MeV]  Tlab[MeV]\n"
    "  14  +2.50e+02  +5.0e+01  +7.5e+01\n"
    "  15  +3.00e+02  +7.0e+01  +1.05e+02\n"
    "  16  +3.50e+02  +9.5e+01  +1.40e+02\n"
    "# BIN MID-POINTS:\n"
    "# idx  q[MeV]  Ecm[MeV]  Tlab[MeV]\n"
    "  14  +2.75e+02  +6.0e+01  +9.0e+01\n"
    "  15  +3.25e+02  +8.0e+01  +1.20e+02\n"
)


class Test190MeVSolverValidationHelpers(unittest.TestCase):
    def test_pw_parse_u_file_extracts_jpi_block(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fpath = Path(tmpdir) / "U_PW_elements_Np_20_Nq_20_JP_1_1_Jmax_2_PSI_0.txt"
            fpath.write_text(_MINIMAL_U_CONTENT, encoding="utf-8")

            block = pw_parse_u_file(fpath)
            self.assertIsNotNone(block)
            self.assertEqual(block.two_J, 1)
            self.assertEqual(block.parity, 1)
            self.assertEqual(len(block.channels), 2)
            self.assertEqual((block.channels[0].l, block.channels[0].two_j), (0, 1))
            self.assertEqual((block.channels[1].l, block.channels[1].two_j), (2, 3))
            self.assertEqual(len(block.points), 1)
            point = block.points[0]
            self.assertEqual(point.q_idx, 14)
            self.assertAlmostEqual(point.tlab, 100.0)
            self.assertEqual(point.matrix.shape, (2, 2))
            self.assertEqual(complex(point.matrix[0, 0]), 1 + 0j)
            self.assertEqual(complex(point.matrix[0, 1]), 0 + 1j)
            self.assertEqual(complex(point.matrix[1, 0]), 0 + 1j)
            self.assertEqual(complex(point.matrix[1, 1]), 1 + 0j)

    def test_experiment_data_parsers(self):
        it11 = read_experimental_iT11(ROOT / "data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt")
        dsig = read_experimental_dsigma(ROOT / "data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt")

        self.assertGreater(len(it11["values"]), 10)
        self.assertGreater(len(dsig["values"]), 10)

    def test_single_energy_input_generation_uses_tlab_mev_naming(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "CPP" / "Input").mkdir(parents=True, exist_ok=True)

            cfg = SingleEnergyRunConfig(
                root=root,
                work_dir=root / "output" / "single",
                solver=root / "CPP" / "run",
                target_tlab_mev=190.0,
                two_j_3n_max=1,
                j_2n_max=2,
                np_wp=20,
                nq_wp=20,
                nphi=24,
                nx=24,
                np_per_wp=6,
                nq_per_wp=6,
                chebyshev_s=150.0,
                chebyshev_t=1.0,
                threads=4,
                potential_model="LO_internal",
                timeout_s=1800,
                calculate_p123=True,
            )

            energy_file, config_file, solver_out_dir, log_file = write_single_inputs(cfg)

            self.assertEqual(energy_file.name, "tlab_190MeV_validation_autogen.txt")
            self.assertEqual(config_file.name, "input_tlab_190MeV_solver_validation_autogen.txt")
            self.assertEqual(energy_file.read_text(encoding="utf-8"), "190.000000\n")
            self.assertTrue(solver_out_dir.exists())
            self.assertEqual(log_file.name, "solver_run.log")

            config_text = config_file.read_text(encoding="utf-8")
            self.assertIn("energy_input_file=Input/tlab_190MeV_validation_autogen.txt", config_text)

    def test_multi_energy_input_generation_uses_tlab_mev_schema(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "CPP" / "Input").mkdir(parents=True, exist_ok=True)

            cfg = MultiEnergyRunConfig(
                root=root,
                work_dir=root / "output" / "multi",
                solver=root / "CPP" / "run",
                target_tlabs_mev=[70.0, 135.0, 190.0],
                two_j_3n_max=1,
                j_2n_max=2,
                np_wp=20,
                nq_wp=22,
                nphi=24,
                nx=24,
                np_per_wp=6,
                nq_per_wp=6,
                chebyshev_s=180.0,
                chebyshev_t=1.0,
                threads=4,
                potential_model="LO_internal",
                timeout_s=3600,
                calculate_p123=True,
                reuse_p123_from=None,
                angle_min_deg=20.0,
                angle_max_deg=170.0,
                angle_step_deg=1.0,
            )

            target_file, energy_file, config_file, solver_out_dir, p123_dir = write_multi_inputs(cfg)

            expected_payload = "70.000000\n135.000000\n190.000000\n"
            self.assertEqual(target_file.name, "target_tlabs_mev.txt")
            self.assertEqual(energy_file.name, "tlab_dpol_observables_autogen.txt")
            self.assertEqual(config_file.name, "input_tlab_dpol_observables_autogen.txt")
            self.assertEqual(target_file.read_text(encoding="utf-8"), expected_payload)
            self.assertEqual(energy_file.read_text(encoding="utf-8"), expected_payload)
            self.assertTrue(solver_out_dir.exists())
            self.assertTrue(p123_dir.parent.exists())
            self.assertEqual(format_tlab_dir_name(70.0), "tlab_070MeV")

            config_text = config_file.read_text(encoding="utf-8")
            self.assertIn("energy_input_file=Input/tlab_dpol_observables_autogen.txt", config_text)

    def test_compare_script_emits_tlab_named_summary_schema(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            work_dir = Path(tmpdir) / "work"
            solver_out_dir = Path(tmpdir) / "solver_out"
            solver_out_dir.mkdir(parents=True, exist_ok=True)
            u_file = solver_out_dir / "U_PW_elements_Np_20_Nq_20_JP_1_1_Jmax_2_PSI_0.txt"
            u_file.write_text(_MINIMAL_U_CONTENT, encoding="utf-8")
            (solver_out_dir / "q_kinematics_Nq_20.txt").write_text(_MINIMAL_Q_KINEMATICS, encoding="utf-8")

            cmd = [
                sys.executable,
                str(ROOT / "examples" / "compare_Ay_experiment.py"),
                "--work-dir",
                str(work_dir),
                "--solver-out-dir",
                str(solver_out_dir),
                "--target-tlab-mev",
                "190",
                "--ay-rmse-pass",
                "1",
                "--dsigma-rel-rmse-pass",
                "1000000000",
                "--energy-delta-pass",
                "1000",
            ]
            result = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, msg=f"stdout:\n{result.stdout}\n\nstderr:\n{result.stderr}")

            summary_path = work_dir / "solver_validation_tlab_190MeV.json"
            self.assertTrue(summary_path.exists())

            payload = json.loads(summary_path.read_text(encoding="utf-8"))
            self.assertIn("target_tlab_mev", payload)
            self.assertNotIn("target_tlab", payload)
            self.assertIn("best_energy", payload)
            self.assertIn("tlab_mev", payload["best_energy"])
            self.assertIn("abs_delta_tlab_mev", payload["best_energy"])
            self.assertNotIn("tlab", payload["best_energy"])


if __name__ == "__main__":
    unittest.main()
