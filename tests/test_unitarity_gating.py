"""extract_phase_shifts.py with --markdown-out must mark blocks where
||SS†-1|| > 0.2 as FAIL in the markdown report."""
import subprocess
import sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[1]
EXTRACT = REPO / "examples" / "extract_phase_shifts.py"


@pytest.fixture
def synthetic_solver_out(tmp_path):
    """Hand-craft a U_PW_elements file whose S has unit_def well above the
    0.2 threshold. Use a 1x1 channel block with U deliberately chosen so that
    |S| >> 1 and ||SS†-1|| is large.
    """
    out = tmp_path / "solver_out"
    out.mkdir()
    # NOTE: Real format is verbose. The test should construct the smallest
    # valid file that extract_phase_shifts.py accepts — match the header
    # written by store_U_matrix_elements_txt.
    u_path = out / "U_PW_elements_Nq_2_Np_2_two_J_3N_1_P_3N_+1.txt"
    u_path.write_text(
        "# Elastic Nd-scattering U-matrix elements ...\n"
        "# JP:          1/2+\n"
        "# Np:          2\n"
        "# Nq:          2\n"
        "# Particles:   nd-scattering\n"
        "#\n"
        "# (header truncated for synthetic test)\n"
        "#\n"
        "# ###################################################################################################\n"
        "#      Name   row-idx   col-idx       l'      2*j'        l       2*j\n"
        "        U00         0         0         0         1         0         1\n"
        "# ###################################################################################################\n"
        "#              Tlab [MeV]               Ecm [MeV]     q_idx                                         U00 [MeV]\n"
        "   1.2000000000000000e+01  8.0000000000000000e+00         0   +1.5000000000000000e+03+0.0000000000000000e+00j\n"
        "# ###################################################################################################\n"
    )
    q_kin = out / "q_kinematics_Nq_2.txt"
    q_kin.write_text(
        "# BOUNDARIES:\n"
        "# Index                  q [MeV]                Ecm [MeV]               Tlab [MeV]\n"
        "      0          +0.00000000e+00          +0.00000000e+00          +0.00000000e+00\n"
        "      1          +5.00000000e+01          +2.00000000e+00          +3.00000000e+00\n"
        "      2          +1.00000000e+02          +8.00000000e+00          +1.20000000e+01\n"
        "# MID-POINTS:\n"
        "# Index                  q [MeV]                Ecm [MeV]               Tlab [MeV]\n"
        "      0          +2.50000000e+01          +1.00000000e+00          +1.50000000e+00\n"
        "      1          +7.50000000e+01          +5.00000000e+00          +7.50000000e+00\n"
    )
    return out


def test_unit_def_gating_marks_fail(synthetic_solver_out, tmp_path):
    md_out = tmp_path / "report.md"
    # Tolerate the script erroring on synthetic data — we only need to
    # verify it produces a markdown file with FAIL marker for high unit_def.
    result = subprocess.run(
        [sys.executable, str(EXTRACT),
         "--solver-out-dir", str(synthetic_solver_out),
         "--markdown-out", str(md_out),
         "--unit-def-threshold", "0.2"],
        capture_output=True, text=True
    )
    assert md_out.exists(), f"markdown not written: {result.stderr}"
    text = md_out.read_text()
    assert "FAIL" in text, f"expected FAIL marker in report; got:\n{text}"
