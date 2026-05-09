# tests/test_pade_honesty.py
"""Regression: U_PW_convergence_*.txt sidecar must appear with valid Conv
codes after the Padé-honesty patch. Cheap 13 MeV Nijmegen-I config, ~3 min wall."""
import os
import re
import subprocess
import sys
from pathlib import Path
import pytest

REPO = Path(__file__).resolve().parents[1]
RUN  = REPO / "CPP" / "run"
CFG  = REPO / "CPP" / "Input" / "input_miller_gate1_dbg.txt"

@pytest.fixture(scope="module")
def cheap_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("pade_honesty_run")
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    cfg_text = CFG.read_text()
    cfg_text = re.sub(r"^output_folder=.*$",
                      f"output_folder={out}", cfg_text, flags=re.M)
    cfg_local = out / "input.txt"
    cfg_local.write_text(cfg_text)
    subprocess.run([str(RUN), str(cfg_local)], check=True, env=env, timeout=900)
    return out

def test_conv_column_present(cheap_run):
    sidecars = list(cheap_run.glob("U_PW_convergence_*.txt"))
    assert sidecars, "no U_PW_convergence_*.txt sidecar produced"
    lines = sidecars[0].read_text().strip().splitlines()
    data_rows = [l for l in lines if not l.startswith("#")]
    assert data_rows, "sidecar has no data rows"
    cols = data_rows[0].split()
    assert len(cols) == 5, f"expected 5 columns, got {len(cols)}: {cols}"
    conv_codes = {int(l.split()[3]) for l in data_rows}
    # at least one row should be coded; both 1 and 2 may appear
    assert conv_codes.issubset({0, 1, 2}) and conv_codes != {0}, \
        f"conv codes look wrong: {conv_codes}"
