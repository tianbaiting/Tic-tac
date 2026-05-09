# tests/test_pade_honesty.py
"""Regression: Conv column must appear in U_PW_elements_*.txt after the
Padé-honesty patch. Cheap 13 MeV Nijmegen-I config, ~3 min wall."""
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
    files = list(cheap_run.glob("U_PW_elements_*.txt"))
    assert files, "no U_PW_elements_*.txt produced"
    text = files[0].read_text()
    assert " Conv " in text or "\tConv\t" in text, \
        f"Conv column missing in {files[0].name}"
