"""End-to-end: trace_im_path=true produces a well-formed im_path_trace.txt."""
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
def trace_run(tmp_path_factory):
    out = tmp_path_factory.mktemp("trace_run")
    cfg_text = CFG.read_text()
    cfg_text = re.sub(r"^output_folder=.*$",
                      f"output_folder={out}", cfg_text, flags=re.M)
    # Read-only test: do NOT redirect P123_folder so we reuse the shared cache,
    # but force recovery/store off so the test cannot poison it (the merged-h5
    # integer-truncation bug at make_permutation_matrix.cpp:1310 only triggers
    # when recovery=true + store=true, which we override below).
    cfg_text = re.sub(r"^P123_recovery=.*$", "P123_recovery=false", cfg_text, flags=re.M)
    cfg_text = re.sub(r"^calculate_and_store_P123=.*$",
                      "calculate_and_store_P123=false", cfg_text, flags=re.M)
    if "trace_im_path" not in cfg_text:
        cfg_text += "\ntrace_im_path=true\n"
    else:
        cfg_text = re.sub(r"^trace_im_path=.*$", "trace_im_path=true",
                          cfg_text, flags=re.M)
    cfg_local = out / "input.txt"
    cfg_local.write_text(cfg_text)
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    subprocess.run([str(RUN), str(cfg_local)], check=True, env=env, timeout=900)
    return out

REQUIRED_STAGES = [
    "G0_BC_on_shell_q",
    "G0_CC_straddle_aggregate",
    "K_n0_on_shell_row",
    "AKn_elastic_n0",
    "Pade_best_PA_elastic_U",
    "S_matrix_diagonal_elastic",
]

def test_trace_file_exists(trace_run):
    f = trace_run / "im_path_trace.txt"
    assert f.exists(), f"trace file not produced; dir contents: {list(trace_run.iterdir())}"

def test_trace_has_all_stages(trace_run):
    f = trace_run / "im_path_trace.txt"
    text = f.read_text()
    for stage in REQUIRED_STAGES:
        assert stage in text, f"stage {stage!r} missing in trace; text={text[:500]}"

def test_g0_im_part_nonzero(trace_run):
    """The bound-continuum on-shell row MUST have non-zero ‖Im‖.
    This is the analytic Heaviside guarantee — if it's zero, either the trace
    instrumentation is wrong or the resolvent itself is broken."""
    f = trace_run / "im_path_trace.txt"
    for line in f.read_text().splitlines():
        if line.startswith("G0_BC_on_shell_q"):
            parts = line.split()
            # cols: stage, ‖Re‖, ‖Im‖, Im/Re
            re_norm, im_norm = float(parts[1]), float(parts[2])
            assert im_norm > 1e-6 * max(re_norm, 1.0), \
                f"BC Im part is essentially zero: Re={re_norm}, Im={im_norm}"
            return
    pytest.fail("G0_BC_on_shell_q row not found")
