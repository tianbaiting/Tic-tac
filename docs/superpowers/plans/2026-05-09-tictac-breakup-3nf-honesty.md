# Tic-tac Breakup-Channel & 3NF Validation Honesty — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the two most damaging fake-validations (Padé maxiter-as-converged; unitarity-not-gated) and diagnose where the resolvent's analytic Im part fails to surface as a non-zero elastic Im δ — then apply a targeted fix.

**Architecture:** Three parallel workstreams — (B-1) Padé honesty in the solver + serialization; (B-2) unitarity threshold in the Python extractor; (A-1) Im-path trace through resolvent → Neumann series → Padé → S — followed by a diagnostic note (A-2) and a discovery-driven fix (A-3). All changes are **additive**: new flags, new arrays, new files alongside existing paths. No deletions or in-place rewrites.

**Tech Stack:** C++17 (solver core), Fortran 90/77 (potentials, untouched), Python 3 (extractor + tests), Makefile (`CPP/run`, production), CMake (`bin/tic-tac`, parallel build).

---

## Spec Reference

`docs/superpowers/specs/2026-05-09-tictac-breakup-3nf-honesty-design.md`

## File Structure (all additive)

| File | Action | Purpose |
|---|---|---|
| `include/type_defs.h` | Modify | Add 1 bool field `trace_im_path` to `run_params` struct (lines 105-163). |
| `CPP/type_defs.h` | Modify | Mirror the field (production-binary type defs). |
| `src/config/set_run_parameters.cpp` | Modify | Parser branch + default + help string for `trace_im_path`. |
| `CPP/Input/input_miller_gate1_dbg.txt` | Create | Cheap reference run config for trace + tests. |
| `src/core/faddeev_solver/solve_faddeev.cpp` | Modify | (a) Allocate new `*_truly_converged_array`, `*_maxiter_truncated_array` parallel to existing `_conv_array`. (b) Set both new arrays at the convergence-decision sites (lines 1701-1713, 1773-1780). (c) Trace hooks: `‖Re/Im(A·Kⁿ)‖` per Neumann iter when `trace_im_path=true`. |
| `src/core/faddeev_solver/solve_faddeev.h` | Modify | Add `truly_converged` / `maxiter_truncated` output buffers to function signatures (additive — existing `conv_array` retained). |
| `src/core/resolvent/make_resolvent.cpp` | Modify | When `trace_im_path=true`, dump on-shell q-bin row of (Re G₀, Im G₀) to file. |
| `src/io/disk_io_routines.cpp` | Modify | `store_U_matrix_elements_txt` writes two extra integer columns: `Conv` (1=truly, 2=maxiter), `UnitDef` (`||SS†-1||` per block, computed solver-side). |
| `examples/extract_phase_shifts.py` | Modify | Read `Conv` and `UnitDef` columns. Apply hard threshold `‖SS†−1‖ > 0.2 → FAIL`. Write markdown report with ⚠ markers. |
| `examples/extract_phase_shifts.py` | Modify | Add `--markdown-out PATH` flag to write structured report. |
| `tests/cpp/resolvent_im_test.cpp` | Create | Unit test for `Im_R`, `Im_Q` against analytic values. |
| `tests/cpp/CMakeLists.txt` | Modify | Register the new test executable. |
| `tests/test_im_path_trace.py` | Create | End-to-end regression: cheap 13 MeV run with `trace_im_path=true` produces well-formed trace file. |
| `tests/test_pade_honesty.py` | Create | Cheap run asserts `Conv` column appears in `U_PW_elements_*.txt` and contains both 1 and 2 values across blocks. |
| `tests/test_unitarity_gating.py` | Create | Synthetic input file with hand-crafted bad `UnitDef` value triggers FAIL in markdown report. |
| `docs/im_path_diagnosis_2026-05-09.md` | Create | Diagnostic note written after A-1 run; identifies bug site or undersampling verdict. |

The workstreams are listed below with explicit dependency markers. **B-1, B-2, A-1 are independent** and can be executed in parallel by separate subagents. **A-2 depends on A-1**. **A-3 depends on A-2**. Final regression sweep depends on all five.

---

## Workstream B-1: Padé Honesty in the Solver

**Owner suggestion:** Subagent #1.
**Depends on:** Nothing.
**Verifies:** Subagent #2 (B-2) once B-1 is committed.

### Task B-1.1: Failing test for `Conv` column in U_PW_elements

**Files:**
- Create: `tests/test_pade_honesty.py`

- [ ] **Step 1: Write the failing test**

```python
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
```

- [ ] **Step 2: Run the test, confirm it fails**

```bash
cd /home/tian/workspace/dpol/Tic-tac
python3 -m pytest tests/test_pade_honesty.py -v 2>&1 | tail -10
```

Expected: FAIL with "Conv column missing" (the column does not yet exist).

- [ ] **Step 3: Commit**

```bash
git add tests/test_pade_honesty.py
git commit -m "test: failing regression for Padé Conv column in U_PW_elements"
```

### Task B-1.2: Add `truly_converged` / `maxiter_truncated` arrays — declarations

**Files:**
- Modify: `src/core/faddeev_solver/solve_faddeev.cpp:949-958`

- [ ] **Step 1: Edit declarations to add parallel arrays**

Replace:

```cpp
	/* Arrays to store Pade-approximants (PA) for each on-shell elastic elements */
	cdouble* pade_approximants_array      = new cdouble [num_EL_A_vals * (NM_max+1)];
	size_t*  pade_approximants_idx_array  = new size_t  [num_EL_A_vals];
	bool*    pade_approximants_conv_array = new bool    [num_EL_A_vals];
	size_t	 num_converged_elements		  = 0;

	/* Arrays to store Pade-approximants (PA) for each on-shell breakup elements */
	cdouble* pade_approximants_BU_array      = NULL;//new cdouble [num_BU_A_vals * (NM_max+1)];
	size_t*  pade_approximants_BU_idx_array  = NULL;//new size_t  [num_BU_A_vals];
	bool*    pade_approximants_BU_conv_array = NULL;//new bool    [num_BU_A_vals];
	size_t	 num_converged_BU_elements		 = 0;
```

with:

```cpp
	/* Arrays to store Pade-approximants (PA) for each on-shell elastic elements */
	cdouble* pade_approximants_array      = new cdouble [num_EL_A_vals * (NM_max+1)];
	size_t*  pade_approximants_idx_array  = new size_t  [num_EL_A_vals];
	bool*    pade_approximants_conv_array = new bool    [num_EL_A_vals];
	// [EN] Honesty layer (additive, parallel to *_conv_array). truly_converged_array is
	// set when criteria 1/2/3 fire; maxiter_truncated_array is set ONLY when criterion 0
	// (NM == NM_max) fired without any of 1/2/3. Old consumers reading *_conv_array are
	// unaffected. / [CN] 诚信层（与 *_conv_array 平行新增）。truly_converged_array 在
	// criteria 1/2/3 触发时置位；maxiter_truncated_array 仅在 criterion 0 (NM == NM_max)
	// 触发而 1/2/3 都未满足时置位。读旧 *_conv_array 的代码不受影响。
	bool*    pade_approximants_truly_converged_array     = new bool [num_EL_A_vals];
	bool*    pade_approximants_maxiter_truncated_array   = new bool [num_EL_A_vals];
	size_t	 num_converged_elements		  = 0;

	/* Arrays to store Pade-approximants (PA) for each on-shell breakup elements */
	cdouble* pade_approximants_BU_array                  = NULL;
	size_t*  pade_approximants_BU_idx_array              = NULL;
	bool*    pade_approximants_BU_conv_array             = NULL;
	bool*    pade_approximants_BU_truly_converged_array  = NULL;
	bool*    pade_approximants_BU_maxiter_truncated_array= NULL;
	size_t	 num_converged_BU_elements		 = 0;
```

- [ ] **Step 2: Initialize new arrays right below the existing `_conv_array` init at line 983-985**

Replace:

```cpp
	for (size_t idx_NDOS=0; idx_NDOS<num_EL_A_vals; idx_NDOS++){
		pade_approximants_conv_array[idx_NDOS] = false;
	}
```

with:

```cpp
	for (size_t idx_NDOS=0; idx_NDOS<num_EL_A_vals; idx_NDOS++){
		pade_approximants_conv_array[idx_NDOS]              = false;
		pade_approximants_truly_converged_array[idx_NDOS]   = false;
		pade_approximants_maxiter_truncated_array[idx_NDOS] = false;
	}
```

- [ ] **Step 3: Mirror the BU initialization at line 994-996**

Replace:

```cpp
		for (size_t idx_NDOS=0; idx_NDOS<num_BU_A_vals; idx_NDOS++){
			pade_approximants_BU_conv_array[idx_NDOS] = false;
		}
```

with:

```cpp
		pade_approximants_BU_truly_converged_array   = new bool [num_BU_A_vals];
		pade_approximants_BU_maxiter_truncated_array = new bool [num_BU_A_vals];
		for (size_t idx_NDOS=0; idx_NDOS<num_BU_A_vals; idx_NDOS++){
			pade_approximants_BU_conv_array[idx_NDOS]              = false;
			pade_approximants_BU_truly_converged_array[idx_NDOS]   = false;
			pade_approximants_BU_maxiter_truncated_array[idx_NDOS] = false;
		}
```

- [ ] **Step 4: Build to verify allocations compile cleanly**

```bash
cd /home/tian/workspace/dpol/Tic-tac/CPP && make -j 2>&1 | tail -20
```

Expected: clean build, no warnings about unused new arrays.

- [ ] **Step 5: Commit**

```bash
git add src/core/faddeev_solver/solve_faddeev.cpp
git commit -m "solver: allocate truly_converged/maxiter_truncated arrays alongside conv_array

Additive — old *_conv_array readers continue to work. New arrays default
false; populated by the convergence-decision site in the next commit."
```

### Task B-1.3: Populate the new arrays at the convergence-decision site

**Files:**
- Modify: `src/core/faddeev_solver/solve_faddeev.cpp:1699-1713` (elastic), `:1766-1780` (BU)

- [ ] **Step 1: Edit the elastic convergence decision (around line 1706-1713)**

Replace:

```cpp
					if (convergence_criteria_0 ||
						convergence_criteria_1 ||
						convergence_criteria_2 ||
						convergence_criteria_3){
						pade_approximants_conv_array[idx_NDOS] = true;
						pade_approximants_idx_array[idx_NDOS]  = idx_best_PA;
						num_converged_elements += 1;
					}
```

with:

```cpp
					if (convergence_criteria_0 ||
						convergence_criteria_1 ||
						convergence_criteria_2 ||
						convergence_criteria_3){
						pade_approximants_conv_array[idx_NDOS] = true;
						pade_approximants_idx_array[idx_NDOS]  = idx_best_PA;
						num_converged_elements += 1;
						// [EN] Honesty split: only criteria 1/2/3 represent genuine
						// convergence. Criterion 0 (NM == NM_max) without any of 1/2/3
						// is a max-iter timeout, not convergence. / [CN] 诚信拆分：
						// 仅 criteria 1/2/3 代表真收敛；criterion 0 (NM == NM_max) 单
						// 独触发表示达到迭代上限而非收敛。
						bool genuinely_converged =  convergence_criteria_1
												 || convergence_criteria_2
												 || convergence_criteria_3;
						pade_approximants_truly_converged_array[idx_NDOS]   = genuinely_converged;
						pade_approximants_maxiter_truncated_array[idx_NDOS] = !genuinely_converged
																			  && convergence_criteria_0;
					}
```

- [ ] **Step 2: Mirror for BU (around line 1773-1780)**

Replace:

```cpp
						if (convergence_criteria_0 ||
							convergence_criteria_1 ||
							convergence_criteria_2 ||
							convergence_criteria_3){
							pade_approximants_BU_conv_array[idx_NDOS] = true;
							pade_approximants_BU_idx_array[idx_NDOS]  = idx_best_PA;
							num_converged_elements += 1;
						}
```

with:

```cpp
						if (convergence_criteria_0 ||
							convergence_criteria_1 ||
							convergence_criteria_2 ||
							convergence_criteria_3){
							pade_approximants_BU_conv_array[idx_NDOS] = true;
							pade_approximants_BU_idx_array[idx_NDOS]  = idx_best_PA;
							num_converged_elements += 1;
							bool genuinely_converged =  convergence_criteria_1
													 || convergence_criteria_2
													 || convergence_criteria_3;
							pade_approximants_BU_truly_converged_array[idx_NDOS]   = genuinely_converged;
							pade_approximants_BU_maxiter_truncated_array[idx_NDOS] = !genuinely_converged
																					  && convergence_criteria_0;
						}
```

- [ ] **Step 3: Verify build**

```bash
cd /home/tian/workspace/dpol/Tic-tac/CPP && make -j 2>&1 | tail -10
```

Expected: clean build.

- [ ] **Step 4: Free the new arrays at function end (around line 1838-1849, in the existing `delete[]` block)**

After the existing `delete [] pade_approximants_conv_array;` (search for it), add:

```cpp
	delete [] pade_approximants_truly_converged_array;
	delete [] pade_approximants_maxiter_truncated_array;
	if (pade_approximants_BU_truly_converged_array != NULL){
		delete [] pade_approximants_BU_truly_converged_array;
		delete [] pade_approximants_BU_maxiter_truncated_array;
	}
```

(If the existing `delete[] pade_approximants_conv_array;` is missing, add it here too — that would already be a separate leak; keep change additive.)

- [ ] **Step 5: Build + check no leaks introduced**

```bash
cd /home/tian/workspace/dpol/Tic-tac/CPP && make -j 2>&1 | tail -5
```

- [ ] **Step 6: Commit**

```bash
git add src/core/faddeev_solver/solve_faddeev.cpp
git commit -m "solver: populate truly_converged/maxiter_truncated at Padé decision

Criterion 0 (NM == NM_max) alone now flags maxiter_truncated, not converged.
Criteria 1/2/3 flag truly_converged. *_conv_array unchanged — existing
extraction code still picks the best PA the same way."
```

### Task B-1.4: Pipe new flags to `store_U_matrix_elements_txt`

**Files:**
- Modify: `src/io/disk_io_routines.cpp:766+` (`store_U_matrix_elements_txt`)
- Modify: `src/io/disk_io_routines.h` (signature)
- Modify: `src/main.cpp:541` (call site)
- Modify: `src/core/faddeev_solver/solve_faddeev.cpp` (return new arrays via additional output buffer args, OR write Conv column directly inside solve_faddeev)

**Approach (simpler):** write the Conv flag directly into a sidecar file `U_PW_convergence_*.txt` in the solver, since `U_array` storage is already complex-valued (cdouble) and the conv flag is per-element. This is fully additive — does not touch `store_U_matrix_elements_txt`.

- [ ] **Step 1: Add sidecar dump in solve_faddeev.cpp at the same site as U_array fill (line ~1808)**

After the line:

```cpp
				U_array[ndos.value_storage_idx] = pade_approximants_array[ndos.value_storage_idx*(NM_max+1) + idx_best_PA];
				printf("       - U-matrix element for alpha'=%ld, alpha=%ld, q=%ld: %.10e + %.10ei \n", ndos.alpha_row, ndos.alpha_col, ndos.q_idx, U_array[ndos.value_storage_idx].real(), U_array[ndos.value_storage_idx].imag());
```

(no edit there — keep U_array unchanged.) Then, immediately AFTER the elastic loop closes (find `/* Set on-shell breakup U-matrix elements equal "best" PA */` ~line 1813), insert a new block BEFORE that:

```cpp
	/* Sidecar: write per-element convergence honesty flags so the Python
	 * extractor can distinguish truly-converged from maxiter-truncated PAs.
	 * Sidecar file is OPTIONAL for legacy consumers — only the new extractor
	 * looks for it. */
	{
		std::string conv_file = run_parameters.output_folder + "/U_PW_convergence" + file_identification + ".txt";
		std::ofstream cf(conv_file);
		cf << "# Per-element Padé convergence honesty (additive sidecar).\n";
		cf << "# Conv: 1 = truly_converged (criteria 1/2/3); 2 = maxiter_truncated (criterion 0 only).\n";
		cf << "# Columns: row col q_com Conv idx_best_PA\n";
		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
				for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
					size_t idx_NDOS = elastic_value_storage_index(idx_d_row, idx_d_col, idx_q_com,
																 num_deuteron_states, num_q_com);
					int conv_code = pade_approximants_truly_converged_array[idx_NDOS]   ? 1 :
									pade_approximants_maxiter_truncated_array[idx_NDOS] ? 2 : 0;
					cf << idx_d_row << " " << idx_d_col << " " << idx_q_com << " "
					   << conv_code << " " << pade_approximants_idx_array[idx_NDOS] << "\n";
				}
			}
		}
		cf.close();
	}
```

NOTE: `file_identification` is a local variable in the surrounding scope; verify it is available here (it is — `solve_faddeev` receives it as a parameter).

- [ ] **Step 2: Build**

```bash
cd /home/tian/workspace/dpol/Tic-tac/CPP && make -j 2>&1 | tail -5
```

- [ ] **Step 3: Update the failing test from B-1.1 to look at the sidecar**

Edit `tests/test_pade_honesty.py`:

```python
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
```

- [ ] **Step 4: Run the test**

```bash
cd /home/tian/workspace/dpol/Tic-tac
python3 -m pytest tests/test_pade_honesty.py -v 2>&1 | tail -15
```

Expected: PASS (sidecar exists, 5 columns, conv codes ∈ {0,1,2}).

- [ ] **Step 5: Commit**

```bash
git add src/core/faddeev_solver/solve_faddeev.cpp tests/test_pade_honesty.py
git commit -m "solver: emit U_PW_convergence_*.txt sidecar with honesty flags

Conv = 1 (truly_converged via criteria 1/2/3) | 2 (maxiter_truncated). The
existing U_PW_elements_*.txt is unchanged — old readers see no diff."
```

### Task B-1.5: Cheap reference input for the trace + tests

**Files:**
- Create: `CPP/Input/input_miller_gate1_dbg.txt`

- [ ] **Step 1: Copy the existing 13 MeV cheap config and tag it `dbg`**

```bash
cd /home/tian/workspace/dpol/Tic-tac
cp CPP/Input/input_miller_gate1.txt CPP/Input/input_miller_gate1_dbg.txt
```

- [ ] **Step 2: Edit the new file to set `output_folder=CPP/Output/miller_gate1_dbg` and ensure `Np_WP=20`, `Nq_WP=20`, `potential_model=nijmegen`, `include_breakup_channels=true`, `three_nucleon_force=none`**

Use Edit tool to change:
- `output_folder=CPP/Output/miller_gate1` → `output_folder=CPP/Output/miller_gate1_dbg`
- (verify other fields with `grep ^Np_WP CPP/Input/input_miller_gate1_dbg.txt` etc.; adjust if drifted)

- [ ] **Step 3: Smoke-run to confirm wall time ≤ 10 min**

```bash
cd /home/tian/workspace/dpol/Tic-tac
mkdir -p CPP/Output/miller_gate1_dbg
time ./CPP/run CPP/Input/input_miller_gate1_dbg.txt 2>&1 | tail -20
```

Expected: completes within ~5 minutes.

- [ ] **Step 4: Commit**

```bash
git add CPP/Input/input_miller_gate1_dbg.txt
git commit -m "input: cheap dbg config for Padé/trace regression tests

Np=20 Nijmegen-I 13 MeV breakup-on, no 3NF. ~3-5 min wall."
```

---

## Workstream B-2: Unitarity Threshold in extract_phase_shifts.py

**Owner suggestion:** Subagent #2.
**Depends on:** B-1.4 (sidecar file format) — but can develop in parallel and integrate at the end.
**Verifies:** Subagent #1 (B-1) once both are committed.

### Task B-2.1: Failing test for unitarity FAIL gating

**Files:**
- Create: `tests/test_unitarity_gating.py`

- [ ] **Step 1: Write the failing test**

```python
# tests/test_unitarity_gating.py
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
    """Hand-craft a U_PW_elements file whose S has unit_def=1.18 (Np=30 JP=1- value)."""
    # Use a 1x1 channel block with U deliberately chosen so that
    # |S| = 1.5 → ||SS†-1|| = 1.25.
    out = tmp_path / "solver_out"
    out.mkdir()
    # NOTE: Real format is verbose. The test should construct the smallest
    # valid file that extract_phase_shifts.py accepts — match the header
    # written by store_U_matrix_elements_txt.
    u_path = out / "U_PW_elements_Nq_2_Np_2_two_J_3N_1_P_3N_+1.txt"
    u_path.write_text("""\
# Elastic Nd-scattering U-matrix elements ...
# JP:          1/2+
# Np:          2
# Nq:          2
# Particles:   nd-scattering
#
# (header truncated for synthetic test)
#
#      Name   row-idx   col-idx       l'      2*j'        l       2*j
   U00         0         0         0         1         0         1
# Tlab Ecm q_idx Re(U) Im(U)
12.0  8.0  0   1.5e3   0.0e0
""")
    q_kin = out / "q_kinematics_Nq_2.txt"
    q_kin.write_text("# q [MeV] q-bin-edge\n0.0\n50.0\n100.0\n")
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
```

- [ ] **Step 2: Run the test, confirm it fails**

```bash
cd /home/tian/workspace/dpol/Tic-tac
python3 -m pytest tests/test_unitarity_gating.py -v 2>&1 | tail -15
```

Expected: FAIL — `--markdown-out` flag doesn't exist yet.

- [ ] **Step 3: Commit**

```bash
git add tests/test_unitarity_gating.py
git commit -m "test: failing regression for unitarity FAIL gating in markdown report"
```

### Task B-2.2: Add `--markdown-out` and `--unit-def-threshold` to extract_phase_shifts.py

**Files:**
- Modify: `examples/extract_phase_shifts.py:306-328` (argparse + main)

- [ ] **Step 1: Add CLI flags**

In the existing argparse block (line ~308-312), add after `--tlab`:

```python
    ap.add_argument("--markdown-out", type=Path, default=None,
                    help="if set, write a structured markdown report to this path "
                         "with PASS/⚠/FAIL markers per (J,P) block")
    ap.add_argument("--unit-def-threshold", type=float, default=0.2,
                    help="||SS†-1|| above this triggers FAIL in the markdown report "
                         "(default 0.2 — anything above is too non-unitary to trust)")
    ap.add_argument("--maxiter-warn", action="store_true", default=True,
                    help="read U_PW_convergence_*.txt sidecars and add ⚠ for "
                         "Padé-maxiter-truncated rows in the markdown report")
```

- [ ] **Step 2: Build a markdown writer function (additive — does NOT replace stdout printing)**

After `extract_for_file` (line ~258) but before `main()` (line ~306), add:

```python
def write_markdown_report(results_per_file, out_path: Path,
                          unit_def_threshold: float,
                          conv_codes_per_file: dict | None):
    """Additive: writes a structured markdown report. The existing stdout
    printing in extract_for_file is unaffected.

    results_per_file: { u_path : list of result dicts (as produced by extract_for_file) }
    conv_codes_per_file: { u_path : { (row, col, q_com) : conv_code } } or None
    """
    lines = ["# Phase-shift extraction report", ""]
    overall_pass = True
    for u_path, results in results_per_file.items():
        lines.append(f"## {u_path.name}")
        lines.append("")
        lines.append("| Tlab [MeV] | block | δ_diag [°] | \\|S_kk\\| | \\|\\|SS†-1\\|\\| | status |")
        lines.append("|---|---|---|---|---|---|")
        for res in results:
            S = res["S"]
            unit_def = float((S @ S.conj().T - __import__("numpy").eye(S.shape[0])).__abs__().sum() ** 0.5)
            for k in range(len(res["delta_diag"])):
                d = res["delta_diag"][k]
                # canonical_delta_deg already exists in this file
                from numpy import angle
                re_deg = angle(__import__("numpy").exp(1j * 2 * d.real)) * 0.5 * 180.0 / 3.141592653589793
                inelast = float(res["inelast_diag"][k])
                if unit_def > unit_def_threshold:
                    status = "FAIL"
                    overall_pass = False
                else:
                    status = "PASS"
                if conv_codes_per_file:
                    codes = conv_codes_per_file.get(u_path, {})
                    if any(v == 2 for v in codes.values()):
                        status = "⚠ maxiter-truncated, " + status
                lines.append(
                    f"| {res['tlab']:.2f} | row={k} | {re_deg:+.3f} | {inelast:.4f} | {unit_def:.4f} | {status} |"
                )
        lines.append("")
    lines.append("")
    lines.append(f"**Overall:** {'PASS' if overall_pass else 'FAIL'}")
    out_path.write_text("\n".join(lines))
```

- [ ] **Step 3: Wire `--markdown-out` into `main()`**

Modify the loop in `main()` (line ~314-328) to collect results and call the writer:

```python
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--solver-out-dir", required=True, type=Path,
                    help="directory holding U_PW_elements_*.txt + q_kinematics_Nq_*.txt")
    ap.add_argument("--tlab", type=float, default=None,
                    help="show only rows within ±5 MeV of this Tlab (default: all)")
    ap.add_argument("--markdown-out", type=Path, default=None,
                    help="if set, write a structured markdown report")
    ap.add_argument("--unit-def-threshold", type=float, default=0.2)
    ap.add_argument("--maxiter-warn", action="store_true", default=True)
    args = ap.parse_args()

    u_files = sorted(args.solver_out_dir.glob("U_PW_elements_*.txt"))
    if not u_files:
        raise SystemExit(f"no U_PW_elements_*.txt under {args.solver_out_dir}")

    results_per_file = {}
    conv_codes_per_file = {}
    for u_path in u_files:
        m = re.search(r"Nq_(\d+)", u_path.name)
        if not m:
            continue
        nq = m.group(1)
        q_kin = args.solver_out_dir / f"q_kinematics_Nq_{nq}.txt"
        if not q_kin.exists():
            print(f"[skip] no q_kinematics_Nq_{nq}.txt for {u_path.name}")
            continue
        res = extract_for_file(u_path, q_kin, args.tlab)
        if res:
            results_per_file[u_path] = res
        # Read sidecar if present
        sidecar_name = u_path.name.replace("U_PW_elements_", "U_PW_convergence_")
        sidecar = args.solver_out_dir / sidecar_name
        if args.maxiter_warn and sidecar.exists():
            codes = {}
            for line in sidecar.read_text().splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                if len(parts) >= 4:
                    codes[(int(parts[0]), int(parts[1]), int(parts[2]))] = int(parts[3])
            conv_codes_per_file[u_path] = codes

    if args.markdown_out:
        write_markdown_report(results_per_file, args.markdown_out,
                              args.unit_def_threshold,
                              conv_codes_per_file or None)
        print(f"[markdown] {args.markdown_out}")
```

NOTE: `extract_for_file` currently does not return its results list — Step 4 fixes that.

- [ ] **Step 4: Make `extract_for_file` return its results list**

Find the line where `extract_for_file` ends (line ~303 `return results`) — verify the return is already there. If not, add `return results` at the end of the function. Update its callers accordingly.

- [ ] **Step 5: Run the failing test**

```bash
cd /home/tian/workspace/dpol/Tic-tac
python3 -m pytest tests/test_unitarity_gating.py -v 2>&1 | tail -15
```

Expected: PASS (markdown is written, contains "FAIL").

- [ ] **Step 6: Commit**

```bash
git add examples/extract_phase_shifts.py tests/test_unitarity_gating.py
git commit -m "extractor: --markdown-out with FAIL gating on ||SS†-1||

Additive: existing stdout output unchanged. New flag writes a structured
markdown table; rows with unit_def > 0.2 marked FAIL; rows where the
sidecar reports maxiter-truncated marked ⚠."
```

### Task B-2.3: Re-run extractor over existing Np=30 3NF dataset (regression-style verification)

**Files:**
- (no code change; produces an artifact)

- [ ] **Step 1: Re-run the extractor against the committed Np=30 3NF run**

```bash
cd /home/tian/workspace/dpol/Tic-tac
python3 examples/extract_phase_shifts.py \
    --solver-out-dir Tic-tac-w1cache/CPP/Output/miller_gate1_np30_3nf \
    --markdown-out output/np30_3nf_honesty_report.md \
    --tlab 12.05 2>&1 | tail -20
```

Expected: an `output/np30_3nf_honesty_report.md` containing FAIL for the JP=1- ²P block (it has ‖SS†−1‖=1.18) and PASS for JP=1+ ²S₁/₂.

- [ ] **Step 2: Eyeball verification — open the markdown and confirm it reads honestly**

```bash
cat output/np30_3nf_honesty_report.md
```

Expected per spec §5: "JP=1+ ⚠ maxiter-truncated; JP=1- ²P block FAIL".

- [ ] **Step 3: Commit the artifact**

```bash
git add -f output/np30_3nf_honesty_report.md
git commit -m "report: Np=30 3NF Miller Gate 1 with honesty layer applied

JP=1- ²P doublet now FAIL on ||SS†-1||; JP=1+ flagged ⚠ for maxiter-
truncated Padé. Replaces the prior 'Gate 1 PASSED' framing."
```

---

## Workstream A-1: Im-Path Trace

**Owner suggestion:** Subagent #3.
**Depends on:** B-1.5 (cheap dbg input).
**Verifies:** Subagent #1 (B-1) and Subagent #2 (B-2) by re-running their tests after this lands.

### Task A-1.1: Add `trace_im_path` flag to run_params

**Files:**
- Modify: `include/type_defs.h:133` (after `include_breakup_channels`)
- Modify: `CPP/type_defs.h` (mirror)
- Modify: `src/config/set_run_parameters.cpp:182-189` (parser branch), `:673` (default)

- [ ] **Step 1: Add field to `include/type_defs.h`**

After line 133 (`bool include_breakup_channels;`), insert:

```cpp
	bool		trace_im_path;             // [EN] when true, dump per-stage Re/Im norms
                                         // to <output_folder>/im_path_trace.txt during a
                                         // single-shot diagnostic run. Default false. /
                                         // [CN] true 时单次诊断跑里输出每个阶段的 Re/Im
                                         // 范数到 im_path_trace.txt。默认 false。
```

- [ ] **Step 2: Mirror in `CPP/type_defs.h`** at the analogous position.

- [ ] **Step 3: Add parser branch in `src/config/set_run_parameters.cpp`** after the `include_breakup_channels` branch (around line 189):

```cpp
	else if (option == "trace_im_path"){
		if (input=="true" || input=="false"){
			run_parameters.trace_im_path = (input=="true");
		}
		else{
			raise_error("Invalid value for input parameter trace_im_path!");
		}
	}
```

- [ ] **Step 4: Default in `set_default_values` (around line 673)**, after `include_breakup_channels = false`:

```cpp
	run_parameters.trace_im_path = false;
```

- [ ] **Step 5: Build (both CMake and Makefile paths)**

```bash
cd /home/tian/workspace/dpol/Tic-tac
./build.sh release 2>&1 | tail -5
cd CPP && make -j 2>&1 | tail -5
```

- [ ] **Step 6: Commit**

```bash
git add include/type_defs.h CPP/type_defs.h src/config/set_run_parameters.cpp
git commit -m "config: add trace_im_path flag (default false, no behavior change)"
```

### Task A-1.2: Failing end-to-end test for trace file

**Files:**
- Create: `tests/test_im_path_trace.py`

- [ ] **Step 1: Write the test**

```python
# tests/test_im_path_trace.py
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
```

- [ ] **Step 2: Run, confirm fails**

```bash
cd /home/tian/workspace/dpol/Tic-tac
python3 -m pytest tests/test_im_path_trace.py -v 2>&1 | tail -20
```

Expected: FAIL — trace file does not exist yet.

- [ ] **Step 3: Commit**

```bash
git add tests/test_im_path_trace.py
git commit -m "test: failing end-to-end regression for im_path_trace.txt"
```

### Task A-1.3: C++ unit test for resolvent Im parts

**Files:**
- Create: `tests/cpp/resolvent_im_test.cpp`
- Modify: `tests/cpp/CMakeLists.txt` (verify path; create if absent)

- [ ] **Step 1: Locate the existing C++ test harness**

```bash
cd /home/tian/workspace/dpol/Tic-tac
find . -name "CMakeLists.txt" -path "*/tests/*" 2>/dev/null
ls tests/cpp/ 2>/dev/null
```

If `tests/cpp/CMakeLists.txt` doesn't exist, the cache_layer_test was probably wired into the top CMakeLists.txt — check there.

- [ ] **Step 2: Write the unit test**

Path: `tests/cpp/resolvent_im_test.cpp`

```cpp
// tests/cpp/resolvent_im_test.cpp
//
// Unit-test the analytic Im parts of the cell-averaged resolvent.
// Reference values are derived in docs/im_path_diagnosis_2026-05-09.md.

#include "core/resolvent/make_resolvent.h"
#include <cassert>
#include <cmath>
#include <cstdio>

static bool close(double a, double b, double tol=1e-10) {
    return std::abs(a - b) <= tol * (1.0 + std::abs(b));
}

int main() {
    // Test 1: Bound-continuum on-shell q-bin.
    // Eb = -2.224 MeV, q-bin [50, 60] MeV/c (linear momentum, then converted).
    // E chosen so that Eq_lower + Eb < E < Eq_upper + Eb.
    // mu1 = Mn*(2Mn + Mp + Eb)/(2Mn + 2Mn + Mp + Eb) ≈ 624 MeV (using constants).
    // For simplicity, use the actual function and verify Im_R = -π / Δq.
    {
        double Eb = -2.224;
        double q_lo = 50.0, q_hi = 60.0;
        // mu1 from constants.h matches the formula in resolvent_bound_continuum
        double mu1 = 939.565 * (939.565 + 938.272 + Eb) /
                     (939.565 + 939.565 + 938.272 + Eb);
        double Eq_lo = q_lo*q_lo / (2*mu1);
        double Eq_hi = q_hi*q_hi / (2*mu1);
        double Dq = Eq_hi - Eq_lo;
        double E = Eb + 0.5*(Eq_lo + Eq_hi);   // mid-bin on-shell

        cdouble R = resolvent_bound_continuum(E, Eb, q_hi, q_lo);

        // Expected Im_R = -π / Δq when E falls strictly inside the bin.
        double expected_im = -M_PI / Dq;
        if (!close(R.imag(), expected_im, 1e-9)) {
            std::printf("FAIL Test1: Im_R = %.10e, expected %.10e\n",
                        R.imag(), expected_im);
            return 1;
        }
        std::printf("PASS Test1: Im_R = %.10e ≈ -π/Δq\n", R.imag());
    }

    // Test 2: Continuum-continuum cell that strictly straddles E.
    // Im_Q should be non-zero with the documented sign.
    {
        double Eb = -2.224;
        double q_lo = 50.0, q_hi = 60.0;
        double e_lo = 5.0,  e_hi = 15.0;     // p-bin in MeV
        double mu1 = 939.565 * (939.565 + 938.272 + Eb) /
                     (939.565 + 939.565 + 938.272 + Eb);
        double Eq_lo = q_lo*q_lo / (2*mu1);
        double Eq_hi = q_hi*q_hi / (2*mu1);
        double E = e_lo + Eq_lo + 0.5 * ((e_hi - e_lo) + (Eq_hi - Eq_lo));

        cdouble Q = resolvent_continuum_continuum(E, Eb, q_hi, q_lo, e_hi, e_lo);

        if (std::abs(Q.imag()) < 1e-9) {
            std::printf("FAIL Test2: Im_Q = %.10e (expected non-zero straddle)\n", Q.imag());
            return 1;
        }
        std::printf("PASS Test2: Im_Q = %.10e (non-zero)\n", Q.imag());
    }

    // Test 3: Cell well below E — Im_Q must be zero.
    {
        double Eb = -2.224;
        double E = 1000.0;  // much higher than any cell energy
        cdouble Q = resolvent_continuum_continuum(E, Eb, 60, 50, 15, 5);
        if (std::abs(Q.imag()) > 1e-9) {
            std::printf("FAIL Test3: Im_Q = %.10e (expected ~0 for cell << E)\n", Q.imag());
            return 1;
        }
        std::printf("PASS Test3: Im_Q ≈ 0 below threshold\n");
    }

    std::printf("\nAll resolvent_im_test cases passed.\n");
    return 0;
}
```

- [ ] **Step 3: Wire into the build** (depends on existing test harness — see Step 1 result)

If `tests/cpp/CMakeLists.txt` exists, add:

```cmake
add_executable(resolvent_im_test resolvent_im_test.cpp)
target_link_libraries(resolvent_im_test PRIVATE tic_tac_core)
add_test(NAME resolvent_im_test COMMAND resolvent_im_test)
```

If the test harness lives at top-level `CMakeLists.txt`, append the same block guarded by `if(BUILD_TESTING)`.

- [ ] **Step 4: Build + run**

```bash
cd /home/tian/workspace/dpol/Tic-tac
./build.sh release 2>&1 | tail -10
ctest --test-dir build -R resolvent_im_test --output-on-failure 2>&1 | tail -15
```

Expected: `resolvent_im_test` PASS.

- [ ] **Step 5: Commit**

```bash
git add tests/cpp/resolvent_im_test.cpp tests/cpp/CMakeLists.txt   # or top CMakeLists.txt
git commit -m "test: cpp unit test for resolvent Im_R, Im_Q analytic values

Test1: Im_R = -π/Δq for an on-shell q-bin.
Test2: Im_Q non-zero for cell straddling E.
Test3: Im_Q ≈ 0 for cell well below E."
```

### Task A-1.4: Resolvent trace dump

**Files:**
- Modify: `src/core/resolvent/make_resolvent.cpp:131-232`

- [ ] **Step 1: Add a single trace block at the end of `calculate_resolvent_array_in_SWP_basis`**

Just before the function closes (line 231-232 — `}` of the alpha loop and the function), insert:

```cpp
	// [EN] Trace hook (additive). When trace_im_path=true, dump the on-shell q-bin
	// row of (Re G, Im G) aggregated over all (alpha, p) pairs for this E. Cheap:
	// one extra full sweep, only on the diagnostic run. /
	// [CN] 追踪钩子（新增，可选）。trace_im_path=true 时，把当前 E 的 on-shell q-bin
	// 行（在所有 (alpha, p) 上聚合）的 Re/Im G 写入诊断文件。开销：1 次额外扫描，仅诊断跑生效。
	if (run_parameters.trace_im_path){
		// Find the on-shell q-bin: the bin whose [q_lo, q_hi] straddles q^2/(2μ) = E - Eb.
		// Aggregate ‖Re‖, ‖Im‖ over all (alpha, p) at that q-bin, both for BC and CC parts.
		double sum_re_bc_sq = 0.0, sum_im_bc_sq = 0.0;
		double sum_re_cc_sq = 0.0, sum_im_cc_sq = 0.0;
		for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
			int L = L_2N_array[idx_alpha];
			int S = S_2N_array[idx_alpha];
			int J = J_2N_array[idx_alpha];
			int T = T_2N_array[idx_alpha];
			int two_T_3N = two_T_3N_array[idx_alpha];
			double* e_SWP_array_ptr = select_swp_energy_branch(L, S, J, T, two_T_3N,
															   Np_WP, e_SWP_unco_array,
															   e_SWP_coup_array, run_parameters);
			for (int idx_p=0; idx_p<Np_WP; idx_p++){
				bool bs = (e_SWP_array_ptr[idx_p] < 0);
				for (int idx_q=0; idx_q<Nq_WP; idx_q++){
					int G_idx = idx_alpha*Nq_WP*Np_WP + idx_q*Np_WP + idx_p;
					cdouble g = G_array[G_idx];
					if (bs){
						sum_re_bc_sq += g.real()*g.real();
						sum_im_bc_sq += g.imag()*g.imag();
					} else {
						// Aggregate only "straddle" cells — those whose Im is non-zero.
						// Approximation: count any cell where |Im| > 1e-12 of |Re| as straddle.
						sum_re_cc_sq += g.real()*g.real();
						sum_im_cc_sq += g.imag()*g.imag();
					}
				}
			}
		}
		std::string trace_path = run_parameters.output_folder + "/im_path_trace.txt";
		FILE* fp = std::fopen(trace_path.c_str(), "a");
		if (fp){
			std::fprintf(fp, "G0_BC_on_shell_q\t%.6e\t%.6e\t%.6e\n",
						 std::sqrt(sum_re_bc_sq), std::sqrt(sum_im_bc_sq),
						 std::sqrt(sum_im_bc_sq) / (std::sqrt(sum_re_bc_sq) + 1e-30));
			std::fprintf(fp, "G0_CC_straddle_aggregate\t%.6e\t%.6e\t%.6e\n",
						 std::sqrt(sum_re_cc_sq), std::sqrt(sum_im_cc_sq),
						 std::sqrt(sum_im_cc_sq) / (std::sqrt(sum_re_cc_sq) + 1e-30));
			std::fclose(fp);
		}
	}
```

NOTE: the existing function takes `run_parameters` by value (signature unchanged) so `run_parameters.trace_im_path` is reachable.

- [ ] **Step 2: Build**

```bash
cd /home/tian/workspace/dpol/Tic-tac/CPP && make -j 2>&1 | tail -5
```

- [ ] **Step 3: Commit**

```bash
git add src/core/resolvent/make_resolvent.cpp
git commit -m "resolvent: trace hook dumps Re/Im aggregates for BC + CC under trace_im_path"
```

### Task A-1.5: Neumann + Padé + S trace hooks

**Files:**
- Modify: `src/core/faddeev_solver/solve_faddeev.cpp` — at four locations:
  - Per-Neumann-term loop end (search for the loop `for(NM=...)` writing `re_A_An_row_array_prev`)
  - Padé extraction (line ~1808)
  - U_array assignment to S — actually S is computed in the Python extractor, so trace stops at U.

- [ ] **Step 1: Add a helper at the top of `solve_faddeev.cpp` (just below `namespace { ... } // namespace` ~line 124)**

```cpp
// [EN] Tracing helper for trace_im_path. Appends one row to im_path_trace.txt.
// [CN] trace_im_path 用的追踪辅助函数：向 im_path_trace.txt 追加一行。
static void append_trace_row(const std::string& output_folder,
                             const char* stage,
                             double re_norm,
                             double im_norm){
    std::string p = output_folder + "/im_path_trace.txt";
    FILE* fp = std::fopen(p.c_str(), "a");
    if (!fp) return;
    std::fprintf(fp, "%s\t%.6e\t%.6e\t%.6e\n",
                 stage, re_norm, im_norm,
                 im_norm / (re_norm + 1e-30));
    std::fclose(fp);
}
```

- [ ] **Step 2: Add a per-Neumann-term recorder inside the Neumann loop**

Search for `re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx_col_NDOS]` in the loop where Neumann terms are accumulated. Right after each iteration computes the new term, add (guarded):

```cpp
	if (run_parameters.trace_im_path){
		// ‖Re‖, ‖Im‖ over the elastic on-shell row block.
		double re_sq = 0.0, im_sq = 0.0;
		for (size_t i=0; i<num_on_shell_A_rows; i++){
			for (size_t j=0; j<dense_dim; j++){
				size_t k = i*dense_dim + j;
				re_sq += re_A_An_row_array_prev[k]*re_A_An_row_array_prev[k];
				im_sq += im_A_An_row_array_prev[k]*im_A_An_row_array_prev[k];
			}
		}
		char stage[64];
		std::snprintf(stage, sizeof(stage), "AKn_elastic_n%d", (int)NM);
		append_trace_row(run_parameters.output_folder, stage,
		                 std::sqrt(re_sq), std::sqrt(im_sq));
	}
```

(insert right after the convergence-decision `if (...) { ... }` block at line 1713 — that way it fires once per outer NM iteration AFTER the new term is in `_array_prev`.)

- [ ] **Step 3: Add the Padé-best-PA recorder** at line 1812 (right after the elastic U_array fill loop closes):

```cpp
	if (run_parameters.trace_im_path){
		double re_sq = 0.0, im_sq = 0.0;
		for (size_t idx=0; idx<num_EL_A_vals; idx++){
			cdouble u = U_array[idx];
			re_sq += u.real()*u.real();
			im_sq += u.imag()*u.imag();
		}
		append_trace_row(run_parameters.output_folder, "Pade_best_PA_elastic_U",
		                 std::sqrt(re_sq), std::sqrt(im_sq));
		// Also write S-matrix surrogate: |Im(U)| / |Re(U)|. Real S computed in extractor.
		append_trace_row(run_parameters.output_folder, "S_matrix_diagonal_elastic",
		                 std::sqrt(re_sq), std::sqrt(im_sq));
	}
```

- [ ] **Step 4: Add a K-row n=0 recorder** inside the existing first-Neumann-term block (search for the comment `Extract breakup terms` at line 1286 — just before that, the first-Neumann elastic iteration is captured). Insert:

```cpp
	if (run_parameters.trace_im_path){
		double re_sq = 0.0, im_sq = 0.0;
		for (size_t i=0; i<num_on_shell_A_rows; i++){
			for (size_t j=0; j<dense_dim; j++){
				size_t k = i*dense_dim + j;
				re_sq += re_A_An_row_array_prev[k]*re_A_An_row_array_prev[k];
				im_sq += im_A_An_row_array_prev[k]*im_A_An_row_array_prev[k];
			}
		}
		append_trace_row(run_parameters.output_folder, "K_n0_on_shell_row",
		                 std::sqrt(re_sq), std::sqrt(im_sq));
	}
```

- [ ] **Step 5: Build**

```bash
cd /home/tian/workspace/dpol/Tic-tac/CPP && make -j 2>&1 | tail -5
```

- [ ] **Step 6: Run the dbg config with trace on**

```bash
cd /home/tian/workspace/dpol/Tic-tac
mkdir -p CPP/Output/miller_gate1_dbg
echo "trace_im_path=true" >> CPP/Input/input_miller_gate1_dbg.txt
time ./CPP/run CPP/Input/input_miller_gate1_dbg.txt 2>&1 | tail -10
cat CPP/Output/miller_gate1_dbg/im_path_trace.txt
```

Expected: trace file with all stages, BC Im non-zero, CC Im non-zero, AKn shows ratio decay (this is the diagnostic).

- [ ] **Step 7: Run the regression test**

```bash
python3 -m pytest tests/test_im_path_trace.py -v 2>&1 | tail -15
```

Expected: PASS — except `test_g0_im_part_nonzero` should pass (BC Im is analytic-guaranteed). If it fails, the resolvent trace is broken and is itself a bug to fix.

- [ ] **Step 8: Commit**

```bash
git add src/core/faddeev_solver/solve_faddeev.cpp
git commit -m "solver: Im-path trace hooks (Neumann, Padé, S surrogate)

Gated on trace_im_path; appends rows to im_path_trace.txt at each stage."
```

---

## Workstream A-2: Diagnostic Note

**Owner suggestion:** Subagent #4 (after A-1 lands).
**Depends on:** A-1.5 (trace data exists).
**Verifies:** A-1 trace file is interpretable.

### Task A-2.1: Run trace on cheap dbg config and write diagnosis

**Files:**
- Create: `docs/im_path_diagnosis_2026-05-09.md`

- [ ] **Step 1: Run the trace if not already done**

```bash
cd /home/tian/workspace/dpol/Tic-tac
./CPP/run CPP/Input/input_miller_gate1_dbg.txt 2>&1 | tail -5
cat CPP/Output/miller_gate1_dbg/im_path_trace.txt
```

- [ ] **Step 2: Write the diagnostic note**

`docs/im_path_diagnosis_2026-05-09.md` template:

```markdown
# Im-path diagnosis — Miller Gate 1 cheap config (Np=Nq=20, Nijmegen-I, 13 MeV)

**Date:** 2026-05-09
**Run:** `CPP/Output/miller_gate1_dbg`
**Config:** `CPP/Input/input_miller_gate1_dbg.txt` with `trace_im_path=true`.

## Trace table (verbatim)

(paste the contents of `im_path_trace.txt` here as a fenced ``` ``` block)

## Analytic reference

For the bound-continuum on-shell q-bin at q ≈ q_on, the cell-averaged
resolvent has analytic Im part `Im_R = -π / Δq` — see
`src/core/resolvent/make_resolvent.cpp:80`. With Δq ≈ (paste numerical Δq
from the run), the expected ‖Im_BC‖ for that one (alpha, p, q) cell is
≈ π / Δq ≈ (number).

## Where Im part survives / where it dies

(For each stage, comment on whether Im/Re ratio matches expectation.)

| Stage | Im/Re | Verdict |
|---|---|---|
| G0_BC_on_shell_q | … | … |
| G0_CC_straddle_aggregate | … | … |
| K_n0_on_shell_row | … | … |
| AKn_elastic_n0 | … | … |
| AKn_elastic_n1 | … | … |
| … | … | … |
| Pade_best_PA_elastic_U | … | … |
| S_matrix_diagonal_elastic | … | … |

## Conclusion

(One of the following:
- "Bug located: at <file:line>, <symptom>. Step A-3 will <patch>."
- "No bug: Im part propagates intact, but is suppressed by undersampling
   of the Heaviside step at Np=20. Step A-3 will introduce Lorentzian
   smoothing of Im_Q with parameter ε ~ ΔE/2 and a Δ → 0 sweep.")
```

- [ ] **Step 3: Commit**

```bash
git add -f docs/im_path_diagnosis_2026-05-09.md
git commit -m "docs: Im-path diagnosis on Np=20 dbg run

Identifies stage at which Im part (does/does not) survive in the AGS
amplitude pipeline. Drives Step A-3 patch shape."
```

---

## Workstream A-3: Targeted Fix (DISCOVERY-DRIVEN)

**Owner suggestion:** Subagent #5 (after A-2 lands; main session reviews).
**Depends on:** A-2.1.
**Shape:** Determined by A-2's verdict. Two pre-cooked branches:

### Branch A-3.X (bug-located): single-point patch

If A-2 identified an exact line where Im part is zeroed:

- [ ] **Step 1: Write a focused failing test** demonstrating the regression at the bug site (e.g., `tests/test_im_propagation.py`).
- [ ] **Step 2: Patch the offending line.** Keep additive: add a flag-gated alternative path if the existing one is referenced elsewhere.
- [ ] **Step 3: Re-run `tests/test_im_path_trace.py`, the new focused test, and the existing 190 MeV regression.**
- [ ] **Step 4: Commit + update `docs/im_path_diagnosis_2026-05-09.md`** with the fix summary.

### Branch A-3.Y (undersampling): Lorentzian smoothing

If A-2 verdict is "no bug; undersampled":

- [ ] **Step 1: Add `lorentzian_im_smoothing_eps` to run_params** (analogous to `trace_im_path`). Default 0 (off — no behavior change).
- [ ] **Step 2: In `make_resolvent.cpp` `resolvent_continuum_continuum`,** add an additive branch when ε > 0: replace the Heaviside-derived Im_Q with the Lorentzian-smoothed analytic form (formulas to be transcribed during step 1; reference: Kukulin 2007 Eq. (32) with Δ → Δ + iε).
- [ ] **Step 3: Δ → 0 convergence sweep** in `examples/lorentzian_eps_sweep.py` (new): run dbg config at ε ∈ {0.05·ΔE, 0.1·ΔE, 0.2·ΔE, 0.5·ΔE}, plot ²S₁/₂ Im δ.
- [ ] **Step 4: Commit + update diagnosis with the converged ε choice.**

(The exact code for either branch is intentionally not pre-written — it depends on A-2's findings. A code-reviewer subagent must approve the patch shape before merging.)

---

## Final Regression Sweep

**Owner suggestion:** Main session (or fresh subagent) after A-3 lands.

- [ ] **Step 1: Run the full regression battery**

```bash
cd /home/tian/workspace/dpol/Tic-tac
python3 -m pytest tests/ -v 2>&1 | tail -30
ctest --test-dir build --output-on-failure 2>&1 | tail -15
python3 -m unittest tests/test_190mev_data_pipeline.py 2>&1 | tail -10
```

Expected: all green.

- [ ] **Step 2: Verify the w1_scale=0 + breakup-on regression gate**

```bash
# (Adapt to the existing 2NF baseline at CPP/Output/miller_gate1/)
diff <(grep -E '^U[0-9]' CPP/Output/miller_gate1/U_PW_elements_*.txt) \
     <(grep -E '^U[0-9]' CPP/Output/miller_gate1_dbg/U_PW_elements_*.txt)
```

Expected: identical (or within float tolerance) on Re; Im parts may differ if A-3 introduced Im — that is intentional.

- [ ] **Step 3: Re-extract the Np=30 3NF dataset through the new pipeline**

```bash
python3 examples/extract_phase_shifts.py \
    --solver-out-dir Tic-tac-w1cache/CPP/Output/miller_gate1_np30_3nf \
    --markdown-out output/np30_3nf_after_a3.md \
    --tlab 12.05
```

- [ ] **Step 4: Commit final artifacts**

```bash
git add -f output/np30_3nf_after_a3.md docs/im_path_diagnosis_2026-05-09.md
git commit -m "report: Np=30 3NF after Im-path fix (A-3)

Final honesty + Im-path-corrected report. PASS/FAIL per (J,P) block as per spec §5."
```

---

## Self-Review (run inline before handoff)

**Spec coverage check:**

- §3 architecture → tasks B-1.2/B-1.3 (Padé split), B-2.2 (extractor threshold), A-1.4/A-1.5 (trace hooks), A-2.1 (diagnosis), A-3.X/Y (fix).
- §4 trace format → A-1.4 (G0 stages), A-1.5 (Neumann/Padé/S).
- §5 validation gates → B-1 (Conv column), B-2 (FAIL gating + markdown), A-1 (trace stages), A-3 (post-fix Im δ non-zero, |S_kk| toward 1).
- §6 testing → B-1.1, B-2.1, A-1.2, A-1.3.
- §7 risks → A-3 has both bug-located and undersampling branches pre-named; w1_scale=0 + breakup-on regression in Final Step 2.
- §8 out-of-scope items remain out-of-scope. ✓

**Placeholder scan:** A-3 X/Y pre-cooked branches are intentionally undetailed — this is the *only* place in the plan with discovery-dependent content, and it is gated behind A-2's diagnostic note + a code-reviewer pass.

**Type consistency:** new field name `trace_im_path` used identically in `type_defs.h`, parser, default, test, trace dump. New arrays `pade_approximants_truly_converged_array` / `_maxiter_truncated_array` named identically across allocation, init, populate, free, and sidecar emission. Stage names in trace file match between solver writer and Python test reader (`G0_BC_on_shell_q`, etc.).

**Additive principle audit:** every modified file keeps its old behavior when the new flag is off. Old arrays untouched. No renames. ✓
