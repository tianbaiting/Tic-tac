# Tic-tac Breakup-Channel & 3NF Validation Honesty — Design

**Date**: 2026-05-09
**Status**: Approved (brainstorming gate)
**Scope**: Combined project — diagnostic honesty pass (B-min) + Im-path diagnosis & targeted fix (A) for the WPCD Faddeev solver.

## 1. Project Definition

Make Tic-tac able to **honestly** report 3NF + breakup-channel results. Two coupled deliverables, one spec, one plan, dependency-ordered execution:

- **B-min** — remove the two most blocking pieces of fake-validation (Padé max-iter silently labelled "converged"; unitarity defect not gated).
- **A** — diagnose why the resolvent's analytic Im parts (which DO exist in `make_resolvent.cpp`) fail to surface as a non-zero elastic Im δ, then apply a targeted fix.

Out of scope (deferred to a follow-up spec): silent LEC injection logging, float32 cache regression test, deprecation of the `check_3nf_normalization` tool, 190 MeV/u full observable sweep, three-body bound-state driver.

## 2. Background and Motivation

A critical audit of the post-3NF state surfaced seven concerns. The two most damaging:

| # | Concern | Citation |
|---|---|---|
| 6 | `solve_faddeev.cpp:1701` sets `pade_approximants_conv_array[idx] = true` whenever `NM == NM_max`, with no distinction between true convergence and max-iteration timeout. The Np=30 3NF JP=1+ run reported in memory used a `P[14,14]` *forced* extraction. | `src/core/faddeev_solver/solve_faddeev.cpp:1701-1710` |
| 7 | `extract_phase_shifts.py:287` prints `||SS†−1||` as a diagnostic but applies no threshold. The Np=30 JP=1- block (where the Ay puzzle physics lives) has `||SS†−1|| = 1.18` (100%+ unitarity violation) yet contributes equally to "Gate 1 PASS". | `examples/extract_phase_shifts.py:283-287` |

The breakup-channel investigation revealed an unexpected fact: the resolvent already contains the analytic +iε prescription as a cell-averaged Heaviside step:

- `resolvent_bound_continuum` (lines 80-81): `Im_R = ±π/Δq` for the on-shell q-bin.
- `resolvent_continuum_continuum` (lines 123-127): `Im_Q ∝ π/(Δp·Δq)` for cells that straddle the threshold.

Yet `extract_phase_shifts.py` consistently reports `Im δ = 0` at Tlab = 13 MeV (well above the np breakup threshold). This is a **pipeline pathology**, not a missing-physics gap: somewhere between G₀ and the final phase-shift extraction the Im part is being silenced or undersampled. Diagnosing where is the core technical work of A.

## 3. Architecture and Components

```
B-min (diagnostic honesty layer)
├─ src/core/faddeev_solver/solve_faddeev.cpp
│   └─ Split pade_approximants_conv_array into:
│        pade_approximants_truly_converged_array
│        pade_approximants_maxiter_truncated_array
│      Same for the BU twin arrays.
└─ examples/extract_phase_shifts.py
    └─ Read the truncated flag from solver output; reject PASS when set.
       Apply ‖SS†−1‖ > 0.2 hard threshold per (J,P) block.

A-step1 (Im-path trace)
├─ src/core/resolvent/make_resolvent.cpp
│   └─ Optional dump of Re/Im G₀ on the on-shell q-bin row when run_parameters.trace_im_path == true
├─ src/core/faddeev_solver/solve_faddeev.cpp
│   └─ Per-Neumann-term ‖Re(A·Kⁿ)‖ / ‖Im(A·Kⁿ)‖ recording on the elastic on-shell row
├─ src/core/faddeev_solver/solve_faddeev.cpp (Padé extraction)
│   └─ Verify the chosen best PA preserves Im part (compare to truncated double of the same)
└─ examples/extract_phase_shifts.py
    └─ Verify Im δ = ½ arg(S) reads the Im part of the on-shell U correctly

A-step2 (targeted fix; shape determined by step 1)
└─ Single-point patch (5–50 lines).
   Most likely candidates, ordered by current suspicion:
     (a) Re-only projection somewhere in the elastic extraction path
     (b) Im part dropped during a symmetrization or transpose step
     (c) Heaviside step undersamples at Np=20-30; remediation = Lorentzian smoothing of the
         cell-averaged Im_Q with parameter ε ~ ΔE/2
   Step 2's exact shape is intentionally NOT pre-specified — it depends on what step 1 finds.
```

## 4. Data Flow — A-step1 Trace Format

A new boolean input parameter `trace_im_path` (default `false`) gates a single-shot diagnostic file `im_path_trace.txt` written to the run output folder. The file is a flat tab-separated table with one row per stage:

```
stage                                            ‖Re‖           ‖Im‖           Im/Re
G0_BC_on_shell_q                                 …              …              …
G0_CC_straddle_aggregate                         …              …              …
K_n=0_on_shell_row                               …              …              …
A·K^0_elastic                                    …              …              …
A·K^1_elastic                                    …              …              …
…
A·K^N_elastic                                    …              …              …
Pade_best_PA_elastic_U                           …              …              …
S_matrix_diagonal_elastic                        …              …              …
```

The location where `Im/Re` collapses from ~1 to ~0 is the bug site.

Reference run for trace: Np=Nq=20 Nijmegen-I, Tlab=13 MeV, JP=1+, no 3NF, breakup-flag-on. ~3 minute wall time on the cheap config (`CPP/Input/input_miller_gate1.txt`).

## 5. Validation Gates (Definition of Done)

### B-min

1. `truly_converged` and `maxiter_truncated` flags are emitted in solver output and consumed by `extract_phase_shifts.py`.
2. The phase-shift markdown report explicitly marks max-iter rows ⚠ and rejects them from PASS.
3. `||SS†−1|| > 0.2` per (J,P) block triggers FAIL in the markdown report (numerical value still recorded).
4. Re-running the existing Np=30 3NF dataset through the new extractor produces a markdown report that explicitly says: "JP=1+ Padé maxiter-truncated; JP=1- ²P block ‖SS†−1‖=1.18 FAIL".

### A-step1

1. `im_path_trace.txt` exists for the reference run and contains every listed stage.
2. A diagnostic note is written to `docs/im_path_diagnosis_2026-05-09.md` identifying either:
   - the precise line where Im part is dropped (with file:line citation), OR
   - confirmation that the Im part propagates intact and the physical Im δ is genuinely zero at Np=20-30 due to undersampling.

### A-step2 (conditional on step 1's verdict)

1. After patching, the reference run produces non-zero ²S₁/₂ Im δ at Tlab=13 MeV.
2. `‖SS†−1‖` for JP=1+ at Np=30 3NF drops below 0.124 (the current value) — i.e. the patch reduces the unitarity defect in the physically-correct direction (sub-unitary, indicating BU absorption).
3. Regression: w1_scale=0 + breakup-on remains bit-for-bit identical (within float tolerance) to the 2NF-only baseline on elastic Re δ.

## 6. Testing Strategy

- **Existing regression**: `python3 -m unittest tests/test_190mev_data_pipeline.py` continues to pass.
- **New regression test**: `tests/test_im_path_trace.py` runs `CPP/run` on the cheap 13 MeV input with `trace_im_path=true`, asserts that:
  - the trace file is created and non-empty;
  - every named stage row exists;
  - `Im/Re` ratios at the resolvent stages match the analytic Heaviside expectation within 5%.
- **New C++ unit test**: `tests/cpp/resolvent_im_test.cpp` constructs a known straddle configuration and asserts the sign and magnitude of `Im_R` and `Im_Q` against hand-derived values.
- **Hand-calculation oracle**: a ~half-page derivation in `docs/im_path_diagnosis_2026-05-09.md` showing the analytic Im_R = −π/Δq for a 13 MeV q_on bin, used as the trace's first-row reference.

## 7. Risk and Error Handling

| Risk | Mitigation |
|---|---|
| Step 1 finds Im part propagates intact and Im δ really is zero at Np=20-30 | Step 2 escalates to Lorentzian smoothing of `Im_Q` (parameter ε ~ ΔE/2 with a Δ → 0 convergence sweep). |
| Step 2 patch shifts Re δ enough to break the Np=30 3NF Gate 1 +104° claim | This was already a forced-PA result; B-min flags it as ⚠ regardless. Use w1_scale=0 + breakup-on regression to certify the 2NF baseline is unchanged. Any genuine 3NF-elastic shift is a *physical improvement*, not a regression. |
| Existing BU U-matrix extraction assumes Im=0 in its normalization or naming | Step 1 trace will surface the assumption. Step 2 scope expands to include the BU extraction path. |
| Production binary (`CPP/run`) and CMake binary (`bin/tic-tac`) drift after the patch | Both build paths must be exercised by the regression tests; the existing two-build-system invariant is documented in `Tic-tac/CLAUDE.md` and respected. |

## 8. Explicitly Out of Scope

- Logging silent LEC injection in `three_nucleon_force_model.cpp:26-29` (item 1 of the audit).
- Renaming the `include_breakup_channels` input flag (item 2).
- Float32 vs double regression for the W1 cache (item 4).
- Deprecating the `check_3nf_normalization` reference table (item 5).
- 190 MeV/u dataset re-run with iT11/T20/T22 sweep.
- Triton bound-state driver (Gate 3, BLOCKED).

These will be folded into a follow-up spec once the present project's diagnostic infrastructure is in place.

## 9. Dependency Graph for Implementation

```
B-min (Padé honesty)        ──┐
B-min (unitarity threshold) ──┤── independent of each other; can run in parallel
                              │
A-step1 (trace + diagnosis)   ── independent of B-min; can run in parallel
                              │
                              ▼
A-step2 (targeted fix)        ── blocked by A-step1; shape determined by step1's verdict
                              │
                              ▼
Final regression sweep        ── blocked by all of the above
```

Three of the four work-streams are independent and parallelizable. Subagent dispatch is appropriate after the writing-plans skill produces a per-task breakdown.
