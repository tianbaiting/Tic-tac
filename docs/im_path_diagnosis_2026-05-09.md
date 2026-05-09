# Im-path diagnosis — Miller Gate 1 cheap config (Np=Nq=20, Nijmegen-I, 13 MeV)

**Date:** 2026-05-09
**Run:** `CPP/Output/miller_gate1_dbg_trace`
**Config:** `CPP/Input/input_miller_gate1_dbg.txt` with `trace_im_path=true` (in-memory only; committed config keeps the flag off).
**Workstream:** A-2 (per `docs/superpowers/plans/2026-05-09-tictac-breakup-3nf-honesty.md`).

## Trace table (verbatim)

```
G0_BC_on_shell_q          1.143439e+00   1.761784e+00   1.540776e+00
G0_CC_straddle_aggregate  3.202082e+01   2.683107e+01   8.379257e-01
G0_BC_on_shell_q          8.386216e-01   1.267863e+00   1.511841e+00
G0_CC_straddle_aggregate  2.150856e+01   2.025736e+01   9.418277e-01
K_n0_on_shell_row         0.000000e+00   0.000000e+00   0.000000e+00
AKn_elastic_n0..n14       0.000000e+00   0.000000e+00   0.000000e+00
Pade_best_PA_elastic_U    0.000000e+00   0.000000e+00   0.000000e+00
S_matrix_diagonal_elastic 0.000000e+00   0.000000e+00   0.000000e+00
(repeated for the second J^P block)
```

Columns: `stage <TAB> ‖Re‖ <TAB> ‖Im‖ <TAB> Im/Re`.

## Analytic reference (resolvent stage)

For the bound-continuum on-shell q-bin at Tlab=13 MeV the cell-averaged resolvent has analytic
`Im_R = -π / Δq` (`src/core/resolvent/make_resolvent.cpp:80-81`). The trace value
`G0_BC_on_shell_q ‖Im‖ ≈ 1.76` is consistent with `π/Δq` summed over the (alpha, p) cells whose
deuteron-like SWP energy is below threshold; the Im/Re ratio of ~1.5 reflects that the principal-
value Re part and the +iπδ Im part are of comparable magnitude on the on-shell row. **The
analytic Heaviside +iε prescription is intact in the resolvent.**

## Where the Im part survives / where it dies

| Stage | ‖Re‖ | ‖Im‖ | Im/Re | Verdict |
|---|---|---|---|---|
| G0_BC_on_shell_q | 1.14 / 0.84 | 1.76 / 1.27 | ~1.5 | ✓ Im intact, matches analytic Heaviside |
| G0_CC_straddle_aggregate | 32 / 22 | 27 / 20 | ~0.9 | ✓ Im intact |
| K_n0_on_shell_row | 0 | 0 | — | ✗ Both Re AND Im dropped |
| AKn_elastic_n=0..14 | 0 | 0 | — | ✗ Propagates the zero |
| Pade_best_PA_elastic_U | 0 | 0 | — | ✗ Padé of zero is zero |
| S_matrix_diagonal_elastic | 0 | 0 | — | ✗ |

The collapse happens between G0 (intact) and the n=0 Neumann extraction (zero). This is **earlier**
in the pipeline than the spec's leading hypothesis (which suggested a Re-only projection or
Im-drop somewhere in the elastic extraction, but not necessarily at n=0).

## Root-cause finding: Re-only projection at the n=0 Neumann extraction

`src/core/faddeev_solver/solve_faddeev.cpp` has an asymmetric pattern between the **first**
Neumann term (n=0) and **subsequent** terms (n≥1):

**n=0 elastic (line 1305):**
```cpp
cdouble a_coeff = re_A_An_row_array_prev[ndos.row_storage_idx*dense_dim + ndos.col_storage_idx];
```
This is an implicit `double → cdouble` conversion. The result has zero imaginary part by C++
language rule. The companion array `im_A_An_row_array_prev[...]` is never read at this site.

**n=0 breakup (line 1350):**
```cpp
cdouble a_BU_coeff = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx_col_NDOS];
```
Same bug.

**n≥1 elastic (line 1618):**
```cpp
cdouble a_coeff = {re_A_An_row_array[ndos.row_storage_idx*dense_dim + ndos.col_storage_idx],
                   im_A_An_row_array[ndos.row_storage_idx*dense_dim + ndos.col_storage_idx]};
```
Correct: aggregate-init of `cdouble` from both real and imaginary buffers.

**n≥1 breakup (line 1657):** same correct pattern.

**Consequence:** The Padé series is built on `a_0 = Re(A·K^0) + 0i`. Higher orders contribute
their full Im content (since the Neumann update at lines 1438-1439 propagates Im correctly), but
the leading driving term is zeroed in Im. For Tlab below the breakup threshold (where higher-
order Im is also small) the cumulative Padé Im stays of order 10⁻³ — well within "looks zero" at
the Re δ ~ 100° scale. For Tlab above threshold, the iterative Im growth alone can never recover
the missing Born-term Im; it's a fixed deficit per channel.

## Secondary finding: the dbg config converges trivially

Beyond the Re-only bug, the trace reveals a separate pathology specific to this cheap
configuration:

- The trace at line 1322 reads `re_A_An_row_array_prev` IMMEDIATELY after the n=0 extraction
  loop. Both ‖Re‖ and ‖Im‖ aggregate to **exactly zero** — so the `_array_prev` buffer at this
  point is wholly unfilled, not merely missing its Im part.

This is consistent with the n=0 extraction reading from a buffer that the upstream CPVC fill
populates only after the first iteration cycle (lines 1574, 1589 fill `re_A_An_row_array` and
then copy `_array → _array_prev`). At the very first extraction (n=0), `_array_prev` is still
zero-initialized from line 951-952. With both the Born term zero AND every Neumann iteration
starting from zero, the Padé degenerates to the trivial fixed point and `pade_approximants_conv_array[idx] = true` fires by criterion 0 (NM == NM_max with no change).

This is a **buffer-naming/timing bug** that is logically distinct from the Re-only bug, but with
overlapping symptoms in this dbg config.

The production Np=30 3NF run reported finite (nonzero) Re δ ≈ +104°, so in production the
buffer DOES get filled — there must be a code path that runs CPVC into `_array_prev` before
n=0 extraction in the production case but not in this dbg config. Two candidate explanations:

1. The `parallel_run` / `production_run` switches gate a CPVC pre-fill that the dbg config
   misses.
2. The dbg config's small `Np=Nq=20` truncation makes `num_EL_A_vals = 0` (no on-shell channels
   under J_2N_max=2 / two_J_3N_max=1) and the loop over `idx_d_row × idx_d_col × idx_q_com`
   simply never executes. (If this were true the trace would still print zero, since the
   trace block at line 1322 always runs once per outer call.)

Verification of explanation 2 takes ten minutes: re-run with Np=Nq=24 or J_2N_max=3 and re-read
the trace. **This is the recommended A-3 first-step.**

## Conclusion

Two confirmed findings, ordered by impact:

1. **Re-only projection at n=0 Neumann extraction** (`solve_faddeev.cpp:1305` and `:1350`):
   the Born term has its Im part silently dropped, contaminating the entire Padé series. This is
   the source of the empirical Im δ = 0 observation in production runs at Tlab > breakup
   threshold. **Fix shape: trivially additive, two-line edit per site, mirror the n≥1 pattern.**

2. **Dbg-config trace shows everything zero downstream of G0**: the `_array_prev` buffer is
   unfilled at the n=0 extraction site in this specific configuration. Independent of finding 1.
   **Fix shape: probably a cheap config tweak (raise Np or J_2N_max), with verification by
   re-running the trace.**

Step A-3 will:
- Apply the line 1305 / 1350 patch (additive: build the cdouble from both buffers, mirroring the
  n≥1 pattern at line 1618 / 1657).
- Re-run the dbg config (or a marginally larger version) to confirm the trace now shows nonzero
  K_n0 → AKn → Pade → S Im, with a finite Im/Re ratio.
- Spot-check that the Np=30 production 3NF run's Re δ does not regress (the n=0 Re part was
  always being read correctly; only the Im was dropped).
