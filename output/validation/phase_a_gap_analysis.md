# Phase A — internal gap analysis (realistic-grid campaign)

**Audit date:** 2026-08-15
**Latest commit:** `093c2fb` (block-level distributed + resumable exact W1 construction)
**Host:** 96-core Xeon, 125 GB RAM, 2.9 TB free disk
**Build:** current (`cmake --build build -j` clean); `ctest` 17/17 pass.

## What already exists (do NOT duplicate)

- **Exact factorized N2LO 3NF** (`chiral_N2LO_full_factorized`): c1/c3/c4/cD/cE, validated
  against independent Python + five-angle reference + Golak Table 2.
- **Exact W1 Hermitian reuse**: 56 evaluate + 45 transpose-fill per 101-block J=1/2 sector.
- **Block-level distributed/resumable W1**: `w1_worker {plan|build|status|assemble|verify}`
  with `--shard K/N` (sector-level) and `--worker-index I --worker-count N` (block-level
  modulo).  Atomic shard publication, schema-v9 full provenance, SHA-256 fingerprint,
  verified bitwise-identical to monolithic `W1_PW_cache::build`.
- **`two_J_3NF_force_max`** cutoff: high-J sectors run pure-2NF, W1 built only for J<=Jc.
- **Realistic candidate grid**: `Np=82, Nq=12, p=(0.88,325), q=(1,90)`,
  `Tlab_midpoint=10.045 MeV`, deuteron binding within 20-keV gate.
- **Realistic 2NF J=1/2 both-parity solve** (feasibility JSON): P123 ~64 s, solve ~152 s,
  16/16 amplitudes converged at [24/24] and [32/32], max |ΔU| = 1.28e-9 MeV.
- **J3NF convergence driver** `examples/run_j3nf_convergence.py` (thin orchestrator).
- **Observable reconstruction** `python/tictac` (amplitudes, observables, io).
- **Deuteron binding p-grid ladder** `output/validation/idaho_n3lo_deuteron_binding_ladder.json`.
- **Stable Padé** resummation (`5e003d6`): denominator-coefficient solve, no tolerance relaxed.

## Gaps this campaign must close

| Gap | Phase | Status before this work |
|---|---|---|
| No way to build a SELECTED SUBSET of evaluate blocks (only modulo partition) | B | **missing** — must extend `w1_worker` |
| Realistic-grid W1 block timing NEVER measured (19.4 h is a linear extrapolation from Np=4,Nq=3) | B | **missing** — the headline new science |
| Load imbalance of static modulo partition unknown at realistic grid | C | **missing** |
| No Nq / J_2N_max / two_J_3N_max convergence ladder at realistic grid (only deuteron binding + single 2NF solve) | D | **partial** |
| No observable-level (dsigma, Ay, T20, T21, T22) convergence summary | E | **missing** |
| W1 quadrature order not tested at realistic grid (only cell-level) | F | **missing** |
| Paired 2NF+3NF realistic solve never run | G | **missing** |
| J3NF ladder driver not exercised at realistic grid | H | **driver exists, not run** |
| Zero-3NF limit and component decomposition not automated | I | **missing** |
| No CI (GitHub Actions) | J | **missing** |

## Realistic W1 cost basis (current, all EXTRAPOLATION)

All current "19.4 h / sector" numbers are linear extrapolations from the measured
(Np=4, Nq=3) 101-block sector × 6724 cell-ratio × Hermitian-reuse × 96-thread scaling.
**No realistic-grid W1 block has ever been timed.**  This is the single most important
gap.
