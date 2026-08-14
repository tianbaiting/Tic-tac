# Distributed, resumable, exact W^(1) construction

This document describes the workflow-level exact parallelism and checkpointability
layer for the complete factorized N2LO 3NF `W^(1)` construction. It introduces
**no physics approximation**: no SVD, no low-rank compression, no sparsification,
no interpolation, no surrogate. The final assembled `W^(1)` is numerically
identical to the existing exact factorized implementation — bitwise identical
whenever execution order permits, proven by a cryptographic SHA-256 fingerprint.

The status of the realistic Ay campaign is unchanged: it is **not** complete.
This layer removes the *operational* blocker (a single 19.4 h, non-resumable,
single-sector build) so that the realistic-grid convergence ladder can be
attempted safely across many processes/machines.

---

## 1. Mathematical decomposition

`W^(1)` is energy-independent at the operator/grid level and block-addressable:
for a fixed basis, grid, regulator, LECs, and quadrature, the dense W^(1) of one
conserved `J^pi` sector factors into independent dense blocks, one per allowed
`(alpha_r, alpha_c)` channel pair:

```
global W1  =  deterministic assembly of independently-computed exact blocks
```

Each block is the four-dimensional radial-bin cell average defined in
`W1_PW_cache::build` (Gauss–Legendre per bin, `1/hbarc^5` conversion,
`1/sqrt(Δ)` WP normalization — see `docs/three_nf_equation_contract.md` §7).
No block depends on another block's *value*; the only cross-block relation is
the **exact Hermiticity** `W(ar,ac; pr,qr,pc,qc) = W(ac,ar; pc,qc,pr,qr)` of the
factorized model, which lets one orientation be the transpose of the other.

---

## 2. Work unit: one `(alpha_r, alpha_c)` channel block

The independent work unit is **one `(alpha_r, alpha_c)` channel block**, with a
role of `evaluate` or `transpose_fill`:

| Role | Action | Who computes it |
|------|--------|------------------|
| `evaluate`        | direct Gauss–Legendre integration via the shared `integrate_w1_channel_blocks` | a worker |
| `transpose_fill`  | exact Hermitian transpose of its evaluated conjugate | the solver, on demand |

**Why this granularity** (not smaller, not larger):

- **Matches existing provenance.** The hash-keyed HDF5 cache (`tictac::cache::*`)
  is already keyed per-block (`W1Key`, schema v9, with all LECs, grid hashes,
  quadrature orders, angular order). Reusing it means the distributed path and
  the in-process path share one provenance pipeline.
- **Matches the Hermitian contract.** Transpose reuse is per-block, so the
  planner can mark one orientation `evaluate` and the other `transpose_fill`
  without ambiguity.
- **Restart granularity.** At the realistic `Np=82, Nq=12` grid one block is
  `82·82·12·12·8 B ≈ 944 KiB` and (extrapolated from the measured 101-block
  sector) ~11 min at 96 threads / 101 blocks. Losing one block costs ~11 min of
  a ~19.4 h sector — excellent granularity, trivial filesystem overhead
  (~101 shards × ~1 MiB per sector).
- **No global lock during integration.** Each worker writes its *own* shard
  file via atomic temp+rename; there is no multi-writer HDF5 contention and no
  central scheduler.
- **MPI/process independence.** A static partition `global_eval_index % worker_count
  == worker_index` assigns units with zero cross-process communication.

**Canonical-orientation rule** (deterministic, cache-state-independent): for
each Hermitian pair `{(a_r,a_c), (a_c,a_r)}` the member with the smaller
sector-local block id is `evaluate`; the other is `transpose_fill`. Diagonals
are `evaluate`. This is the absolute form of the rule already used inside
`W1_PW_cache::build` (`blk < reverse`), so pre-building every `evaluate` block
makes the solver a pure cache-hit + cheap transpose-fill pass with **zero
expensive W^(1) evaluation**. Workers therefore never duplicate work and never
evaluate both orientations of a pair.

---

## 3. Architecture

```
W1WorkPlan      determines the exact work units (pure: pw + hermitian flag)
W1WorkUnit      immutable description of one piece (id, a_r/a_c, role, conjugate)
W1BlockExecutor computes one unit (shared integrate_w1_channel_blocks)
W1Store         atomically persists completed units (tictac::cache::store_w1)
W1Assembler     reconstructs/verifies the full W1 (SHA-256 fingerprint)
W1Manifest      records physics/numerics/schema provenance + completion state
```

Source files (all under `src/interactions/` and `tools/w1_worker/`):

| File | Role |
|------|------|
| `w1_integrate.{h,cpp}` | `build_per_bin_quadrature`, `integrate_w1_channel_blocks`, `W1BlockExecutor` — the single source of truth for the per-cell Gauss–Legendre integration. |
| `w1_work_plan.{h,cpp}` | `W1WorkPlan` (enumerate + classify), `make_channel_view`. |
| `w1_manifest.{h,cpp}` | `W1Signature`, `make_w1_signature`, `signature_hash`, `write_manifest_json`, `W1Assembler::fingerprint`. |
| `w1_pw_cache.{h,cpp}` | The monolithic dense builder; now calls `integrate_w1_channel_blocks` (refactor — bitwise-identical). |
| `tools/w1_worker/w1_worker.cpp` | CLI: `plan | build | status | assemble | verify`. |

**Key design invariant:** both the monolithic `W1_PW_cache::build` and the
distributed `W1BlockExecutor` call the *same* `integrate_w1_channel_blocks`, so
the two paths produce bitwise-identical blocks **by construction** — not by
after-the-fact testing. The exactness test (`tests/cpp/test_w1_distributed.cpp`)
is belt-and-suspenders.

The factorized matrix-element implementation
(`chiral_N2LO_3NF_factorized::W1_element` / `evaluate_factorized_element`) is
**unchanged** and remains usable independently of distributed execution.

---

## 4. Provenance / cache contract

Every persisted block carries the full `W1Key` (schema v9):

```
schema_version, potential_model, tnf_model,
Np_WP, Nq_WP, J_2N_max, two_J_3N_max, two_J_3N, P_3N, a_r, a_c,
c_D, c_E, Lambda_3NF, c_1, c_3, c_4,
g_A, f_pi_MeV, m_pi_MeV, Lambda_chi_MeV, hbarc_MeV_fm, regulator_kind,
chebyshev_s, chebyshev_t, tensor_force, isospin_breaking_1S0,
p_grid_hash, q_grid_hash,            # SHA-256 of the exact boundary arrays
Np_per_WP_W1, Nq_per_WP_W1, Nangle_3NF
```

The manifest adds `two_J_3NF_force_max` (the active-sector selector, not in
`W1Key`) to the signature hash. Two campaigns are interchangeable **only** if
their signature hashes agree. A stale/incompatible shard is rejected: `read_w1_h5`
re-validates `schema_version`, `key_hash_full`, kind, and shape on every read;
a mismatch is a clean miss, never a silent reinterpretation.

---

## 5. Crash safety

Block publication is atomic (Phase E, already in `cache_io_w1.cpp`):

```
write unique temp file in the same dir (same filesystem)
  -> flush
  -> validate metadata + payload
  -> atomic rename() onto the final path
```

A killed process never leaves a partial block at its final path. The worst case
is an orphaned `*.tmp.<pid>.<seq>` file, which no reader consults. Two workers
building the same block each get a distinct temp name; the second rename simply
overwrites the first with another valid file. `read_w1_h5` re-checks the key
hash, so a corrupt file can never be silently accepted (verified by
`test_w1_distributed`'s corrupt-shard test).

**One immutable shard per work unit** (not concurrent writes into one HDF5 file).
Correctness and recoverability are prioritized over having one file during
construction.

---

## 6. Hermitian reuse

Preserved exactly. Workers compute **only** `evaluate` units (the canonical
triangle). `transpose_fill` units are produced at solve time by the existing
`W1_PW_cache::build` Hermitian logic (load evaluate hit → transpose-fill reverse
→ store). Workers never evaluate both orientations. `assemble` verifies the
contract: every reverse block present in the cache equals the exact transpose
of its conjugate (measured: 90/90 pairs, `bad=0`).

For non-Hermitian models (`W1_is_exactly_hermitian() == false`), every block is
`evaluate` (no transpose reuse), matching the `W1_PW_cache::build` fallback.

---

## 7. Resumption semantics

A worker's per-unit loop is:

```
for each evaluate-unit in my partition:
    if lookup_w1(key).hit:  skip          # already complete
    else:                   integrate + store_w1 (atomic)
```

So:

- **already-complete units** → loaded/skipped (cache hit);
- **missing units** → computed;
- **half-written units** → rejected as misses (§5) and recomputed.

No completed expensive integral is recomputed. Demonstrated (§10): Run A builds
7 blocks then stops; Run B resumes with 16 workers and worker-0 reports
`built=0 cache_hits=7`.

---

## 8. Multi-worker usage

Two complementary distribution modes:

### Sector-level (`--shard K/N`) — single-machine, batched
Each worker owns sectors with `(chn % N) == K` and builds each owned sector
monolithically via `W1_PW_cache::build` (Hermitian triangle + transpose-fill +
store). Best single-machine throughput because it batches all of a sector's
channels per cell (cross-channel orbital-cache reuse).

### Block-level (`--worker-index I --worker-count N`) — cross-machine, resumable
Worker I owns the evaluate-units whose global evaluate-index satisfies
`idx % N == I`. Each worker batches *its* owned missing channels into one
`integrate_w1_channel_blocks` call (partial cross-channel reuse) and publishes
each shard atomically. Lets N workers (or N machines) attack one large sector
in parallel.

```bash
# create the work plan + manifest
build/bin/w1_worker plan  CPP/Input/input.txt --manifest run/w1_manifest.json

# run 16 independent block-level workers (set OMP so workers×threads <= cores)
for i in $(seq 0 15); do
  OMP_NUM_THREADS=6 build/bin/w1_worker build CPP/Input/input.txt \
    --worker-index "$i" --worker-count 16 \
    cache_root=run/cache output_folder=run/out P123_folder=run/out \
    > run/w$i.log 2>&1 &
done
wait

# inspect completion (exit 0 = complete, 1 = incomplete)
build/bin/w1_worker status CPP/Input/input.txt cache_root=run/cache

# assemble/verify + cryptographic fingerprint
build/bin/w1_worker assemble CPP/Input/input.txt cache_root=run/cache

# run the normal solver using the completed W1 cache (near-zero W1 build time)
./CPP/run CPP/Input/input.txt cache_root=run/cache
```

If half the workers die, rerun them — completed shards are cache hits and are
not recomputed.

### Slurm array example

```bash
#SBATCH --array=0-15
#SBATCH --cpus-per-task=6
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
build/bin/w1_worker build CPP/Input/input.txt \
  --worker-index $SLURM_ARRAY_TASK_ID --worker-count 16 \
  cache_root=/scratch/$USER/w1run/cache \
  output_folder=/scratch/$USER/w1run/out P123_folder=/scratch/$USER/w1run/out
# then a single serial job checks status + assembles:
build/bin/w1_worker status  CPP/Input/input.txt cache_root=/scratch/$USER/w1run/cache
build/bin/w1_worker assemble CPP/Input/input.txt cache_root=/scratch/$USER/w1run/cache
```

No central scheduler or daemon is required; shards are independent files on
shared storage.

---

## 9. Reproducible validation commands

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release && cmake --build build -j

# C++ contract tests (tiny grid, < 2 s each)
ctest --test-dir build -R 'w1_distributed|chiral_n2lo_w1_cache|w1_block_store'

# Sector-level end-to-end acceptance (~3 min)
bash tools/refactor_harness/run_w1_blockdb_acceptance.sh

# Block-level distributed + resume + exactness acceptance (~4 min)
bash tools/refactor_harness/run_w1_block_dist_acceptance.sh
```

---

## 10. Validation evidence (measured)

All measured on the `input_j3nf_multiblock.txt` grid (`Np=4, Nq=3,
two_J_3N_max=3, two_J_3NF_force_max=1, Nangle_3NF=2`, 2 active J=1/2 sectors,
101 blocks/sector = 56 evaluate + 45 transpose_fill each).

| Check | Result |
|-------|--------|
| **Exactness: block-level 16-worker vs monolithic sector-level W1 fingerprint** | **identical** `611ce32acac0a3395178cb6b5aed5eb3bcaaed2d9550d5d555b08eeebb0456df` |
| Exactness: block-level 2-worker vs monolithic | identical (same fingerprint) |
| Exactness: interrupted(7 blocks)+resumed(16 workers) vs monolithic | identical (same fingerprint) |
| Downstream U: block-level cache → solver vs monolithic → solver | **bitwise-identical** both `J^pi = 1/2±` |
| Hermitian contract (assemble reverse-block check) | 90/90 pairs `bad=0` |
| Resume: worker-0 after interruption | `owned=7 built=0 cache_hits=7` (zero recomputation) |
| Cache-hit solver W1 build time | `build=0.0 s` (zero evaluation) |
| C++ unit test `test_w1_distributed` | 0 failures (all 14 contract items) |
| `ctest` (17 tests) | 17/17 passed |

The fingerprint is a streaming SHA-256 over every evaluate-block payload in
deterministic (sector, `a_r`, `a_c`) order. Identical fingerprints ⇔
bitwise-identical `W^(1)`.

---

## 11. Performance measurements

Block-level scaling on one sector (J=1/2+, 56 evaluate blocks), 96-core Xeon
Gold 5318Y, `workers × OMP_THREADS = 96`:

| workers | threads/worker | wall (s) |
|---------|----------------|----------|
| 1       | 96             | 23 |
| 4       | 24             | 22 |
| 16      | 6              | 19 |

The single-machine scaling is weak because the tiny grid is **memory-bandwidth
bound** above ~16 threads — this is the same behavior the existing monolithic
path exhibits (`output/validation/n2lo_3nf_factorized_performance.json`: 16→96
threads gives only 1.12×). Block-level process parallelism cannot beat the
single-machine memory ceiling.

**Where block-level distribution wins** is across **machines** (independent
memory buses) and in **resumability**. Cost model for the realistic candidate
grid (`Np=82, Nq=12`), extrapolated linearly from the measured 101-block sector
(`output/validation/n2lo_3nf_physical_Ay_feasibility.json`):

```
one J^pi sector, 96 threads:      ~19.4 h   (56 evaluate blocks)
across 16 machines (no memory ceiling): ~19.4 h / 16 ≈ 1.2 h
across 64 machines:                     ~19.4 h / 64 ≈ 18 min
```

This is an operational cost model, not a measured large-grid result. The
single-machine block-level path is ~2× slower per sector than the sector-level
batched path (`W1_PW_cache::build`) because it batches only a worker's subset
of channels; for single-machine runs prefer `--shard`, for cross-machine runs
use `--worker-index/--worker-count`.

---

## 12. Convergence-ladder support

Each parameter point gets its own unmistakable cache/manifest identity because
*all* of `Np_WP`, `Nq_WP`, the boundary-array hashes, `J_2N_max`, `two_J_3N_max`,
`two_J_3NF_force_max`, the LECs, the regulator, `Nangle_3NF`, and
`Np/Nq_per_WP_W1` are in the signature hash. The intended sequence
(p-grid → Nq → W1 quadrature → J_2N → J_3N) is supported by reusing lower-J W1
sectors verbatim when raising `two_J_3NF_force_max` (the J ladder reuses
completed sectors — only the single new `J = Jc` sector is built fresh). See
`docs/w1_block_database_workflow.md` for the ladder workflow.

---

## 13. Limitations

- The realistic `Np=82, Nq=12` W1 build has **not** been executed to completion;
  the timings above are a linear cost model from the measured 101-block sector.
- Single-machine block-level scaling is memory-bound; use cross-machine
  distribution or the sector-level `--shard` path for single-machine throughput.
- Static partition (`idx % N`) is adequate when block costs are roughly uniform;
  the channel-zero scan found no exact-zero blocks to prune and costs are
  channel-dependent. A cost-weighted/dynamic-claim queue is a future upgrade
  (the static partition is the safe, lock-free v1).
- transpose_fill blocks are produced by the solver on demand; an `assemble`
  before any solver run reports them as absent (this is not a defect — they are
  generated lazily and verified when both orientations are present).
