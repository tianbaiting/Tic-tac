# J_3NF truncation: design note (Phase A)

**Status:** design. Production edits follow in Phase B. Physics contracts in
`docs/three_nf_equation_contract.md` and `docs/n2lo_3nf_conventions.md` are
unchanged by this work — only *which J^π sectors carry a 3NF* becomes
configurable.

## 1. Motivation

Exact full-N2LO W1 construction is the computational wall
(`time ∝ Nα²·Np²·Nq²·Nang³`). Physically the 3NF contribution to elastic nd/pd
observables converges in `J` well before the 2NF partial-wave ladder does: the
long-range 2NF drives the high-J tails, while the 3NF is short-ranged. Building
W1 for every J^π sector up to `two_J_3N_max` is therefore wasted work. We want

```
K_J = P V + (1+P) W1 ,   J ≤ J_3NF_max
    = P V              ,   J >  J_3NF_max
```

while the total Faddeev calculation runs to a larger `two_J_3N_max`. The
production convergence question is

```
ΔO_3NF(Jc) = O[2NF + 3NF active for J ≤ Jc] − O[2NF]
```

established as `Jc = 1/2, 3/2, 5/2, …` **without rebuilding lower-J W1 sectors**.

## 2. Current behaviour (audited)

| Concern | Today | File:line |
|---|---|---|
| pair truncation | `J_2N_max` | `run_params.J_2N_max` |
| 3-body total-J cutoff | `two_J_3N_max` (loop bound over `N_chn_3N`) | `solver_pipeline.cpp` channel loop |
| 3NF activation | one `tnf_ptr_` from `three_nucleon_force_model::create`, passed to **every** solved block | `solver_pipeline.cpp`, `three_nucleon_force_model.cpp:39` |
| W1 build gate | `tnf != nullptr && tnf->enabled() && w1_scale != 0` — **no J term** | `solve_faddeev.cpp:1474` |
| W1 cache key | already encodes `(two_J_3N, P_3N, a_r, a_c)` + all physics/grid/quadrature | `cache_keys.h:24` (`W1Key`) |
| W1 cache write | direct HDF5 write, **not atomic** | `cache_layer.cpp:121` (`store_w1`) |
| observable assembly | resums every emitted `U^{Jpi}` file; mixed-J works | `examples/pw_amplitudes.py` |

**Leak:** the code implicitly assumes "3NF active in every J block that exists".
That assumption lives in exactly one place — the channel loop passing the same
`tnf_ptr_` to each block. Basis construction, the cache key, the kernel algebra,
and observable collection are all already per-block and need no change.

## 3. Desired behaviour / exact operator meaning

Introduce an explicit, independently-named parameter

```
two_J_3NF_force_max   (default -1)
```

Semantics:

- `two_J_3NF_force_max < 0` (default) → 3NF active in **all** solved blocks
  (bit-identical to today; preserves every existing run/output/cache).
- `two_J_3NF_force_max ≥ 0` → for a block with conserved `two_J_3N`,
  3NF is active iff `two_J_3N ≤ two_J_3NF_force_max`.

  - active block: build/load/use W1 normally (`tnf` passed through);
  - inactive block (`two_J_3N > cutoff`): `tnf = nullptr` passed to the solver,
    so `K_J = P V` exactly. **No W1 build, no fake dense zero W1, no allocation.**

This is a strict superset of existing behaviour: `three_nucleon_force=none`
remains a whole-run 2NF; `w1_scale=0` still reduces exactly to 2NF; and
`two_J_3NF_force_max ≥ two_J_3N_max` reproduces full-3NF.

## 4. Why this is safe (no physics change)

- The `tnf == nullptr` path is **already** the audited, bitwise-verified 2NF path
  (the `three_nucleon_force=none` baseline). Routing high-J blocks through it
  changes no arithmetic — it selects a different value of `K_J` that is already a
  supported, contract-compliant kernel.
- The AGS kernel module (`cpvc_kernel.cpp`) is the single source of truth for
  `A = C^T·(P·V + (1+P)·W1)·C`. When `tnf==nullptr`, `add_one_plus_P_W1_C_col`
  is never called and the kernel reduces to `A = C^T·P·V·C`. No new algebra.
- `G` (channel resolvent) is 2-body and 3NF-independent regardless; the
  Neumann/Padé solve is unchanged.

## 5. Cache consequences

None required. The `W1Key` already partitions by `(two_J_3N, P_3N, a_r, a_c)`:

- An inactive high-J sector simply never calls `W1_PW_cache::build`, hence never
  emits/loads W1 cache entries. **No high-J W1 file is created.**
- Lower-J sectors built under one `J_3NF_max` are reused verbatim under a larger
  `J_3NF_max` (the J ladder reuses completed sectors without rebuild).
- Cross-sector contamination is impossible: two sectors differ in `two_J_3N` or
  `P_3N`, hence in the key, hence in distinct files. (Phase B test 6 proves this.)
- `two_J_3NF_force_max` is **not** part of `W1Key`: a block's W1 value depends
  only on the block's own physics, never on which other sectors happen to be
  active. This is correct and deliberate.

## 6. Implementation surface (Phase B)

Minimal:

1. Add `int two_J_3NF_force_max` to `run_params` (default -1) and to the semantic
   `RunConfig` (`PhysicsConfig.three_body`), with parsing + adapters.
2. In the channel loop (`solver_pipeline.cpp`), gate the per-block `tnf`:
   ```
   const bool active = (two_J_3N <= rp.two_J_3NF_force_max) || (rp.two_J_3NF_force_max < 0);
   problem.tnf = (active && tnf_ptr_) ? tnf_ptr_.get() : nullptr;
   ```
   (`two_J_3N` is already computed there.)
3. Tests (counting/fault-injection) proving high-J blocks make zero W1 calls.

No change to: `cpvc_kernel`, `W1Key`, `W1_PW_cache::build`, `solve_faddeev`'s W1
gate, the Padé path, or observable reconstruction.

## 7. Convergence workflow (Phase C)

For fixed `Np, Nq, J_2N_max, Nangle_3NF, W1 radial order, LECs, regulator`:

1. build the 2NF baseline once (all sectors 2NF);
2. for `Jc = 1/2, 3/2, …`: set `two_J_3NF_force_max = 2*Jc`, re-solve — lower-J
   W1 sectors are cache hits, only the new `J = Jc` sector is built fresh;
3. report `O_2NF`, `O_2NF+3NF(Jc)`, `ΔO_3NF(Jc)`, and `δ_J O`.

Results are labelled **diagnostic** until the Np/Nq/W1-quadrature/Padé gates are
met; this workflow isolates the *J_3NF* convergence axis, which is independent of
the Np/Nq / W1-quadrature / J_2N axes.
