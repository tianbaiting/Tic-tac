# W1 block database: production workflow

This is the production workflow for full-N2LO 3NF calculations after the
J_3NF-truncation + resumable-W1 work. The goal: converge each independent axis
separately, build the expensive exact W1 operator **once**, and reuse it across
solves and energy sweeps.

The expensive object — exact W1 — satisfies

```
time   ~ Nalpha^2 * Np^2 * Nq^2 * Nang^3
memory ~ Nalpha^2 * Np^2 * Nq^2
```

and is **energy-independent at the operator/grid level**. So it is built once
into a permanent, block-addressable database and reused. The Faddeev solve
itself is cheap and unchanged by the 3NF.

## The five independent convergence axes (do NOT conflate)

| Axis | Knob | What it controls |
|---|---|---|
| 2NF pair partial waves | `J_2N_max` | convergence of `V`, `C`, `G` |
| total 3-body J | `two_J_3N_max` | how many J^pi blocks the Faddeev sum runs over |
| **3NF-active J** | **`two_J_3NF_force_max`** | the highest J^pi sector that carries 3NF/W1 |
| W1 angular order | `Nangle_3NF` | exactness of the 3NF angular projection |
| W1 radial quadrature | `Np_per_WP_W1`, `Nq_per_WP_W1` | exactness of the WP bin integration |
| packet grid | `Np_WP`, `Nq_WP` | discretization of the p,q lattice |

`two_J_3NF_force_max` is the new axis: physically the short-range 3NF converges
in J well before the long-range 2NF, so we set

```
K_J = P V + (1+P) W1 ,  J <= J_3NF_max      (W1 built)
    = P V             ,  J >  J_3NF_max      (pure 2NF, no W1)
```

while the Faddeev sum continues to `two_J_3N_max`. Default
`two_J_3NF_force_max = -1` = active in all blocks (legacy full-3NF behaviour).

## Workflow

### Step 1 — converge the 2NF basis/grid
Run 2NF-only (`three_nucleon_force=none`) over increasing `Np_WP`, `Nq_WP`,
`J_2N_max`, `two_J_3N_max` until observables are stable. This fixes the packet
grid and the 2NF partial-wave space. No W1 is involved.

### Step 2 — choose W1 quadrature
Pick `Nangle_3NF` and `Np_per_WP_W1`/`Nq_per_WP_W1`. These are part of the W1
block identity: changing them invalidates the W1 database (by design — the key
encodes them). Convergence is checked by comparing the W1 database at two
orders.

### Step 3 — build W1 sectors incrementally (resumable, distributable)
Use the worker to populate the W1 database without running the solver:

```bash
# What needs to be built (per-sector block counts, active vs 2NF):
build/bin/w1_worker plan   CPP/Input/your_input.txt

# Build one sector:
build/bin/w1_worker build  CPP/Input/your_input.txt --sector 1 1      # J=1/2, +

# Distribute over a cluster (each job one shard; sectors above the cutoff are
# automatically skipped):
build/bin/w1_worker build  CPP/Input/your_input.txt --shard 0/16
build/bin/w1_worker build  CPP/Input/your_input.txt --shard 1/16
# ...
build/bin/w1_worker verify CPP/Input/your_input.txt
```

Properties guaranteed:
- **exact** — no low-rank, no compression, no sparsity thresholds;
- **resumable** — a killed job leaves no partial block (atomic temp+rename
  publication); rerun `build` and completed blocks are cache hits;
- **distributable** — shards are independent; duplicate work is harmless (the
  final publication is atomic);
- **reusable** — the database is keyed by all physics/grid/quadrature inputs, so
  it is reused across energies and across the J ladder.

### Step 4 — run the J_3NF convergence ladder
For `Jc = 1/2, 3/2, 5/2, …` set `two_J_3NF_force_max = 2*Jc` and re-solve.
Lower-J W1 sectors are already in the database (cache hits, zero evaluation);
only the single new `J = Jc` sector is built fresh. So the marginal cost of one
more ladder rung is one sector's W1 build + a cheap solve.

### Step 5 — establish the required 3NF J cutoff
Inspect `ΔO_3NF(Jc) = O[2NF+3NF(J<=Jc)] − O[2NF]` and the incremental
`δ_J O = O(Jc) − O(Jc−1)`. Choose the smallest `Jc` where `δ_J O` is below
tolerance. That is the production `two_J_3NF_force_max`.

### Step 6 — reuse the W1 database for multiple solves/energies
Once the grid/basis/W1-quadrature are fixed, any number of energies (and the
paired 2NF/3NF runs) reuse the same W1 database — the solver reports
`build=0.0 s` for every cached sector. **Contract:** W1 cache reuse across
energies is valid ONLY if `Np_WP, Nq_WP`, the packet boundaries, `J_2N_max`,
the LECs, the regulator, and the W1 quadrature orders are all unchanged. The
SHA-256 grid hashes + full LEC/quadrature fields in the block key enforce this
automatically — a changed grid cannot silently reuse stale blocks.

## Crash safety / correctness contracts
- Block publication is atomic (temp file in the same dir, then `rename()`). A
  killed worker never exposes a partial block at its final path.
- Every read re-validates the key hash and schema; a corrupt or mismatched file
  is treated as a miss, never silently accepted.
- Two parities and all J are handled; sectors above `two_J_3NF_force_max`
  produce **no** W1 blocks (proven by `run_j3nf_cutoff_test.sh`).

## Commands at a glance
```bash
make -C CPP                                                # solver
cmake -S . -B build && cmake --build build -j              # + w1_worker
build/bin/w1_worker plan   CPP/Input/input.txt             # survey
build/bin/w1_worker build  CPP/Input/input.txt --shard 0/8 # build a shard
build/bin/w1_worker verify CPP/Input/input.txt             # load-check
./CPP/run CPP/Input/input.txt two_J_3NF_force_max=3        # solve, 3NF up to J=3/2
```
