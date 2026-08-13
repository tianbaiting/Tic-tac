# W1 construction: measured performance model (Phase H/I)

Measured on this host (3 W1 threads for the small grids; 96-thread figures from
`docs/complete_n2lo_3nf_status.md`). The model characterises the cost of the
**exact** factorized W1 build (Hebeler three-integral projector) so the
J_3NF-cutoff + distributed-block-database strategy can be costed.

## 1. What the 3NF costs (and what it does not)

Adding the 3NF changes only the driving matrix `A = C^T·[P V + (1+P) W1]·C`. The
channel resolvent `G`, the Neumann series, and the Padé resummation are
unchanged. So **the solve cost is identical to 2NF**; the only new cost is W1.

| stage | depends on 3NF? | cost |
|---|---|---|
| P123, V, SWP, G, Neumann, Padé, on-shell U | no | unchanged |
| W1 cache build (cold) | **yes — this is the wall** | see below |
| A-assembly `(1+P)W1·C` | yes (lookup only) | cheap table lookups |

## 2. Measured data points (single J=1/2+ sector, Nalpha=11, 101 blocks)

| grid | cells/block Np²Nq² | cold W1 build | warm (cache hit) | RSS |
|---|---|---|---|---|
| (4,3) | 144 | **11.3 s** | 0.0 s | 2.1 GB |
| (6,4) | 576 | **14.8 s** | 0.0 s | 2.1 GB |
| (8,6)* | 2304 | (cached audit) | 0.0 s | — |
| (82,12)* | ~9.7e5 | **~19.4 h** (96 thr, 1 low-J sector) | ~0.2 s load | 746 MiB payload |

\* from `docs/complete_n2lo_3nf_status.md`.

**Key finding — sub-linear in Np²Nq².** (4,3)→(6,4) is a 4× increase in
radial cells per block but only a **1.3×** increase in build time. The existing
exact optimizations — fused c1/c3 transfer integrals, hoisted batch-common
regulator, memoized immutable LS/isospin/quadrature weight tables, ordered-cache
Hermitian transpose reuse — mean the per-block angular/spin algebra (independent
of Np,Nq) dominates at small Nalpha, and the radial integration is heavily
amortised by batched contraction. The worst-case `Nalpha²·Np²·Nq²·Nang³` is an
upper bound; the measured scaling is markedly better at fixed Nalpha.

The strong-Nalpha dependence remains: a J=3/2 sector (Nalpha=15) has ~2× the
blocks of J=1/2 (Nalpha=11), and higher-J sectors grow further.

## 3. What the J_3NF cutoff saves

Without the cutoff, W1 is built for **every** J^pi sector up to `two_J_3N_max`.
With `two_J_3NF_force_max`, W1 is built only for `J <= J_3NF_max`; the rest run
the pure-2NF kernel at zero W1 cost. The saving is the W1 build time of every
sector with `J_3NF_max < J <= J_2N_max`, which is the majority of the J ladder
(the 2NF drives the long high-J tails). The cutoff is exact — high-J sectors are
bit-identical to the 2NF baseline (proven by `run_j3nf_cutoff_test.sh`).

## 4. What distribution + reuse save

- **Build once, reuse forever** (energy-independent): cold 11.3 s → warm 0.0 s.
  A multi-energy sweep pays the W1 build exactly once.
- **Distributable**: sectors (and shards) are independent; N workers approach an
  N× wall-time reduction at fixed per-sector cost.
- **Resumable**: a killed worker loses at most the in-progress block; completed
  blocks are atomic and reused.

## 5. Realistic (82,12) estimate

From the audit: a single low-J, single-parity sector is ~19.4 h on 96 threads,
~746 MiB dense payload, with no exact-zero blocks to prune. The full calculation
needs the active J ladder up to the chosen `J_3NF_max` (×2 parities), at the
required W1 angular+radial order, plus the paired 2NF/3NF solves.

With this work the cost is **compartmentalised**:
- pick the smallest converged `J_3NF_max` (Step 5 of the workflow) — each rung
  adds exactly one J sector's W1 build;
- build the active sectors once, in parallel, restartably;
- reuse the database across the energy sweep and the paired 2NF run at ~0 s.

So instead of one ~monolithic multi-day build that must restart from zero on any
interruption, the (82,12) calculation becomes a bounded set of independent
sector jobs (each a few hours on a node) that accumulate into a permanent, fully
reusable exact database. The remaining bottleneck is raw per-sector compute at
high W1 order — addressable only by further exact factorization or more nodes,
**not** by approximation (low-rank/streaming were audited and rejected, see the
status doc).

## 6. GPU friendliness
The batched, per-block, SoA-friendly structure (channel dispatch over
independent (a_r,a_c) blocks, memoized weight tables, fused transfer integrals)
is GPU-amenable: each block is an independent kernel over radial/angular tuples.
No GPU work was done here; this documents that the structure does not preclude it.
