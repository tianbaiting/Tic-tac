# Sean-compatible 2NF convergence campaign — live status

**tmux session:** `tictac-2nf-conv`
**git SHA:** `dfaf583`
**hostname:** `ENPG`
**OMP threads:** 48
**Started:** 2026-08-15 22:08

## How to inspect

```bash
# Attach to tmux session
tmux attach -t tictac-2nf-conv

# Inspect without attaching
tmux capture-pane -pt tictac-2nf-conv | tail -100

# Check campaign log
tail -20 output/sean_2nf_convergence/logs/campaign.log

# Check JSON status
cat output/sean_2nf_convergence/status/status.json | python3 -m json.tool
```

## J_2N_max ladder (Np=82, Nq=12, two_J_3N_max=1, pure 2NF)

| rung | Nalpha | status | wall (s) | P123 time (s) | solve time (s) | notes |
|------|--------|--------|----------|----------------|-----------------|-------|
| J2N=1 | 11 | DONE (reused) | 582 | — | — | from previous campaign |
| J2N=2 | 19 | DONE (reused) | 1932 | — | — | from previous campaign |
| J2N=3 | 27 | DONE | 3437 | 0.2 (reuse) | ~3437 | P123 reused, both parities solved |
| J2N=4 | 35 | RUNNING | — | 427 (build) | — | P123 built, solve in progress |

### Convergence at Tlab≈10 MeV

| transition | JP=1/2+ max\|dU\| (MeV) | JP=1/2- max\|dU\| (MeV) | reduction factor |
|------------|--------------------------|--------------------------|------------------|
| J2N=1→2 | 5.57e-03 | 1.46e-02 | — |
| J2N=2→3 | 8.69e-04 | 1.27e-03 | ~11x (1/2-), ~6x (1/2+) |

The J2N=2→3 correction is ~10x smaller than J2N=1→2. If J2N=3→4 shows another ~10x reduction,
J2N=3 will be converged to ~1e-4 MeV — sufficient for the campaign.

### P123 timing (Phase 8 analysis)

- P123 construction is NOT the bottleneck: J2N=4 P123 build = 427s vs solve ~5700s (est)
- P123 reuse via `calculate_and_store_P123=false` works perfectly (0.2s read vs 427s build)
- P123 files are keyed by (two_J_3N, P_3N, Np, Nq, J_2N_max) — reusable across two_J_3N_max rungs

## two_J_3N_max ladder (pending J2N convergence)

| rung | status |
|------|--------|
| 2J=1 | (reused from J2N ladder) |
| 2J=3 | PENDING |
| 2J=5 | PENDING |
| ... | PENDING |
| 2J=17 | PENDING |

## Nq ladder (deferred)

Pending until angular axes are settled.
