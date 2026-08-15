# Sean-compatible 2NF convergence — status report

**Date:** 2026-08-16
**Git SHA:** `7958b0d`
**tmux session:** `tictac-2nf-conv`

## How to inspect

```bash
# Attach
tmux attach -t tictac-2nf-conv

# Inspect without attaching
tmux capture-pane -pt tictac-2nf-conv | tail -100

# Campaign log
tail -20 output/sean_2nf_convergence/logs/campaign_j3n.log

# JSON status
cat output/sean_2nf_convergence/status/status.json | python3 -m json.tool
```

## Phase 3: J_2N_max ladder — COMPLETE

Grid: Np=82, Nq=12, two_J_3N_max=1, three_nucleon_force=none

| rung | Nalpha | wall (s) | P123 time (s) | solve time (s) | peak RSS (MB) |
|------|--------|----------|----------------|-----------------|----------------|
| J2N=1 | 11 | 582 | — | — | — |
| J2N=2 | 19 | 1932 | — | — | — |
| J2N=3 | 27 | 3437 | 0.2 (reuse) | ~3437 | 2423 |
| J2N=4 | 35 | 6507 | 427 (build) | ~6080 | TBD |

### Convergence at Tlab≈10 MeV

| transition | JP=1/2+ max\|dU\| (MeV) | JP=1/2- max\|dU\| (MeV) | reduction | converged? |
|------------|--------------------------|--------------------------|-----------|------------|
| J2N=1→2 | 5.57e-03 | 1.46e-02 | — | no |
| J2N=2→3 | 8.69e-04 | 1.27e-03 | 6-11× | borderline |
| J2N=3→4 | 2.48e-04 | 5.62e-04 | 2-4× | **yes** (<0.001 MeV) |

**Conclusion:** J2N=4 is converged to <0.001 MeV at Tlab≈10 MeV. J2N=3 is borderline (max|dU|≈0.001 MeV).

### Full energy dependence (J2N=3→4)

| Tlab (MeV) | JP=1/2+ max\|dU\| | JP=1/2- max\|dU\| |
|-------------|---------------------|---------------------|
| 0.75 | 3.27e-05 | 2.28e-05 |
| 1.74 | 5.79e-05 | 9.18e-05 |
| 3.35 | 7.63e-05 | 2.12e-04 |
| 5.90 | 9.45e-05 | 3.79e-04 |
| 10.05 | 2.48e-04 | 5.62e-04 |
| 17.18 | 5.78e-04 | 7.04e-04 |
| 30.80 | 1.34e-03 | 8.07e-04 |
| 61.91 | 2.09e-03 | 1.01e-03 |

Convergence is energy-dependent: at Tlab=62 MeV, max|dU|≈0.002 MeV (still small but not <0.001).

## Phase 4: two_J_3N_max ladder — RUNNING

At J2N=3 (borderline) or J2N=4 (converged). Currently running at J2N=3.

| rung | status | notes |
|------|--------|-------|
| 2J=1 | RUNNING (J3N_1) | same as J2N=3, P123 reused |
| 2J=3 | PENDING | adds J=3/2 sectors |
| 2J=5 | PENDING | adds J=5/2 sectors |
| ... | PENDING | |
| 2J=17 | PENDING | Sean's production truncation |

## Phase 8: P123 vs solve timing

**P123 construction is NOT the bottleneck:**
- J2N=4 P123 build: 427 s (7 min)
- J2N=4 solve: ~6080 s (101 min)
- P123 reuse (read HDF5): 0.2 s

P123 files are keyed by (two_J_3N, P_3N, Np, Nq, J_2N_max) and can be reused:
- Across energy sweeps (energy-independent)
- Across two_J_3N_max rungs (per-sector, lower-J sectors unchanged)

## Resource summary

| rung | wall (min) | core-hours | peak RSS (GB) |
|------|------------|------------|----------------|
| J2N=1 | 10 | 8 | ~2 |
| J2N=2 | 32 | 26 | ~2 |
| J2N=3 | 57 | 46 | 2.4 |
| J2N=4 | 108 | 87 | TBD |
| **J2N total** | **207** | **167** | |
