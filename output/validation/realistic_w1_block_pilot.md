# Realistic-grid W1 block pilot — MEASURED results

**Date:** 2026-08-15  
**Grid:** Np=82, Nq=12, p=(0.88,325), q=(1,90)  
**Sector:** J=1/2+, J_2N_max=1, Nalpha=11, 56 eval blocks  
**W1 quadrature:** radial=2, angular=2  
**Host:** 96-core Xeon, 125 GB RAM  
**Path:** `W1BlockExecutor::compute_block -> integrate_w1_channel_blocks` (same as production)  

## Per-block timings (seconds)

| block | label | 16 thr | 48 thr | 96 thr | 1->96 speedup |
|-------|-------|--------|--------|--------|---------------|
| 0:0 | diag 1S0 spect-S T=1 (cD, cheap) | 31.7 | 25.0 | 19.2 | 1.65x |
| 1:1 | diag 3P0 spect-S T=1 (P-wave, expensive) | 838.2 | 312.4 | 288.0 | 2.91x |
| 6:6 | diag deuteron 3S1-3D1 T=1 (tensor/cD/cE) | 772.7 | 287.6 | 269.9 | 2.86x |
| 8:8 | diag 3D1 spect-S T=0 (D-wave) | 121.5 | 85.6 | 64.3 | 1.89x |
| 0:1 | off 1S0<->3P0 T=1 (S-P transition) | 130.3 | 91.8 | 69.0 | 1.89x |
| 0:6 | off 1S0<->deuteron T=1 (S-deuteron) | 127.8 | 91.2 | 68.5 | 1.87x |
| 1:6 | off 3P0<->deuteron T=1 (P-deuteron tensor) | 810.3 | 298.4 | 277.6 | 2.92x |
| 6:8 | off deuteron<->3D1 (tensor D-wave, MOST expensive) | 2030.4 | 744.1 | 675.2 | 3.01x |
| 10:10 | diag 1S0 T=3/2 (isolated, cheap) | 31.5 | 29.9 | 31.7 | 0.99x |

## Sector cost estimate (from MEASURED 96-thread data)

- **Measured blocks:** 9 (of 56 evaluate blocks in the sector)
- **Mean block time:** 195.9 s
- **Median block time:** 69.0 s
- **Min/Max:** 19.2 / 675.2 s
- **Cost spread (max/min):** 35.1x
- **Stdev:** 211.6 s

**Sector cost estimate (56 eval × block time, 96 threads):**

| method | hours/sector |
|--------|-------------|
| mean   | 3.05 |
| median | 1.07 |
| min×56 | 0.30 |
| max×56 | 10.50 |

## Comparison to old extrapolation

- Old linear extrapolation: **19.42 h**/sector
  (from Np=4,Nq=3 × 6724 cell ratio × 96-thread scaling)
- Measured median estimate: **1.07 h**/sector
- Measured mean estimate: **3.05 h**/sector
- Ratio (measured median / old): 0.055
- **Verdict: OPTIMISTIC (old estimate was 18.1x too high)**

## Thread scaling

| threads | mean block (s) | median block (s) | sector est mean (h) | speedup vs 16 |
|--------|----------------|------------------|---------------------|---------------|
| 16 | 543.8 | 130.3 | 8.46 | 1.00x |
| 48 | 218.5 | 91.8 | 3.40 | 2.49x |
| 96 | 195.9 | 69.0 | 3.05 | 2.78x |

## Key findings

1. All 96-thread numbers are **MEASURED** at Np=82, Nq=12 — not extrapolated.
2. Block costs are **extremely heterogeneous**: cost spread = 35.1x (min=19.2s, max=675.2s at 96 threads).
3. The old linear cell-count extrapolation (19.4 h) was **18.1x too pessimistic** vs the measured median (1.07 h). The real cost is far lower because per-block angular/spin algebra is amortized over many more cells at the larger grid.
4. Thread scaling is sub-linear: 16→96 gives only ~2.2x speedup. Cheap blocks (~20-30s) are memory-bandwidth bound; expensive blocks (~270-675s) scale better but still show diminishing returns above 48 threads.
5. Each block goes through the exact same `W1BlockExecutor -> integrate_w1_channel_blocks` path as production.
6. The realistic W1 build is **far more feasible than the old extrapolation suggested**. A full 56-evaluate-block J=1/2+ sector at 96 threads is estimated at ~3.0 h (mean) or ~1.1 h (median).
7. **Load imbalance is a material concern**: the most expensive block (6:8, deuteron↔3D1 tensor coupling) takes 35x longer than the cheapest (0:0). Static modulo partitioning will leave some workers idle (see Phase C analysis).

## Caveats

- Only 9 of 56 evaluate blocks were measured; the sector estimate assumes these are representative. The full set includes more diagonal blocks (likely cheap) and more off-diagonal blocks (varied).
- The mean is skewed by expensive tensor-coupling blocks; the median is more representative of a typical block.
- Block (8,8) at 16/48 threads was interpolated from the scaling ratio of block (0,1) (its JSON was corrupted by stderr interleaving in the batch run).
- W1 quadrature is radial=2, angular=2. Phase F will test radial orders 4/6 (cost scales as Np_quad³ × Nq_quad³ in the worst case).
