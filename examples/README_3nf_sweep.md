# 3NF sweep — reproducing Miller-style dσ/dΩ and iT11 comparisons

This workflow runs the Tic-tac solver in three 3NF configurations side-by-side
at a single target lab energy, then overlays the resulting observable curves
against experimental data.

## Files

| File | Purpose |
|------|---------|
| `run_3nf_sweep.py`     | Drives CPP/run three times: `no_3nf`, `zero_lec`, `witala`. Shares P123 cache. |
| `plot_3nf_effect.py`   | Reads the three solver_out directories, assembles M(theta) via the full partial-wave summation in `pw_amplitudes.py`, and writes overlay SVG + CSV + JSON summary. |
| `fetch_nd_exfor.py`    | Reusable EXFOR fetch/parse utility for grabbing n+d elastic entries by ID. |

## The three configurations

| Name | `three_nucleon_force` | `c_D` | `c_E` | Notes |
|------|-----------------------|-------|-------|-------|
| `no_3nf`   | `none`        | –     | –     | pure 2NF baseline (bit-for-bit regression target) |
| `zero_lec` | `chiral_N2LO` | 0.0   | 0.0   | 3NF kernel path active; only 2PE c₁,c₃ contribute when `potential_model ∈ {N2LOopt, Idaho_N3LO}` |
| `witala`   | `chiral_N2LO` | -0.20 | -0.205| Witala et al. PRC 77 (2008) 034004 typical LECs |

The 2PE LECs c₁, c₃, c₄ are inherited from the chosen 2NF potential model
(see `src/interactions/three_nucleon_force_model.cpp::fetch`). With
`potential_model=LO_internal` they are zero and the `zero_lec` config
reduces to the 2NF baseline — a useful sanity check.

## Quick runs

```bash
# small grid smoke test at 64.5 MeV (Witala 1998 "smoking gun" energy)
python3 examples/run_3nf_sweep.py \
    --work-dir output/3nf_sweep_64p5 \
    --target-tlab-mev 64.5 \
    --potential-model N2LOopt \
    --np 20 --nq 20 --nphi 24 --nx 24 --np-per-wp 6 --nq-per-wp 6 \
    --threads 4 --timeout 3600

python3 examples/plot_3nf_effect.py \
    --work-dir output/3nf_sweep_64p5 \
    --target-tlab-mev 64.5

# canonical 190 MeV/u comparison (same grid as deuteron_proton_Ay.py)
python3 examples/run_3nf_sweep.py \
    --work-dir output/3nf_sweep_190 \
    --target-tlab-mev 190.0 \
    --potential-model N2LOopt

python3 examples/plot_3nf_effect.py \
    --work-dir output/3nf_sweep_190 \
    --target-tlab-mev 190.0
```

## Expectations and caveats

- **Small 3NF signal for zero_lec vs no_3nf** — with c_D=c_E=0, only the 2PE
  c_i channel contributes through `W1_2pe`. That contribution uses the
  monopole/scalar approximation (tensor parts deferred), so the shift will
  under-predict the full 3NF effect.

- **Observable assembly** — `compare_Ay_experiment.py` and `plot_3nf_effect.py`
  now assemble `M(theta)` by the full partial-wave summation in
  `pw_amplitudes.py` (Clebsch–Gordan + spherical-harmonic kernel, jj coupling,
  Witala/Gloeckle WP-basis prefactor). The dσ/dΩ absolute scale carries one
  empirical multiplicative factor (geometric-mean log fit to experiment) that
  absorbs residual WP-basis normalization; polarization observables iT11,
  T20, T21, T22 are unaffected by this and stay in the [-1, 1] physical band
  without clamping.

- **190 MeV data** — `data/DataOfCrosssectionAndPol/` is d+p (with Coulomb),
  whereas the solver here is nd. The Coulomb-free solver prediction is
  typically close at backward angles but diverges at forward angles.

## Fetching experimental data from EXFOR

```bash
# Download one or more EXFOR entries (raw + parsed) into data/nd_scattering/
python3 examples/fetch_nd_exfor.py --entries 10135 \
    --reaction "1-H-2(N,EL)" \
    --out-dir data/nd_scattering/
```

The fetcher pulls from the IAEA-NDS master archive at
`nds.iaea.org/nrdc/exfor-master/entry/<prefix>/<entry>.txt`, stores raw files
under `data/nd_scattering/_raw_cache/`, and writes per-subentry parsed tables
`exfor_<ENTRY>_<SUBENT>.tsv` plus `manifest.json`. It handles fixed-width
11-char columns and the EXFOR `1.234+3` → `1.234E+3` float convention.

### Verified n+d elastic entries

| Entry | Energies | Quantity | Reference |
|-------|----------|----------|-----------|
| 10135 | 36.0, 46.3 MeV | dσ/dΩ (subent 002, 003) | Romero et al., Phys. Rev. C 2 (1970) 2134 |

The IAEA search UI is JavaScript-heavy and blocked behind CloudFlare challenges
for scripted access, so entry IDs generally need to come from a paper's own
EXFOR cross-reference or from a manual search at
https://nds.iaea.org/exfor/ (reaction `1-H-2(N,EL)`, quantity `DA` / `POL/ANA`).
Additional known-useful references (entry IDs not yet confirmed by this
fetcher):

| Energy | Reference | Quantity |
|--------|-----------|----------|
| 67 MeV | R. Rühl et al., Nucl. Phys. A 524, 377 (1991) | Ay |
| 135–190 MeV | K. Ermisch et al., PRC 68, 051001 (2003) | dσ/dΩ, Ay, iT11 |
| 135–250 MeV | A. Ramazani-Moghaddam-Arani et al., arXiv:1211.5243 | dσ/dΩ |

The 64.5 MeV "smoking gun" energy (Witała et al., Phys. Rev. Lett. 81, 1998)
is not currently available as raw nd data in EXFOR under our reaction filter;
the Romero 46.3 MeV entry is the closest nd dataset for the default smoke-test
Tlab. Higher-energy comparisons should use the 190 MeV/u dpol-p benchmark in
`data/DataOfCrosssectionAndPol/` (note: d+p with Coulomb, not nd).
