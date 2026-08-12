# Tic-tac

forked from https://github.com/seanbsm/Tic-tac. and reconstructed. 

Tic-tac is a three-nucleon Faddeev-equation solver (wave-packet discretization, WPCD) for nucleon-deuteron scattering.
This repository includes the core solver, input/config workflow, and a validated Tlab-aligned dpol-p comparison pipeline for the benchmark often labeled 190 MeV/u.

## Scope

- Solve elastic nd/pd amplitudes (`U_PW_elements_*`).
- Build WP/SWP state spaces and solve AGS/Faddeev equations.
- Validate the maintained `Tlab = 190 MeV` polarized deuteron-proton benchmark against experimental data in `data/DataOfCrosssectionAndPol`.

## Repository Layout

- `src/`: main C++ source (core solver, state space, resolvent, potential, IO, config).
- `CPP/`: legacy-compatible build/run tree used by `CPP/run` workflow.
- `examples/`: end-to-end scripts for solver run, validation, and plotting.
- `data/`: experimental data and auxiliary inputs.
- `docs/`: validation and method notes.
- `tests/`: Python and legacy numeric tests.

## Dependencies

Minimum runtime/build dependencies:

- C++ compiler with C++17 support
- gfortran
- BLAS/LAPACK
- GSL
- HDF5
- OpenMP
- Python 3

Optional for plotting:

- matplotlib (recommended via micromamba env `anaroot-env` in this workspace)

## Build

There are two build entrypoints over the same `src/` + `include/` tree:

- **`make`** (repository root) builds the executable `./tic-tac`. This is the
  canonical CMake-parity build target.
- **`make -C CPP`** builds `./CPP/run`. This is the executable the maintained
  Python workflows (`examples/*.py`) and the regression tests expect, so it is
  the recommended target for day-to-day work:

```bash
make -C CPP
```

The two are the same sources compiled with the same flags; `CPP/makefile` is a
thin wrapper that delegates to the root `Makefile` with `TARGET=CPP/run`. Both
honour the same `TICTAC_USE_NEW_CACHE_LAYER` switch (default `1`, matching the
CMake `option(... ON)`); pass `TICTAC_USE_NEW_CACHE_LAYER=0` to build the legacy
no-cache path.

When using the `anaroot-env` conda environment for HDF5/GSL/BLAS, point the
Makefile at it with `CONDA_PREFIX` (the conda BLAS additionally requires
`-lcblas`, which the Makefile adds automatically in that case):

```bash
CONDA_PREFIX=/data/tian/conda/envs/anaroot-env make -C CPP -j
```

If you need command help from solver:

```bash
./CPP/run -h
```

## Input Profiles

Canonical hand-written input files are kept in `CPP/Input/`:

- `CPP/Input/input.txt`: baseline profile
- `CPP/Input/input_Ay_test.txt`: lower-cost `Tlab = 190 MeV` quick test profile
- `CPP/Input/input_Ay_fixed.txt`: refined `Tlab = 190 MeV` profile
- `CPP/Input/input_Ay_nijmegen_legacy.txt`: legacy Nijmegen profile (kept for reproducibility)

Run from repository root:

```bash
./CPP/run CPP/Input/input.txt
```

## Tlab = 190 MeV dpol-p Workflow

External references often call this benchmark `190 MeV/u`. In this repository, solver-facing inputs are always `Tlab [MeV]` to match the C++ core.

### 1. Run solver

```bash
python3 examples/deuteron_proton_Ay.py --work-dir output/deuteron_proton_Ay --target-tlab-mev 190 --reuse-p123
```

### 2. Validate against experiment

```bash
python3 examples/compare_Ay_experiment.py --work-dir output/deuteron_proton_Ay --solver-out-dir output/deuteron_proton_Ay/solver_out --target-tlab-mev 190
```

`dSigma/dOmega` output unit can be selected with `--dsigma-unit`:
- `mb/sr` (default)
- `fm2/sr`

### 3. Generate comparison figures (matplotlib)

```bash
micromamba run -n anaroot-env python examples/plot_validation_curves.py --work-dir output/deuteron_proton_Ay
```

## Multi-Tlab dpol-p Observables (70/135/190 MeV)

Run solver-driven (Faddeev `U` -> observables) pipeline:

```bash
python3 examples/run_dpol_p_observables.py --work-dir output/dpol_p_observables --target-tlabs-mev 70,135,190
```

Plot all observables (`dSigma/dOmega`, `iT11`, `T20`, `T21`, `T22`) with matplotlib:

```bash
micromamba run -n anaroot-env python examples/plot_dpol_p_observables.py --work-dir output/dpol_p_observables
```

To generate curves in `fm^2/sr` instead of `mb/sr`:

```bash
python3 examples/run_dpol_p_observables.py --work-dir output/dpol_p_observables --target-tlabs-mev 70,135,190 --dsigma-unit fm2/sr
```

Detailed algorithm and output hierarchy:
`docs/dpol_p_multi_energy_observables.md`

Tic-tac tunable parameters and practical tuning order:
`docs/tictac_parameter_tuning.md`

## Key Outputs

Generated under `output/deuteron_proton_Ay/`:

- `solver_validation_tlab_190MeV.txt`
- `solver_validation_tlab_190MeV.json`
- `best_energy_iT11_curve.csv`
- `best_energy_dsigma_curve.csv`
- `best_energy_exp_vs_faddeev_annotated.png` (single annotated comparison figure)
- `best_energy_iT11_comparison.png`
- `best_energy_dsigma_comparison.png`

## dSigma Unit Conversion

Experimental `DSigamaOverDOmega.txt` is in `mb/sr`.

- `1 mb/sr = 0.1 fm^2/sr`
- `1 fm^2/sr = 10 mb/sr`

Convert file units reproducibly:

```bash
python3 examples/convert_dsigma_units.py --to-unit fm2/sr
```

## How Experiment vs Simulation Is Computed

- Experiment data:
  - `data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt` for `iT11(theta)`
  - `data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt` for `dSigma/dOmega(theta)`
- Simulation data:
  - Parse `U00,U01,U10,U11` from solver `U_PW_elements_*.txt` (all available `JP` channels from the latest run family)
  - Combine channels at the same energy with `|U|^2` weights (across `JP` and parity)
  - Build angle-dependent observables from reduced-U invariants with fixed formulas
  - Do not fit model coefficients to experimental curves
  - Report MAE/RMSE/max error and relative RMSE

Details and formulas: `docs/dpol_p_190MeV_validation.md`

## Tests

```bash
python3 -m unittest tests/test_190mev_data_pipeline.py
```

Legacy/extended tests exist in subfolders and may require separate build steps.

## Notes

- Input `Tlab` is mapped to solver WP on-shell bins; output energy points can differ from the exact input energy.
- For new output directories, avoid `--reuse-p123` until P123 files are generated once.

## License and Citation

- License: `LICENSE` (MIT)
- Reference publication: https://journals.aps.org/prc/abstract/10.1103/PhysRevC.106.024001
