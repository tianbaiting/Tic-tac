# Tic-tac

Tic-tac is a three-nucleon Faddeev-equation solver (wave-packet discretization, WPCD) for nucleon-deuteron scattering.
This repository includes the core solver, input/config workflow, and a validated 190 MeV/u dpol-p comparison pipeline.

## Scope

- Solve elastic nd/pd amplitudes (`U_PW_elements_*`).
- Build WP/SWP state spaces and solve AGS/Faddeev equations.
- Validate 190 MeV/u polarized deuteron-proton observables against experimental data in `data/DataOfCrosssectionAndPol`.

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

```bash
make
```

If you need command help from solver:

```bash
./CPP/run -h
```

## 190 MeV/u dpol-p Workflow

### 1. Run solver

```bash
python3 examples/deuteron_proton_Ay.py --work-dir output/deuteron_proton_Ay --target-tlab 190 --reuse-p123
```

### 2. Validate against experiment

```bash
python3 examples/compare_Ay_experiment.py --work-dir output/deuteron_proton_Ay --solver-out-dir output/deuteron_proton_Ay/solver_out --target-tlab 190
```

### 3. Generate comparison figures (matplotlib)

```bash
micromamba run -n anaroot-env python examples/plot_validation_curves.py --work-dir output/deuteron_proton_Ay
```

## Multi-Energy dpol-p Observables (70/135/190 MeV/u)

Run solver-driven (Faddeev `U` -> observables) pipeline:

```bash
python3 examples/run_dpol_p_observables.py --work-dir output/dpol_p_observables --energies 70,135,190
```

Plot all observables (`dSigma/dOmega`, `iT11`, `T20`, `T21`, `T22`) with matplotlib:

```bash
micromamba run -n anaroot-env python examples/plot_dpol_p_observables.py --work-dir output/dpol_p_observables
```

Detailed algorithm and output hierarchy:
`docs/dpol_p_multi_energy_observables.md`

## Key Outputs

Generated under `output/deuteron_proton_Ay/`:

- `solver_validation_190MeV.txt`
- `solver_validation_190MeV.json`
- `best_energy_iT11_curve.csv`
- `best_energy_dsigma_curve.csv`
- `best_energy_exp_vs_faddeev_annotated.png` (single annotated comparison figure)
- `best_energy_iT11_comparison.png`
- `best_energy_dsigma_comparison.png`

## How Experiment vs Simulation Is Computed

- Experiment data:
  - `data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt` for `iT11(theta)`
  - `data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt` for `dSigma/dOmega(theta)`
- Simulation data:
  - Parse `U00,U01,U10,U11` from solver `U_PW_elements_*.txt`
  - Combine parity channels at the same energy with `|U|^2` weights
  - Build angle-dependent observables with Legendre/tanh parameterization
  - Fit coefficients by ridge regression (pure Python normal equations + Gaussian elimination)
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
