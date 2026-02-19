# Tic-tac Functionality

This file summarizes the maintained capabilities in this repository.

## Solver Core
- Three-body Faddeev/AGS solving with wave-packet discretization (WPCD).
- Partial-wave basis construction for 3N channels.
- Free/scattering wave-packet basis generation.
- Sparse permutation matrix (`P123`) construction and HDF5 storage.
- Iterative (Padé/Neumann) and dense debug solve modes.

Main code paths:
- `src/core/faddeev_solver/`
- `src/core/state_space/`
- `src/core/potential/`
- `src/core/resolvent/`
- `src/config/`, `src/io/`

## Potentials and Physics Controls
- Potential models: `LO_internal`, `N2LOopt`, `Idaho_N3LO`, `nijmegen`, `malfliet_tjon`.
- Switches for tensor-force coupling and `1S0` isospin breaking.
- Configurable angular-momentum truncation and momentum-grid discretization.

Parameter reference:
- `docs/tictac_parameter_tuning.md`

## Validated Workflows
- 190 MeV/u dpol-p solver output + experiment comparison:
  - `examples/deuteron_proton_Ay.py`
  - `examples/compare_Ay_experiment.py`
  - `examples/plot_validation_curves.py`
- Multi-energy observables (70/135/190 MeV/u):
  - `examples/run_dpol_p_observables.py`
  - `examples/plot_dpol_p_observables.py`

Workflow documentation:
- `docs/dpol_p_190MeV_validation.md`
- `docs/dpol_p_multi_energy_observables.md`
- `docs/algorithm_flow_and_logic.md`

## Build/Run Entry Points
Primary maintained runtime path:
- `cd CPP && make -j`
- `./CPP/run CPP/Input/input.txt`

Automated check:
- `python3 -m unittest tests/test_190mev_data_pipeline.py`

## Notes on Legacy Content
Legacy experiments and one-off exploratory scripts are intentionally minimized.
Use the `examples/` scripts above for reproducible outputs and documentation-aligned behavior.
