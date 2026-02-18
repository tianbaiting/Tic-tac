# Repository Guidelines

## Project Structure & Module Organization
- `src/` is the primary codebase: `core/` (solver/state-space/resolvent/potential), `interactions/` (nuclear potentials, including Fortran bridges), `config/`, `io/`, and `utils/`.
- `include/` contains shared headers (`constants.h`, `type_defs.h`).
- `data/` stores runtime inputs and reference datasets; default runtime input is `data/input.txt`.
- `examples/` contains runnable workflows (`run_examples.sh`, `deuteron_proton_Ay.py`, `quick_Ay_test.py`).
- `tests/` and `Test/` contain legacy numerical test harnesses (`Cont_Faddeev`, `Free_energy`).
- `CPP/` is an older parallel implementation; prefer changes in `src/` unless intentionally maintaining both trees.

## Build, Test, and Development Commands
- `make` builds the main executable `./tic-tac` from `src/`.
- `./build.sh release` runs CMake + parallel build in `build/` and copies `bin/tic-tac` to project root.
- `./build.sh debug` produces a debug build.
- `./tic-tac --input data/input.txt` runs a configured calculation.
- `python3 config.py save data/input.txt` regenerates a baseline input file.
- `./examples/run_examples.sh test` runs a fast smoke test configuration.
- `cd tests/Cont_Faddeev && make && ./run` (and similarly in `tests/Free_energy`) runs legacy numerical checks.

## Coding Style & Naming Conventions
- Use C++17-compatible code (current build config) and preserve mixed C++/Fortran interoperability.
- Match local file style: existing code uses braces on the same line and mixed tabs/spaces in legacy files.
- Prefer descriptive `snake_case` for variables/functions and keep file pairs aligned (e.g., `foo.h` + `foo.cpp`).
- No repo formatter is configured (`.clang-format` absent); keep formatting changes minimal and local.

## Testing Guidelines
- There is no reliable top-level `make test`; use example smoke tests plus relevant `tests/` harnesses.
- For solver or potential changes, run at least one `examples/run_examples.sh test` case and one domain-specific harness.
- For physics-facing changes, record key output deltas (energies, amplitudes, Ay curves) in the PR.

## Commit & Pull Request Guidelines
- Follow existing history: short, imperative commit subjects (e.g., `core: improve Neumann convergence printout`).
- Keep commits focused to one technical change.
- PRs should include: problem statement, files touched, exact reproduce commands, and numerical impact.
- Link related issues and include plots/tables when analysis scripts or observable outputs change.
