# Tic-tac architecture: a reader's map from mathematics to code

This document explains Tic-tac in the same order as the mathematics of the
three-nucleon scattering calculation, so a few-body physicist who has never seen
the code can locate every mathematical object. It is the entry point for
navigation; the physics contracts live in `docs/three_nf_equation_contract.md`
and `docs/n2lo_3nf_conventions.md` and are **not** restated here.

The math flow is:

```
RunConfig
   │
   ▼
PartialWaveBasis + PacketGrid            (the discretized Hilbert space)
   │
   ▼
V ,  P ,  W^(1) ,  G                      (the operators on that space)
   │
   ▼
AGS kernel   K = P V + (1 + P) W^(1)      (single source of truth)
   │
   ▼
packet-space driving matrix  A = C^T K C
   │
   ▼
Neumann series  ->  Pade resummation      (solve  U = A + A G U)
   │
   ▼
on-shell elastic U^{Jpi}
   │
   ▼
scattering amplitude M(theta)             (partial-wave resummation, Python)
   │
   ▼
spin-1 observables  dSigma/dOmega, iT11, T20, T21, T22
```

Each box below names the file(s) that implement it and the type that exposes it
semantically. The labels **[storage]** and **[semantic]** distinguish the
performance-critical raw-array representation from the read-only semantic
interface introduced by the structural refactor; both alias the same memory.

---

## 0. Run configuration

The run is described by a flat `key=value` record parsed into `run_params`
(`include/type_defs.h`). The structural refactor partitions those fields by
physical meaning into a typed `RunConfig` (`src/config/run_config.h`):

| Math / concern | Sub-config | Representative fields |
|----------------|-----------|------------------------|
| interaction model | `PhysicsConfig` | `two_body.potential_model`, `three_body.{three_nucleon_force,c_D,c_E,Lambda_3NF,w1_scale}` |
| basis cutoffs | `BasisTruncation` | `two_J_3N_max`, `J_2N_max` |
| (p,q) packet mesh | `PacketGridConfig` | `Np_WP`, `Nq_WP`, per-axis Chebyshev controls |
| quadrature | `QuadratureConfig` | `Nphi`, `Nx`, `N*_per_WP`, `N*_per_WP_W1`, `Nangle_3NF` |
| solver | `SolverConfig` | `pade_max_order`, `solve_faddeev`, `solve_dense`, `include_breakup_channels` |
| execution | `RuntimeConfig` | `channel_idx`, `parallel_run`, `P123_omp_num_threads` |
| IO | `IoConfig` | `output_folder`, `cache_root`, `energy_input_file` |

`run_params` remains the authoritative record the solver consumes;
`make_run_config()` / `apply_to()` are lossless adapters between the two
(`tests/cpp/test_run_config.cpp` proves a faithful 1:1 round-trip).

Physical constants (`hbarc`, `gA`, `f_pi`, `m_pi`, `Lambda_chi`) are locked in
`include/constants.h`.

---

## 1. The discretized Hilbert space: basis + packet grid

Tic-tac uses three nested bases for the `(alpha, p, q)` discretization
(`docs/three_nf_equation_contract.md` §4.1):

1. **Channel basis** — the partial-wave label
   `alpha = (L_pair, S_pair, J_pair, T_pair, lambda, 2j_spectator; 2J, 2T, parity)`,
   built channel-by-channel in contiguous ranges so each conserved `J^pi` block
   solves independently.
   - **[storage]** `pw_3N_statespace` (`include/type_defs.h`).
   - **[semantic]** `ThreeBodyChannel` + `PartialWaveBasisView`
     (`src/core/state_space/three_body_channel.h`,
     `partial_wave_basis_view.h`): `basis.channel(alpha)`, `basis.block_range(b)`.
   - Builder: `make_pw_symm_states` (`src/core/state_space/make_pw_symm_states.*`).
2. **Free wave-packet (WP) basis** — momentum bins on the Jacobi `p` (pair) and
   `q` (spectator) axes.
   - **[storage]** `fwp_statespace` (`include/type_defs.h`).
   - **[semantic]** `PacketAxisView` + `PacketGridView`
     (`src/core/state_space/packet_grid_view.h`): `grid.p()`, `grid.q()`.
   - Builder: `make_wp_states` (`src/core/state_space/make_wp_states.*`).
3. **Scattering wave-packet (SWP) basis** — the WP basis diagonalizing the pair
   Hamiltonian `H_pair = H0 + V_pair`; the rotation is `C`, the eigenvalues
   `e_SWP`.
   - **[storage]** `swp_statespace` (`C_SWP_*_array`, `e_SWP_*_array`).
   - Builder: `make_swp_states` (`src/core/state_space/make_swp_states.*`).

The pipeline that sequences these is `tictac::app::run_solver`
(`src/app/solver_pipeline.cpp`); `main` is a one-line thunk (`src/main.cpp`).

---

## 2. The operators V, P, W^(1), G

| Math object | What it is | Where |
|-------------|-----------|-------|
| `V` | 2NF pair potential in the WP basis | `V_WP_{unco,coup}_array`, built by `make_potential_matrix` (`src/core/potential/`); model selected via `src/interactions/potential_model.*` |
| `P = P_123 + P_132` | cyclic particle permutation, sparse on the packet lattice | `P123_sparse_{val,row,col}_array`; built by `make_permutation_matrix` (`src/core/state_space/make_permutation_matrix.*`) <br> **[semantic]** `PermutationOperatorView` (`src/core/state_space/permutation_operator_view.h`) |
| `W^(1)` | bare spectator-1 3NF component | abstract `three_nucleon_force_model::W1_element()` (`src/interactions/three_nucleon_force_model.*`); concrete models: `chiral_N2LO_3NF_factorized` (Hebeler-style exact factorized), `chiral_N2LO_3NF_full_reference` (independent 5D reference); cached WP matrix elements via `w1_pw_cache` (`src/interactions/w1_pw_cache.*`) |
| `G` | channel resolvent, diagonal in the SWP basis, `G = (E - e_SWP)^-1` | `G_array`, built by `calculate_resolvent_array_in_SWP_basis` (`src/core/resolvent/make_resolvent.*`) |

The 3NF operator content is the **complete factorized N2LO** set
`c1, c3, c4, cD, cE`; the LEC/regulator/Fourier conventions are locked in
`docs/n2lo_3nf_conventions.md`. `W^(1)` is the *bare* spectator component; the
pair dressing `(1+tG0)` is carried by the channel resolvent `G`, not inserted
into `W1_element` (`three_nf_equation_contract.md` §3.2).

---

## 3. The AGS kernel — the one place K = P V + (1+P) W^(1)

The packet-space driving matrix is

```
A = C^T · K · C ,    K = P·V + (1 + P)·W^(1) ,
```

with `(1+P)` fixed on the **left** of `W^(1)` by Deltuva Eq. (7a)
(`three_nf_equation_contract.md` §3). This algebra has a **single source of
truth**: `src/core/faddeev_solver/cpvc_kernel.{h,cpp}`. Every code path that
claims to compute `A` — the sparse column builder, the on-shell row builder, the
dense solver, and the Neumann/Padé iteration — goes through one of:

- `calculate_PVC_col` — the `P·(V·C)` column.
- `add_one_plus_P_W1_C_col` — the `(1+P)·W^(1)·C` column (used by **both** the
  column and row builders, so the 3NF algebra is never duplicated).
- `calculate_CPVC_col` — the full `A[:,col] = [C^T·(P·V + (1+P)·W^(1))·C][:,col]`.
- `calculate_all_CPVC_rows` — the on-shell row form of the same `A`.

The name `CPVC` is the historical `C^T·P·V·C` (2NF-only) name; with a non-null
3NF context it assembles the full AGS kernel. The module is deliberately free of
HDF5 / 2NF-model / solver-loop dependencies so it can be unit-tested in
isolation. The finite-dimensional operator-level oracle
(`tests/cpp/test_3nf_operator_oracle.cpp`) and the noncommuting operator-order
discriminator (`tests/cpp/test_faddeev_operator_order.cpp`) guard this algebra.

The SWP pair rotation `C` enters as `CT_RM[a*Nalpha+b][i*Np+j]` (row-major `C^T`)
and `VC_CM` (column-major `V·C`); the coupled-channel layout is assembled in
`coupled_channel_transform` (`src/core/faddeev_solver/`).

---

## 4. Solving U = A + A G U: Neumann series and Padé resummation

Once `A` is formed, the AGS equation becomes the matrix iteration

```
a_0 = A ,   a_{n+1} = a_n · G · A ,   U ≈ Σ_{n=0}^{N} a_n ,
```

then a `[N/N]` Padé resums the `a_n` to the on-shell elastic `U`. Both live in
`src/core/faddeev_solver/`:

- the Neumann iteration and on-shell extraction: `solve_faddeev.{h,cpp}`;
- the Padé approximant and convergence honesty: `pade_approximant.{h,cpp}`.

The 3NF enters **only** through `A` (i.e. through `K`); the resolvent `G` is
independent of the 3NF LECs.

---

## 5. From on-shell U to the solver's high-level call

The expensive per-`J^pi`-block solve is exposed as a problem/result API
(`src/core/faddeev_solver/faddeev_solver_facade.{h,cpp}`, Phase 4):

```cpp
tictac::core::FaddeevProblem   problem{ /* P, V, W1, G, bases, on_shell */ };
tictac::core::FaddeevSolveOptions options{ run_parameters, file_identification };
tictac::core::FaddeevResult result = solve_faddeev_block(problem, options);
// result.elastic_U  -> U^{Jpi}_{alpha' alpha} at the on-shell q bins
```

The facade forwards to `solve_faddeev_equations()` with argument order preserved
verbatim, so it is bit-identical to the legacy direct call. On-shell `U` is
written by `store_U_matrix_elements_txt` (`src/io/disk_io_routines.cpp`) as
`U_PW_elements_*.txt`.

> **Public vs. reference headers.** `solve_faddeev.h` mixes the production entry
> (`solve_faddeev_equations`) with brute-force/reference utilities
> (`PVC_col_brute_force`, `CPVC_col_brute_force`, `*_calc_test`) kept for
> diagnostics and the operator oracle. The *production* path is the facade above
> plus `cpvc_kernel`; the `*_brute_force` / `*_calc_test` entries are
> reference/diagnostic only and are not on the hot path.

---

## 6. From U to observables (Python)

The on-shell `U^{Jpi}` becomes physical observables through the partial-wave
reconstruction (not a heuristic reduced-`U` combination), in
`examples/pw_amplitudes.py`:

1. `U^{Jpi} -> M(theta)`: the partial-wave resummation
   `M = -(2pi)^2 mu_pd/q_on · Σ_{Jpi} Σ_{alpha'alpha} U^{Jpi}_{alpha'alpha} G_{alpha'alpha}(theta)`,
   with `G` the `jj`-coupling geometric factor (CG + `Y_{l,m}`).
2. `M(theta) -> observables`: the spin-1 spin-density observables
   `T_{kq}(theta) = Tr[M tau_{kq} M^dagger] / Tr[M M^dagger]`, giving
   `dSigma/dOmega`, `iT11`, `T20`, `T21`, `T22`.

Thin runnable drivers under `examples/` (`run_dpol_p_observables.py`,
`compare_Ay_experiment.py`) call this library; they do not reimplement the
physics.

---

## 7. Cache layer (P123, W1)

The energy-independent, expensive-to-build artifacts are persisted in a
hash-keyed cache (`src/io/cache_layer/`): the permutation matrix `P123` and the
W1 block payloads. Cache identity includes all physics/numerical inputs
(model, LECs, grid, truncation, quadrature, regulator, `Nangle_3NF`) so stale or
cross-convention blocks are rejected by schema versioning.

---

## Where to look for each "success criterion" question

| A physicist looking for... | ...finds it in |
|---------------------------|----------------|
| the basis definition | `make_pw_symm_states.*`, view `PartialWaveBasisView` |
| where `P` is constructed | `make_permutation_matrix.*`, view `PermutationOperatorView` |
| where `V` and `W1` are defined | `make_potential_matrix.*`; `three_nucleon_force_model::W1_element` |
| the one module implementing `K = P V + (1+P) W1` | `cpvc_kernel.{h,cpp}` (`add_one_plus_P_W1_C_col`) |
| the channel resolvent `G` | `make_resolvent.*` (`calculate_resolvent_array_in_SWP_basis`) |
| the Neumann and Padé steps | `solve_faddeev.*`, `pade_approximant.*` |
| where on-shell `U` becomes `M(theta)` | `examples/pw_amplitudes.py` |
| where `M(theta)` becomes `Ay` / tensor observables | `examples/pw_amplitudes.py` (`Tr[M tau_kq M^dagger]/Tr[M M^dagger]`) |

## Regression safety net

The structural refactor is guarded byte-for-byte by
`tools/refactor_harness/run_and_compare.sh` (two deterministic reduced-grid
solves: 2NF-only and complete-factorized-3NF) plus the 15 CTest targets and the
Python `unittest` suites documented in `docs/complete_n2lo_3nf_status.md`.
