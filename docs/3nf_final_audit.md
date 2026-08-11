# Final 3NF physics audit — `fix/3nf-physics-contract`

**Audit date:** 2026-08-11

**Audited code through:** `2e53680`
**Scope:** AGS/WPCD kernel algebra, coupled-channel basis rotations, W1 cache,
enabled N2LO-3NF pieces, WP integration, dense/Padé solution paths, and small
numerical runs.

This report supersedes `docs/3nf_physics_audit_report.md` where the two differ.
It deliberately separates demonstrated statements from approximations and
uncompleted validation. Passing tests is evidence for the property tested; it
is not evidence that the current interaction is a complete N2LO 3NF or that a
physical neutron--deuteron calculation is converged.

## Physics-status matrix

| Item | Status | Independent evidence and scope |
|---|---|---|
| AGS operator structure | **VERIFIED** | Direct compound-index derivation and noncommuting dense-matrix oracle give `A=C^T[PV+W1(1+P)]C`; the oracle rejects `(1+P)W1`. The channel resolvent, rather than `G0`, is used. Literature traceability is recorded in `three_nf_equation_contract.md`. |
| `W1*C` coupled contraction | **VERIFIED** | Regression first failed the old code with 64 mismatches (max error 1.378), then passed after adding `sum_alpha_j W1_(r,j) C_(j,c)`. The present operator oracle has 26 checks and zero failures. |
| `W1*P*C` contraction | **VERIFIED** | Independent dense assembly with nonsymmetric, alpha-off-diagonal `C`, `W1`, and noncommuting `P`; production row and column paths both agree. `P=2P123` is exercised only in the full permutation test, not in the isolated `W1*C` discriminator. |
| `C`/`C^T` convention | **VERIFIED** | Pointer-layout oracle and a nonsymmetric real 3S1--3D1 coupled-block restructuring test verify the transpose and explicitly reject the untransposed orientation. |
| 2NF-only and `W1=0` | **VERIFIED** | Python regression compares 40 complex U elements; current output equals the validated baseline. A separate end-to-end run gives byte-identical U and Neumann files for `three_nucleon_force=none` and `w1_scale=0`. |
| Operator order `W1(1+P)` | **VERIFIED** | `test_faddeev_operator_order` uses `[W1,P] != 0`; five checks pass and the wrong operator differs by a nonzero norm. |
| W1 cache equivalence | **VERIFIED** | Production cache consumer equals direct `W1_element` fallback for every toy element, including intermediate `alpha_j`. Cache keys include row/column alpha block, all LECs, cutoff/regulator, model, grid hashes, channel metadata, and W1 quadrature orders. |
| Row/column kernel builders | **VERIFIED** | `calculate_CPVC_col` and `calculate_all_CPVC_rows` call the same W1 helper and agree with each other and independent dense algebra for every selected toy element. |
| Dense vs Padé/Neumann | **VERIFIED, bounded scope** | Production Padé `[2/2]` agrees with a two-alpha dense toy solve to `1e-11`. On the Np=Nq=5 2NF grid, the selected truncations agree with dense to max `6.57420273e-9 MeV`, mean `4.82005695e-10 MeV`. This does not certify the final Padé tail: all those elements remain honestly labelled `Conv=2`. |
| Padé convergence reporting | **VERIFIED** | The previous code had order-minus-one out-of-bounds reads and could call a moving sequence converged at `P[0/0]`. The solver now evaluates the full configured history and requires the final three adjacent orders to satisfy `1e-7 + 1e-5*scale`; otherwise it reports max-order truncation. Ten focused checks pass. |
| c_E | **VERIFIED within repository conventions** | Independent Pauli enumeration gives `tau2.tau3=-3,+1`; the spin-scalar 3S1/1S0 ratio is `-3`. Closed-form and production matrix elements agree, and the independent angular oracle gives production/oracle `1.0000`. |
| c_D | **VERIFIED only for enabled rank-0 spectator S-wave sector** | Independent eight-state spin/isospin m-scheme enumeration, analytic angular integral, three production/oracle momentum points (`1.00000000`), reverse-matrix Hermiticity, and selection rules agree. Higher-l and rank-2 pieces are fail-closed and therefore **NOT VERIFIED/NOT IMPLEMENTED**. |
| c1/c3 2PE | **APPROXIMATE** | Production is restricted to rank-0 monopole/azimuthal approximation. A full-vector angular oracle finds pointwise production/full ratios 0.5821--1.3407 (3.9--41.8% absolute relative error over the six-point scan). Rank-2 is fail-closed. This is not an exact PWD. |
| c4 | **NOT IMPLEMENTED; fail-closed** | The full model name is rejected, nonzero c4 is rejected by the exact class, and the approximation explicitly reports the inherited c4 drop. The isospin Pauli structure is tested, but the spin--momentum cross-product PWD is absent. |
| W1 WP integration | **CONDITIONAL** | Default per-cell order is now N=2, not midpoint N=1. In the cell oracle, normal c1/c3 bins have N2/N8 error `1.16e-4`, while wide bins have `1.01e-2`; N4/N8 is `1.28e-5`. On a coarse end-to-end Np=Nq=5 run, N2/N4 changes U by up to 28.6%, while N4/N8 changes it by `3.67e-4` relative. Per-run N2/N4 (and when needed N8) convergence is mandatory. |
| Real tensor-coupled channel | **VERIFIED for transform/kernel indexing** | The 3S1--3D1 block test uses the production coupled-matrix restructuring with a deliberately nonsymmetric matrix; the full operator oracle uses genuine alpha-off-diagonal S/D-like blocks. This is not a verification of a physical c1/c3 or cD tensor PWD, because those unverified pieces are disabled. |
| Physical nd benchmark | **NOT VERIFIED** | The Np=Nq=5 diagnostics have an unphysical `E_bound=-0.020026 MeV` and are algebra/cache/solver checks only. The quick 190-MeV Ay workflow selected a 161.145-MeV grid point and failed its observable comparison. No converged comparison with Witała/Miller data was established in this task. |

## Verified correct

### Kernel and index contract

The implemented equation is

```text
U = K + K G U,
K = P V + W^(1) (1+P),
A = C^T [P V + W^(1) + W^(1)P] C.
```

`G` is the pair-channel resolvent constructed in the SWP basis. `W^(1)` is the
bare spectator-1 3NF component; no extra `(1+tG0)` is inserted. The 3NF remains
to the left of `(1+P)`.

For compound pair indices, the identity contribution is

```text
(C^T W1 C)_(alpha_r p_r q_r, alpha_c p_c q_c)
 = sum_(alpha_x,p_x) (C^T)_(alpha_r p_r,alpha_x p_x)
   sum_(alpha_j,p_j) W1_(alpha_x p_x q_r,alpha_j p_j q_c)
                     C_(alpha_j p_j,alpha_c p_c).
```

The actual pointer-table convention is

```text
CT_RM_array[a*Nalpha+b][i*Np+j]
    = (C^T)_(a i,b j)
    = C_(b j,a i),

C_(alpha_j p_j,alpha_c p_c)
    = CT_RM_array[alpha_c*Nalpha+alpha_j][p_c*Np+p_j].
```

The intermediate alpha is therefore required. A tensor-coupled pair
Hamiltonian makes the 3S1--3D1 blocks of C nonzero, so omitting `alpha_j !=
alpha_c` is not an uncoupled-channel optimization.

### Independent and adversarial checks

The central evidence is executable dense algebra with deliberately
nonsymmetric matrices, not agreement between two code paths sharing an
assumption. `test_3nf_operator_oracle` separately checks pointer layout,
right-hand C orientation, real S/D block restructuring, isolated `C^T W1 C`,
the complete operator, row/column equality, zero W1, cache/fallback equality,
operator order, C versus C^T, the dense equation residual, Neumann convergence,
and the production Padé routine.

OpenCode `paratera/GLM-5.2` supplied additional read-only adversarial reviews:

- Task A, session `ses_01906b545ffeR65TtdVS8PX305`, independently derived the
  alpha contraction before production was edited and identified the old
  oracle's alpha-diagonal blind spot.
- Task B, session `ses_018ed54e9ffeKUQ4TIHfelJfPM`, challenged the transpose
  convention; after reduction to the pointer-layout oracle it accepted the
  failing regression as an independent discriminator.
- Task C, session `ses_018cad73bffeOPvAqW6DIUT0Gp`, tried to falsify the kernel
  patch across indices, ordering, cache/fallback, row/column, dense, and Padé
  paths and found no required change.
- Session `ses_0188411aaffedlO6SsIHTMex2C` separately reviewed the Padé repair,
  including order-zero bounds, NaN handling, convergence labels, and C++/Python
  parity, and returned `ACCEPT` with no required change.

These model reviews are cross-checks, not substitutes for the executable
oracles or a published physical benchmark.

## Fixed in this task

1. A regression with two alpha channels, nonsymmetric nonzero `C_01`/`C_10`,
   and alpha-off-diagonal W1 was committed before the production fix. Against
   the old implementation it produced 64 mismatches and max error 1.378.
2. Both identity-path realizations now contract `sum_alpha_j W1_(r,j)C_(j,c)`.
   Cache lookup uses `(alpha_r,alpha_j)`, not `(alpha_r,alpha_c)`.
3. Duplicated W1 identity/permutation application was consolidated into one
   helper used by the row and column builders. Dense and Padé/Neumann routes
   therefore share the same kernel algebra.
4. The C/C^T storage convention was isolated in a coupled-channel helper and
   locked with a nonsymmetric 3S1--3D1 test.
5. Unverified c1/c3 rank-2 and cD higher-l/rank-2 production paths were removed
   from the advertised calculation and made fail-closed. The remaining model
   keeps its honest `chiral_N2LO_c1c3cDcE_approx` name.
6. W1 quadrature defaults were raised from N=1 midpoint to N=2, with an explicit
   warning that N=1 is diagnostic only and N=2 is not a convergence certificate.
7. Padé false convergence was fixed: two order-minus-one reads were removed,
   breakup uses its own completion counter, no element freezes before the
   configured maximum order, and convergence requires a stable final tail.
   When that test fails, the closest finite adjacent pair is only an explicitly
   labelled optimal-truncation heuristic.

No global scale, fitted normalization, manual sign flip, or post-hoc Hermitian
average was used. `w1_scale` remains solely a loud diagnostic/fault-injection
knob and physical runs require `w1_scale=1`.

## Still approximate

### c1/c3

The production rank-0 formula replaces angular dependence inside nonlinear
pion propagators by a monopole/azimuthal average. The independent full-vector
angular oracle gives:

| `(p,q,p',q')` fm^-1 | production/full | absolute relative error |
|---|---:|---:|
| `(0.5,0.4,0.6,0.7)` | 0.6584 | 34.2% |
| `(0.5,0.5,0.5,0.5)` | 0.5821 | 41.8% |
| `(1.0,0.5,0.8,0.6)` | 1.0386 | 3.9% |
| `(0.3,0.3,0.7,0.7)` | 0.8044 | 19.6% |
| `(1.0,1.0,1.5,1.5)` | 0.8847 | 11.5% |
| `(0.2,0.2,0.2,0.2)` | 1.3407 | 34.1% |

The sign of the error changes with momentum, so no constant factor can repair
it. A trustworthy replacement requires retaining the actual q2 and q3 vectors
through a full angular/PW projection and validating selected matrix elements
against an independent 5D oracle. That research-sized implementation was not
invented here.

### c_D beyond the enabled sector

Only rank-0 spectator S-wave cD is enabled. Existing candidate rank-2
recoupling tests are useful groundwork, but they do not establish phases,
normalization, and all angular projections in production. Higher partial waves
remain zero by design until such an oracle exists.

## Not implemented

c4 remains absent. Its `tau1.(tau2 x tau3)` isospin matrix is independently
checked to be off-diagonal and imaginary in the pair-isospin basis, but the
required spin/momentum cross-product and full PWD have not been implemented.
Nonzero c4 therefore fails closed; cloning c1/c3 or guessing a sign would not be
acceptable.

No converged, publication-quality physical nd benchmark was completed. In
particular, the task did not establish Np/Nq, two-body J, total three-body J,
energy interpolation, deuteron binding, unitarity, and observable convergence
together on one calculation.

## Numerical convergence requirements

- Use at least `Np_per_WP_W1=Nq_per_WP_W1=2` as a starting point. Compare N=2
  with N=4; use N=8 when N2/N4 is not stable. N=1 is a legacy diagnostic.
- Refine `Np_WP` and `Nq_WP` until U/S/observables are stable. The Np=Nq=5
  runs in this report are oracles and debugging calculations, not physics.
- Increase `J_2N_max` and `two_J_3N_max` until the observable of interest is
  stable. No J-truncation certificate was produced here.
- Require the Padé sidecar to report `Conv=1` for a physical claim. `Conv=2`
  means max-order truncated even if an optimal truncation happens to match a
  dense small-grid solve. For persistent `Conv=2`, compare directly with dense
  where feasible or increase/change the resummation strategy.
- Require cache-on/cache-off equality at the chosen W1 quadrature. Changing
  LECs, cutoff, regulator, grids, channel truncations, or quadrature order must
  produce a distinct cache key.
- Check the deuteron bound-state energy, elastic S-matrix/unitarity diagnostics,
  energy placement, and observable stability before comparison with experiment.

## Reproducible commands

All commands start in the repository root.

### Build and automated tests

```bash
cmake --build build -j4
ctest --test-dir build --output-on-failure

python3 -m unittest \
  tests/test_190mev_data_pipeline.py \
  tests/test_2nf_miller_baseline.py \
  tests/test_3nf_matrix_elements.py \
  tests/test_3nf_physics.py \
  tests/test_3nf_regression.py \
  tests/test_coupling_coefficients.py

python3 tools/3nf_oracle/angular_oracle.py
python3 tools/3nf_oracle/wp_quadrature_convergence.py
python3 -m py_compile examples/reconstruct_u_from_neumann.py
```

Observed result: CTest 10/10; Python 25 tests passed and 5 conditionally
skipped; angular and WP reports completed without failure. The core operator
binary reports 26 passes, zero failures; the matrix-element binary reports 431
passes, zero failures.

### End-to-end 2NF and diagnostic W1=0 equivalence

```bash
mkdir -p /tmp/tictac-final-none /tmp/tictac-final-zero

common=(
  two_J_3N_max=1 Np_WP=5 Nq_WP=5 J_2N_max=1 Nphi=48 Nx=48
  Np_per_WP=8 Nq_per_WP=8 Np_per_WP_W1=2 Nq_per_WP_W1=2
  P123_omp_num_threads=2 tensor_force=true isospin_breaking_1S0=true
  midpoint_approx=false calculate_and_store_P123=false
  include_breakup_channels=false solve_faddeev=true solve_dense=false
  production_run=true potential_model=N2LOopt
  energy_input_file=CPP/Input/lab_energies.txt P123_folder=cache/p123
  cache_root=cache
)

OMP_NUM_THREADS=2 ./build/bin/Tic-tac "${common[@]}" \
  three_nucleon_force=none output_folder=/tmp/tictac-final-none

OMP_NUM_THREADS=2 ./build/bin/Tic-tac "${common[@]}" \
  three_nucleon_force=chiral_N2LO_c1c3cDcE_approx \
  c_D=-0.2 c_E=-0.205 Lambda_3NF=500 w1_scale=0 \
  output_folder=/tmp/tictac-final-zero

cmp /tmp/tictac-final-none/U_PW_elements_Np_5_Nq_5_JP_1_-1_Jmax_1_PSI_0.txt \
    /tmp/tictac-final-zero/U_PW_elements_Np_5_Nq_5_JP_1_-1_Jmax_1_PSI_0.txt
cmp /tmp/tictac-final-none/U_PW_elements_Np_5_Nq_5_JP_1_1_Jmax_1_PSI_0.txt \
    /tmp/tictac-final-zero/U_PW_elements_Np_5_Nq_5_JP_1_1_Jmax_1_PSI_0.txt
cmp /tmp/tictac-final-none/neumann_terms_Np_5_Nq_5_JP_1_-1_Jmax_1.txt \
    /tmp/tictac-final-zero/neumann_terms_Np_5_Nq_5_JP_1_-1_Jmax_1.txt
cmp /tmp/tictac-final-none/neumann_terms_Np_5_Nq_5_JP_1_1_Jmax_1.txt \
    /tmp/tictac-final-zero/neumann_terms_Np_5_Nq_5_JP_1_1_Jmax_1.txt
```

All four `cmp` commands returned success in the audited run. The warning from
`w1_scale=0` is intentional; this setting is not a physical 3NF calculation.

### Small dense cross-check

Repeat the first solver command with

```text
solve_dense=true calculate_and_store_P123=false
output_folder=/tmp/tictac-final-dense
```

and compare its two `U_PW_elements_*.txt` files with the Padé files. The audited
Np=Nq=5 2NF run gave max `6.57420273e-9 MeV` and mean
`4.82005695e-10 MeV` over 40 complex elements. The checked-in baseline and its
provenance are in `tests/data/baseline_2nf_only/README.md`.

### Physical-workflow negative check

```bash
python3 examples/quick_Ay_test.py
python3 -m unittest tests/test_190mev_data_pipeline.py
```

The data-pipeline test passes, but the quick Ay calculation did not establish a
190-MeV physical benchmark: its available solver grid selected 161.145 MeV and
the observable comparison failed. This command is retained as a reproducible
open validation target, not a successful benchmark.

## Remaining blockers for a trustworthy full N2LO result

1. Implement and independently verify the full angular/PW projection of c1,
   c3, and the disabled tensor/higher-wave cD pieces.
2. Implement c4 only with its complete isospin, spin, momentum cross-product,
   phases, and an independent oracle.
3. Establish WP, W1-cell, J-truncation, energy, Padé, and unitarity convergence
   simultaneously on a physical nd benchmark.
4. Compare converged observables and selected amplitudes with primary published
   references, rather than tuning any global factor to a curve.

Until these are complete, the scientifically correct statement is: the AGS
operator assembly and its coupled-channel W1 contractions are verified, the
enabled cE and S-wave rank-0 cD pieces have independent matrix-element evidence,
and the advertised interaction remains an explicitly approximate c1/c3/cD/cE
model with no c4 and no completed physical nd validation.
