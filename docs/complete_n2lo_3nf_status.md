# Complete chiral N2LO 3NF: implementation and validation status

**Audit date:** 2026-08-13
**Branch:** `fix/3nf-physics-contract`
**Audited baseline:** `46b7308`, plus the validation artifacts listed below

## Scope and publication gate

The target is a convention-locked, independently validated implementation of the
complete local chiral N2LO three-nucleon force (3NF), including the
`c1`, `c3`, `c4`, `cD`, and `cE` operators, its partial-wave projection, wave-packet
(WP) matrix elements, and its use in the three-body scattering equation.  A
low-energy neutron-deuteron analyzing-power calculation is an end-to-end
benchmark, not a substitute for operator-level validation.

**Current disposition: not publication-ready.**  The legacy fast model is
deliberately named `chiral_N2LO_c1c3cDcE_approx`.  A complete, solver-selectable
direct-Jj model now exists under the explicit name
`chiral_N2LO_full_5d_reference`; it is an exact-operator reference with
`O(Nangle_3NF^5)` cost.  The complete three-integral implementation is now
solver-selectable as `chiral_N2LO_full_factorized`, but remains provisional
until the WP, cache, solver, and observable convergence gates are complete.
The unqualified name `chiral_N2LO` therefore remains fail-closed.

## Evidence ledger

The labels used below are:

- **Verified result:** reproduced in this repository by the listed test or direct
  source inspection.
- **Literature fact:** read from the cited primary paper.
- **Inference/open item:** a conclusion that still needs an independent derivation,
  implementation, or numerical test.

### Repository and regression state

| Item | Evidence | Status |
|---|---|---|
| C++ build | `cmake --build build -j4` | **Verified result:** successful; existing compiler/HDF5 warnings remain. |
| C++ tests | `ctest --test-dir build --output-on-failure` | **Verified result:** 13/13 passed. |
| Direct 3NF oracle test | `build/tests/test_3nf_operator_oracle` | **Verified result:** 26 passed, 0 failed. |
| Faddeev ordering discriminator | `build/tests/test_faddeev_operator_order` | **Superseded audit-start result:** the former test only proved conformity to the repository's incorrect assumed ordering.  It has been replaced by a primary-equation discriminator described below. |
| 3NF matrix-element tests | `build/tests/test_3nf_matrix_elements` | **Verified result:** 432 passed, 0 failed; the complete `c4` requirement remains an expected failure/skip. |
| Python regressions | `python3 -m unittest tests/test_190mev_data_pipeline.py tests/test_2nf_miller_baseline.py tests/test_3nf_matrix_elements.py tests/test_3nf_physics.py tests/test_3nf_regression.py tests/test_coupling_coefficients.py` | **Verified result:** 25 passed, 5 skipped.  The 40 checked 2NF amplitudes reproduce the stored baseline with `max|delta| = 0`. |
| Full-vector/PWD oracle | `python3 -m unittest tests/test_full_vector_n2lo_oracle.py` | **Verified result:** 13/13 passed.  The generic five-angle projector agrees with the independently transcribed Golak integrands at `N=4` to 10 decimal places, reproduces the published Table 2 values at `N=12` within `3e-4` relative, and gives identical direct-Jj and unitary-9j-transformed LS projections. |
| Complete C++ reference PWD | `build/tests/cpp/test_chiral_n2lo_full_reference` | **Verified result:** all signed component and combined checks passed.  Direct Jj `c1`, `c3`, `c4`, `cD`, and `cE` values agree with the independent Python projector at the same `Nangle_3NF=2` to roughly `1e-17` absolute.  Both low-order forward and reverse values are frozen separately; no post-hoc Hermitian averaging is used. |
| Factorized Python PWD | `python3 -m unittest tests/test_factorized_scalar_pwd.py tests/test_factorized_n2lo_pwd.py` | **Verified result:** 5/5 passed.  The Hebeler Eq. (6) three-integral scalar kernel gives the exact `(4pi)^2` contact limit, agrees with the independent five-angle projector for P-wave and `l=0<->2` finite-rank kernels, reproduces Golak Table 2 for `c1/c3/c4` within `3e-4` relative, and matches complete `cD` S-wave and spectator-D transitions. |
| Factorized C++ PWD | `build/tests/cpp/test_chiral_n2lo_factorized` | **Verified result:** signed `c1/c3/c4/cD/cE` values agree with the independently tested Python factorization to about `1e-17` absolute for S, P, pair-D, `J=3/2`, and `T=3/2` channels.  Exact forward/reverse Hermiticity holds, and the independent five-angle result approaches it from `N=4` to `N=6`. |
| Broad N2LO matrix table | `tools/3nf_oracle/generate_n2lo_matrix_table.py`; `output/validation/n2lo_3nf_matrix_elements.json`; `python3 -m unittest tests/test_n2lo_matrix_validation.py` | **Verified result:** 26/26 signed matrix elements, 130/130 separate component checks, and 13/13 forward/reverse pairs pass.  The table contains 58 nonzero and 72 selection-rule-zero comparisons, both signs for every `c1/c3/c4/cD/cE` component, three momentum regimes, S/P/pair-D/spectator-D channels, both parities, and `J=1/2,3/2`, `T=1/2,3/2`.  Production C++ versus the separate Python factorization differs by at most `1.17e-15 fm^5`; the independent full-vector five-angle result differs by at most `2.83e-6 fm^5`, within its explicit order-6-to-8 error estimate.  Production order 8 versus 10 differs by at most `3.10e-8 fm^5`. |
| Factorized W1 performance | `output/validation/n2lo_3nf_factorized_performance.json` | **Verified implementation result:** a fixed single-thread 200-channel-pair pass fell from `233.55 s` at `208159a` to `65.50 s` at `b39901e` (3.57x).  Subsequent magnetic factorization reduced the production `Np=4,Nq=3`, 101-block radial-order-one W1 build from `102.3 s` to `31.0 s`; channel batching reduced the radial-order-two build from `162.4 s` to `87.8 s`.  Exact-Hermitian reuse at `d14d2bb` directly integrates 56 blocks, transpose-fills 45, and reduces the same three-thread build to `48.4 s` (another 1.81x).  Commit `2563e8c` prevents the P123-specific three-thread limit from leaking into W1.  Commit `9eaa5f8` fuses the shared `c1/c3` reduced integration and removes repeated momentum powers, ordered-cache lookup, and weight-table mutex traffic.  Commit `964d4d7` additionally hoists the batch-common regulator and memoizes immutable LS/isospin/quadrature metadata: fresh 1/16/96-thread builds now take `53.9/11.7/10.4 s`, versus `126.5/16.9/12.7 s` before the two hot-path passes.  All 101 W1 payloads, U files, and convergence files are bitwise identical. |
| W1 metadata-only warm-up | `three_nucleon_force_model::prepare_W1_channel`; `output/validation/n2lo_3nf_factorized_performance.json` | **Verified neutral optimization:** the factorized model now populates LS expansions and Cartesian/spin weight tables without evaluating and discarding 56 complete transfer integrals.  The 1/16/96-thread builds take `127.8/16.9/13.8 s`, versus `130.4/17.3/13.0 s` for the former full-element warm-up.  The outputs are byte-identical, but the high-thread timing is within run variation, so no material speedup is claimed.  Removing warm-up entirely was rejected because it increased the 96-thread time to `22.9 s`. |
| Complete W1 cell integration | `python3 -m unittest tests/test_factorized_wp_quadrature.py`; `output/validation/n2lo_3nf_wp_quadrature.json` | **Verified result:** for ordinary and deliberately wide cells, the legacy midpoint has relative errors from `2.28e-3` to `1.90e-1`.  The hardest wide `c4+cD` transition still has `3.84e-3` error at radial order two, while order four agrees with order six to `7.89e-7`. |
| W1 cache off/miss/hit parity | `build/tests/cpp/test_chiral_n2lo_w1_cache` | **Verified result:** at radial order two, an uninitialized-cache direct build, an HDF5-cache miss/store build, and a subsequent all-block disk hit are bitwise identical.  The two-bin-per-axis test exhausts the swapped `p/q` layout and proves that every reverse slab is an exact transpose.  A counting wrapper proves triangular evaluation on a miss and zero calls on a hit; the wide `c4+cD` transition agrees with the independent Python cell integral to `5.2e-11` relative.  Schema v7 excludes older independently accumulated reverse blocks. |
| Reduced solver interface smoke | `build/bin/Tic-tac ... three_nucleon_force=chiral_N2LO_full_factorized Np_WP=1 Nq_WP=1 ...` | **Verified interface result only:** factory construction, basis generation, permutation, W1 cache construction, AGS assembly, and Padé code completed without an interface error.  The selected energy file yielded no on-shell energies and the radial grids used one midpoint, so this is not solver-physics or convergence evidence. |
| Historical reduced complete-3NF dense/Padé cross-check | `tools/3nf_oracle/compare_reduced_solver_outputs.py`; `output/validation/n2lo_3nf_reduced_solver_crosscheck.json` | **Superseded Padé-stability baseline:** before commit `5e003d6`, one `J=1/2`, positive-parity block at `Np=4`, `Nq=3`, radial order one, and angular order two agreed with dense inversion to `max|delta U|=2.12e-9 MeV`, but the determinant-ratio evaluator marked all 12 elements maximum-order truncated.  The grid and unphysical `-0.0114 MeV` deuteron energy remain useful negative WPCD evidence; the Padé status has been superseded by the next row. |
| Stable Padé resummation | `output/validation/n2lo_3nf_pade_stability.json`; commit `5e003d6` | **Verified numerical result:** replacing the high-order ratio of determinants with a row-scaled denominator-coefficient solve leaves the Neumann coefficients byte-identical and does not change the fixed three-step, `1e-5` relative plus `1e-7` absolute convergence test.  At radial order two, all 12 amplitudes are genuinely converged at both `[14/14]` and `[24/24]`; `[14/14]` agrees with dense inversion to `max|delta U|=5.79e-15 MeV`, and orders 14 versus 24 differ by at most `1.34e-14 MeV`.  Cached `(4,3)`, `(6,4)`, `(8,6)`, and `(10,8)` reruns give `76/76` code-1 amplitudes at the unchanged default order 14. |
| Reduced WPCD convergence diagnostic | `tools/3nf_oracle/summarize_factorized_solver_convergence.py`; `output/validation/n2lo_3nf_reduced_wpcd_convergence.json` | **Verified negative grid result:** with `J_2N_max=1`, transfer order two, and midpoint W1 integration, increasing `(Np,Nq)` through `(4,3)`, `(6,4)`, `(8,6)`, and `(10,8)` moves the deuteron energy from `-0.0114` to `-1.2586 MeV`, so the 2NF/WP basis is still far from converged.  Different q grids have different on-shell bin midpoints and their U values are not compared directly.  At identical `(4,3)`, radial order one versus two changes the 12 amplitudes by as much as `0.1335 MeV` and `67.4%` relative.  The archived JSON contains pre-`5e003d6` code-2 sidecars; cached post-fix reruns resolve those Padé codes without resolving the WPCD drift. |
| Idaho-N3LO deuteron p-grid ladder | `deuteron_binding_only=true`; `output/validation/idaho_n3lo_deuteron_binding_ladder.json`; `python3 -m unittest tests/test_deuteron_binding_mode.py` | **Verified production-path result:** the fail-closed mode calls the ordinary `V_WP` builder and SWP diagonalizer but returns before P123.  Independent `p/q` Chebyshev controls preserve legacy inheritance exactly and isolate binding from `q` tuning.  A two-dimensional shape scan followed by the fixed `p_t=0.88,p_s=325` ladder gives binding magnitudes `2.14681,2.18846,2.20378,2.21107,2.21510,2.21844 MeV` at `Np=40,60,80,100,120,150`.  A common `1/Np^2` fit has `E_infinity=2.2242145 MeV`, `0.360 keV` below the reference, with at most `0.048 keV` residual.  `Np=82` is the lowest tested point inside the present 20-keV gate (`19.828 keV` error). |
| Realistic-grid Ay feasibility | `output/validation/n2lo_3nf_physical_Ay_feasibility.json` | **Verified mixed result:** the candidate `Np=82,Nq=12,p=(0.88,325),q=(1,90)` places a bin midpoint at `Tlab=10.0453 MeV`.  An actual 2NF `J=1/2`, both-parity run builds schema-v2 exact-grid P123 caches in about 64 s and solves in about 152 s; all 16 amplitudes converge at `[24/24]` and `[32/32]`, with maximum order-to-order difference `1.28e-9 MeV`.  Exact-Hermitian reuse, restored 96-thread execution, and two bitwise-equivalent projection optimizations lower the complete-3NF W1 projection from about 164 to 19.4 h for only one 101-block low-order J-parity sector.  All 100 ordered channel pairs in each parity are nonzero at the scan points, so exact-zero pruning is unavailable.  This is a feasibility rejection, not an Ay result. |
| W1 low-rank/streaming feasibility | `output/validation/n2lo_3nf_w1_low_rank_feasibility.json` | **Verified rejection for production use:** W1 is consumed element-wise while forming `W1*C`, so streaming changes memory rather than the number of costly W1 evaluations, and a post-build SVD pays the dense construction cost first.  At block tolerance `1e-4`, the schema-v9 `(4,3)` probe changes U by `7.49e-5` in relative L2 and up to `1.23e-3` pointwise.  At `1e-6`, U changes by `8.86e-7`, but factor storage is larger than dense at `(4,3)` and only about 5% smaller in the archived `(10,8)` structural probe.  No low-rank production path is enabled. |
| Distributed/resumable exact W1 construction | `docs/w1_distributed_construction.md`; `tests/cpp/test_w1_distributed.cpp`; `tools/refactor_harness/run_w1_block_dist_acceptance.sh`; commits on `main` | **Verified implementation result (exact, no approximation):** the expensive W1 construction factors into independently-computable exact blocks (one per `(alpha_r, alpha_c)` channel pair, classified `evaluate`/`transpose_fill`), shared by a single `integrate_w1_channel_blocks` helper that both the monolithic `W1_PW_cache::build` and a per-block `W1BlockExecutor` call (bitwise-identical by construction).  Each block publishes via an atomic temp+rename shard with full `W1Key` (schema v9) provenance; a stale/corrupt shard is a clean miss.  Block-level `--worker-index/--worker-count` partitioning lets N workers (or N machines) attack one sector in parallel; resume skips completed blocks with zero recomputation.  Measured on the `input_j3nf_multiblock.txt` grid: the block-level 16-worker, 2-worker, and interrupted(7-block)+resumed builds all produce the **identical** W1 fingerprint `611ce32acac0a3395178cb6b5aed5eb3bcaaed2d9550d5d555b08eeebb0456df`, matching the monolithic sector-level build; the downstream solver U is bitwise-identical across both parities; 90/90 Hermitian reverse-block transpose pairs verify `bad=0`; `ctest` 17/17 pass including the new `w1_distributed` contract test.  Single-machine block-level scaling (1/4/16 workers: 23/22/19 s) is memory-bound (consistent with the existing 16→96-thread 1.12× ratio); the operational value is cross-machine distribution (cost model: ~19.4 h/sector ÷ N machines) and resumability.  This removes the operational blocker; the realistic Ay campaign is still not complete. |
| Low-energy nd Ay preflight | `examples/audit_low_energy_Ay.py`; `output/validation/n2lo_3nf_low_energy_Ay_preflight.json`; `python3 -m unittest tests/test_low_energy_ay_audit.py` | **Verified negative acceptance result:** the fail-closed audit computes `Ay_max`, its angle, peak deficits, RMSE, uncertainty-weighted chi-square, maximum residual, and `R_Ay`, but permits an interpretation only after all physics gates pass.  The existing paired `Np=Nq=10`, `J_3N<=9/2` Idaho-N3LO artifacts match in numerical settings and have all requested J-parity U blocks, but use the old approximate 3NF, have `E_d=-0.108641 MeV`, contain 400 truncated 2NF and 51 truncated 3NF amplitudes, and provide none of the required `Np/Nq/W1/J2N/J3N/Pade` Ay ladders.  Their diagnostic `R_Ay=1.1595` is therefore explicitly withheld as a physical conclusion. |
| Full Python discovery | `python3 -m unittest ...` including `tests/test_pade_honesty.py` | **Environment limitation:** collection fails because the local Python environment does not provide `pytest`; this is separate from the clean `unittest` regression set above. |

At audit start the worktree also contained three unrelated untracked user files:
`CPP/Input/input_j9_3nf.txt`, `examples/plot_latest_3nf_comparison.py`, and
`scripts/launch_j9_3nf_parallel.sh`.  They are outside this audit's edit scope and
are being preserved.

### Primary-source formula audit

1. **Literature fact:** Epelbaum et al., *Phys. Rev. C* **66**, 064001
   ([arXiv:nucl-th/0208023](https://arxiv.org/abs/nucl-th/0208023),
   [DOI](https://doi.org/10.1103/PhysRevC.66.064001)) give the N2LO two-pion,
   one-pion-contact, and contact operators in Eqs. (2.2) and (2.10), with
   `D = cD/(fpi^2 Lambda_chi)` and `E = cE/(fpi^4 Lambda_chi)` in Eq. (2.12).
   Their Eq. (2.10) has a **positive** `+1/2 E` coefficient for each ordered
   `j != k` term.  Each unordered pair occurs twice, so the spectator-1
   component is `u1=E tau2·tau3`, not `(E/2) tau2·tau3`.  Production and the
   independent Pauli enumeration now use this unit coefficient.

2. **Literature fact:** Appendix A of the same paper contains the full partial-wave
   expressions.  In particular, the `cD` term is not restricted to the present
   spectator-S-wave rank-zero piece.  Visual inspection of the rendered source
   page confirms that Eq. (A-4) reads `6 E (4pi)^2`, not `E/2`.  The `6j`
   values reduce this to `(4pi)^2 E tau23`.

3. **Literature fact:** Golak et al., *Eur. Phys. J. A* **43**, 241
   ([arXiv:0911.4173](https://arxiv.org/abs/0911.4173),
   [DOI](https://doi.org/10.1140/epja/i2009-10903-6)) reduce the direct
   three-nucleon-force projection from eight angular integrations to five, not
   three.  Their Table 2 supplies numerical `c1/c3/c4` matrix elements suitable
   for an independent fixed-momentum oracle.

4. **Literature fact:** Hebeler et al., *Phys. Rev. C* **91**, 044001
   ([arXiv:1502.02977](https://arxiv.org/abs/1502.02977),
   [DOI](https://doi.org/10.1103/PhysRevC.91.044001)) provide a factorized exact
   projection of local three-nucleon forces.  Their Jacobi spectator convention
   differs from Tic-tac's and must be translated explicitly before reuse.

5. **Verified result:** Tic-tac's 2NF convention applies `(2pi)^-3` to its one
   relative coordinate.  A 3NF has the independent `p` and `q` coordinates,
   so the raw PWD requires `(2pi)^-6`.  The WP cache later supplies only radial
   measures and `hbarc` conversion.  Combining this with Epelbaum A-4 gives the
   exact contact normalization `1/(4pi^4)`; the former `1/(8pi^3)` coefficient
   was too large by `pi/2`.

6. **Verified citation defect:** Witala et al., *Phys. Rev. C* **77**, 034004
   ([arXiv:0801.0367](https://arxiv.org/abs/0801.0367),
   [DOI](https://doi.org/10.1103/PhysRevC.77.034004)) do not establish the
   repository's asserted 3NF AGS kernel.  Their Eq. (3) is a basis-coupling
   formula, and their scattering equation is the 2NF-only Eq. (50).  The present
   citation in `docs/three_nf_equation_contract.md` is therefore invalid.

7. **Literature fact and corrected implementation:** Deltuva, *Phys. Rev. C*
   **80**, 064002 ([arXiv:0912.0240](https://arxiv.org/abs/0912.0240),
   [DOI](https://doi.org/10.1103/PhysRevC.80.064002)), Eq. (7a), gives the
   symmetrised elastic AGS equation.  Using `tG0=vG` and
   `G0(1+tG0)=G` reduces its iteration kernel exactly to
   `[P v + (1+P)W1]G`.  This is the same elastic transition operator and
   channel resolvent used in Miller et al., *Phys. Rev. C* **106**, 024001,
   Eqs. (1) and (10).  The former code order `W1(1+P)` mixed the elastic AGS
   `P v` ordering with the distinct Faddeev breakup-component equation and was
   wrong.  Production now forms `W1*C` once and adds both `W1*C` and
   `P*W1*C`.

## What the code currently implements

| Component | Current implementation | Assessment |
|---|---|---|
| `cE` contact | Spin scalar, pair-isospin eigenvalue, exact `alpha_r=alpha_c`, `l_pair=lambda=0` selection, regulator, unit spectator-component coefficient, and closed `tau23*cE/(4*pi^4*fpi^4*Lambda_chi)` normalization | Verified by explicit Pauli states, the generic five-angle projector, visual A-4 inspection, signed C++ golden values, and all 26 cE rows of the broad table to `4.18e-16 fm^5`.  The closed form removes transfer-order-two leakage into higher orbital channels; W1 cache schema v8 rejects the old blocks. |
| `cD` one-pion-contact | Pair-contact/spectator-S-wave rank-zero subset | Incomplete.  Required rank-two and higher-orbital structures are absent. |
| `c1`, `c3` two-pion exchange | Diagonal rank-zero azimuthal/monopole approximation | Incomplete.  Off-diagonal and higher-rank angular-momentum couplings are absent. |
| `c4` in legacy fast model | Constructor rejects nonzero `c4` | Missing by design; the approximate model cannot represent the complete N2LO force. |
| Complete five-angle reference | Explicit Jj angular-spin states, 8-state Pauli spin and isospin algebra, all five operator components, regulator, rotational volume, and `(2pi)^-6` | Implemented as `chiral_N2LO_full_5d_reference` and signed-oracle tested.  Correctness-first only: cost scales as `Nangle_3NF^5`, so it is not accepted for converged WPCD production grids. |
| Factorized three-integral implementation | Hebeler Eq. (6) scalar kernel plus an automated Cartesian-vector/spherical-harmonic finite-rank expansion and explicit spin/isospin matrices | Complete in Python and C++ for `c1`, `c3`, `c4`, `cD`, and `cE`.  It retains three nontrivial integrals, caches momentum-independent angular/spin weights, applies the nonlocal regulator after PWD, and performs a unitary LS-to-Jj recoupling.  It is solver-selectable as `chiral_N2LO_full_factorized`; realistic WP cost and convergence remain open gates. |
| Regulator | Squared nonlocal Gaussian associated with Epelbaum Eq. (3.19) | Present; convention and cutoff pairing must remain explicit in every benchmark. |
| WP cache | Four-dimensional radial-bin quadrature with model/coupling/grid/truncation fields and schema-v9 operator versioning in the key | `Nangle_3NF`, the distinct reference model name, and numerical `gA`, `fpi`, `mpi`, `Lambda_chi`, and `hbarc` values prevent cross-projector/order/convention reuse.  Independent production-driver integration shows radial order four is converged below `8e-7` for the tested difficult cells.  Cache-off, miss/store, and hit paths are bitwise identical at order two. |
| Scattering insertion | Code builds `W1*C` and then its sparse left-permuted `P*W1*C`, including the complete intermediate-channel contraction | Corrected to the primary-source kernel `(1+P)W1`.  The noncommuting test now derives the reference matrix directly from Deltuva Eq. (7a). |

The implementation also exposes `w1_scale`.  It is correctly marked as a
diagnostic fault-injection control and must remain exactly one in physical runs.

## Independent-oracle status

The former `tools/3nf_oracle/angular_oracle.py` remains a diagnostic for the
legacy approximate production kernels.  It is not the authoritative full
operator oracle.

The independent reference stack is now:

- `full_vector_n2lo_oracle.py`: explicit 64-state spin-isospin matrices and
  full Cartesian `c1,c3,c4,cD,cE` spectator-1 operators;
- `full_vector_five_angle_pwd.py`: generic LS-coupled five-angle projection
  with explicit magnetic-substate sums and no production recoupling reuse;
- `golak_table2_benchmark.py`: a separate transcription of Golak Eq. (25)
  that reproduces `G(1,1)=443.618 fm^5` and `G(2,1)=1200.219 fm^5`.
- `factorized_scalar_pwd.py` and `factorized_n2lo_pwd.py`: an Eq. (6)
  three-integral implementation in which Cartesian Jacobi-vector factors are
  converted into finite `l-1,l+1` spherical-harmonic expansions.  The
  momentum-independent magnetic/spin algebra is summed separately from the
  transfer-magnitude kernels.

Verified discriminators include Jacobi momentum conservation, all five LEC
switches, reverse-kernel Hermiticity, `2<->3` symmetry, the contact eigenvalues,
`c4` zero/nonzero cross products, and
`<t23'=0|tau1.(tau2 x tau3)|t23=1>=-2 sqrt(3) i`.  The generic projector also
returns the raw contact factor `(4pi)^2 E tau23` and the normalized
`E tau23/(4pi^4)` result.  Its direct Jj angular states agree pointwise with
the unitary 9j transformation of the LS basis, including a multicomponent
P-wave channel; direct and recoupled matrix elements agree for `c1/c3` and an
off-diagonal `c4/cD` discriminator.

The legacy oracle limitations are retained here to prevent accidental reuse:

- its `cE` coefficient and Fourier factor duplicate production conventions;
- its `cD` calculation covers only the rank-zero S-wave subset;
- its advertised "full" `c1/c3` integral fixes both incoming Jacobi directions
  and integrates only three angular variables, whereas the primary projection
  requires five after rotational reduction;
- it contains no `c4` operator and no complete magnetic-quantum-number summation.

Consequently, only the new full-vector/five-angle stack may serve as the
production PWD reference.  The C++ class
`chiral_N2LO_3NF_full_reference` is a separate transcription of that algorithm:
it uses GSL spherical harmonics and Clebsch--Gordan coefficients plus explicit
Pauli actions, and matches the Python signed values without importing its
implementation.  The scalable `chiral_N2LO_3NF_factorized` class is a second
C++ algorithm, validated against both this reference and the independent
Python Hebeler factorization; its matrix-element implementation is complete,
while its cache, solver, and observable convergence gates remain open.

At low angular order, reverse matrix elements need not yet agree because the
two fixed-axis representations have different quadrature errors.  For the
off-diagonal `c4+cD` discriminator at `Nangle_3NF=2`, both C++ and Python give
`0.0109108167629193 fm^5` forward and `0.0135450686980236 fm^5` reverse.  Their
difference is retained as an explicit non-convergence signal.  Hermiticity is
therefore a successive-order gate, not an operation applied to the result.

## Numerical benchmark status

Existing 10 MeV analyzing-power results are diagnostic only:

- `output/ay_3nf_current/Ay_2nf_vs_3nf_10MeV.json` reports an Ay RMSE change
  from about `0.1903` (2NF) to `0.2150` (approximate 3NF), with zero fully
  converged 2NF Padé elements and an unphysical deuteron binding energy near
  `-0.109 MeV`.
- `output/ay_diagnosis_current_binary/Ay_j9_2nf_vs_3nf_10MeV.json` reports an
  Ay RMSE change from about `0.0829` to `0.0961`; it is also explicitly marked
  diagnostic and has incomplete Padé convergence.

These calculations use coarse WP/angular truncations and the incomplete force.
They do not satisfy the requested low-energy nd Ay publication benchmark.

The dedicated P123-free deuteron ladder isolates the leading problem in those
artifacts.  It exactly reuses the full solver's potential-packet and SWP
diagonalization stages, reproduces the old `Np=10` header, and shows monotone
`1/Np^2` convergence toward the physical Idaho-N3LO binding energy.  The
default `Np=10` value is therefore a finite-packet artifact, not a different
physical prediction.  The solver now exposes independent `p` and `q` Chebyshev
controls; zero preserves the former common-control behavior exactly.  The
tested candidate `Np=82,p_t=0.88,p_s=325` is inside the 20-keV binding window,
while `Nq=12,q_t=1,q_s=90` places a scattering midpoint at `Tlab=10.0453 MeV`.
This establishes a feasible starting grid, not `Nq` or Ay convergence.

The candidate was exercised through an actual 2NF `J=1/2`, both-parity solve.
P123 construction, memory, and Padé are viable: all 16 amplitudes are honestly
converged at both `[24/24]` and `[32/32]`, and the two orders differ by at most
`1.28e-9 MeV`.  P123 schema v2 keys the exact boundary arrays rather than only
the common Chebyshev controls, preventing false hits for independent or custom
grids.

Complete-3NF W1 construction remains the launch blocker.  The candidate has
6,724 times as many radial cells per channel block as the measured `(4,3)`
case.  Exact factorized-model Hermiticity reduces the measured 101-block
radial-order-two build from 87.8 to 48.4 seconds by evaluating 56 orientations
and transpose-filling 45 at the previously leaked three-thread setting.  After
restoring independent W1 thread control, 96 threads reduce it to 13.0 seconds;
the first projection hot-path optimization reduces the same build to 11.1 seconds.
Hoisting the batch-common regulator and memoizing immutable LS, isospin, and
quadrature metadata reduces it again to 10.4 seconds.  Linear cell-count
scaling still gives about 19.4 hours for a single low-order
J-parity sector, before higher W1 orders, the J ladder, both parities, and the
paired 2NF/3NF runs.  Its roughly 746-MiB dense
W1 payload is manageable; compute is not.  A full ordered-channel scan found
no exact-zero blocks to prune.  A quantified SVD audit also rejects immediate
production low-rank or streaming work: streaming leaves construction unchanged,
while useful compression at tolerance `1e-4` visibly perturbs U and the safer
`1e-6` tolerance provides little factor-storage benefit.  Further exact
factorization or distributed resumable construction is required before a
responsible full launch.

The first reduced-grid solve with the complete factorized force adds an
algorithmic check, not a physics benchmark.  On this machine, constructing 101
W1 blocks for `Np=4`, `Nq=3`, radial order one, angular order two, and four
OpenMP threads took roughly 6.5 minutes; a later all-hit W1 load took about
0.2 seconds.  The cache-hit Padé solve completed in about 7.4 seconds and the
dense solve in about 0.31 seconds.  The large, strongly channel-dependent first
build cost is an explicit scalability issue for higher-order realistic grids.

The former code-2 Padé tails were traced to numerical evaluation rather than to
the discrete Faddeev equation.  The old implementation formed two complex
determinants and divided them; by order 24, the individual products overflowed,
underflowed, or became `0/0` even though the rational approximant was finite.
Commit `5e003d6` instead solves directly for the denominator coefficients,
uses a rank-revealing fallback for defective tables, and evaluates numerator
and denominator by Horner's rule.  No convergence tolerance was relaxed.  The
same stored Neumann histories now produce stable final tails and dense-level
agreement, while the separate WPCD and radial-quadrature failures remain.

## Open convention ledger

No physical kernel or normalization should be changed until each item below has
an explicit derivation, primary-source anchor, and noncommuting or numerical
discriminator.

1. **Resolved 2026-08-12:** the elastic channel-resolvent kernel is
   `P V + (1+P)W1`, derived from Deltuva Eq. (7a).  The distinct
   `W1(1+P)` order remains valid inside the Faddeev breakup-component equation,
   not the elastic `U` equation solved by Tic-tac.
   A GLM-5.2 adversarial audit argued for the opposite order by identifying
   `U_array` with Deltuva's unsymmetrised `X` (Eq. 10).  That finding is rejected:
   production solves `(I-A G)U=A` with the Miller 2NF kernel `A=P v`, whereas
   the `X` equation has a different 2NF order.  The same audit also treated
   Witala's reconstruction of elastic `U` from breakup-component `T` as the
   independent equation being solved.  The primary Eq. (7a) reduction and the
   noncommuting discriminator therefore remain authoritative.
2. **Resolved 2026-08-13:** the Jacobi, transfer, Pauli, Fourier, regulator,
   radial-state, spherical-harmonic, and phase conventions are locked in
   `docs/n2lo_3nf_conventions.md`.
3. **Resolved 2026-08-13:** the `cE` ordered-pair, A-4 angular, Fourier, sign,
   and regulator factors are exact-golden tested.  Old caches are invalidated.
4. **Resolved 2026-08-13:** the full-vector and generic five-angle oracles
   reproduce independent closed integrands and Golak Table 2.
5. **Resolved at the matrix-element implementation and broad-table level:** a second C++ direct-Jj
   five-angle implementation contains every `c1`, `c3`, `c4`, `cD`, and `cE`
   structure and is available through the factory with explicit slow-reference
   naming.  `Nangle_3NF` is printed, parsed, and hashed.  The Hebeler
   three-integral factorization is translated and independently validated in
   Python for all five components, then ported to a solver-selectable C++ model.
   The machine-readable table now covers 26 signed forward/reverse matrix
   elements and 130 component comparisons across diagonal/off-diagonal, S/P/D,
   tensor, both-parity, multiple-momentum, `J=1/2,3/2`, and `T=1/2,3/2`
   sectors.  Production, a separately implemented finite-rank PWD, and the
   independent full-vector five-angle projector are all retained in every row,
   including their last quadrature steps.  Realistic-grid WP/solver convergence
   remains required; no unsupported dimensional reduction is accepted.
6. **Resolved for representative-cell cache parity and quadrature:** W1 schema
   v9 hashes all chiral constants and the angular order, locks the exact-
   Hermitian transpose construction contract, and uses the exact cE contact
   selection rule plus the analytic cD relative-angle integral rather than
   finite transfer quadrature for those angle-independent pieces.  Direct integration
   through the production C++ driver demonstrates
   that the midpoint is not reliable and that radial order four agrees with
   order six below `8e-7` in the hardest tested wide transition.  Cache-off,
   miss/store, and hit paths are bitwise equal at order two, with the disk hit
   making no operator calls.  Repeat the ladder on realistic production grids
   before fixing a default.
7. **Padé numerical stability resolved at diagnostic grids; physics convergence
   still open:** the unchanged default `[14/14]` sequence now gives `76/76`
   genuinely converged amplitudes across the cached grid ladder.  On the radial-
   order-two `(4,3)` case, it agrees with dense inversion below `5.8e-15 MeV`
   and with `[24/24]` below `1.4e-14 MeV`.  The coarse-grid deuteron binding is
   still unphysical and strongly grid dependent.  Demonstrate stable two-body
   binding, WP/rank/J convergence, symmetry/Hermiticity/permutation checks,
   zero-LEC and 2NF-only limits, and a reproducible low-energy nd Ay comparison
   with uncertainty/convergence tables.  The realistic candidate passes the
   2NF `[24/24]` versus `[32/32]` test, but the dense W1 builder is not yet
   computationally viable for the complete paired calculation.

## Acceptance decision

The repository now satisfies the equation-ordering, full-vector-oracle,
published raw-PWD, exact-contact-normalization, complete slow-C++-reference,
factorized matrix-element implementation, and broad machine-readable matrix-table
gates.  Representative W1 cell integration is also converged at radial order
four.  Independent grids now supply a binding-qualified realistic starting
point, and its 2NF P123/Padé path is measured to be viable.  The
**operational** blocker to the complete-3NF realistic-grid W1 build is removed:
W1 construction now factors into independently-computable exact blocks with
atomic-shard persistence, full-provenance cache keys, block-level distribution
across N workers/machines, and verified resumability (completed blocks are never
recomputed; all construction paths produce an identical SHA-256 W1 fingerprint;
downstream U is bitwise-identical) — see
`docs/w1_distributed_construction.md`.  No low-rank, SVD, compression, or other
physics approximation was introduced.  The repository still fails the decisive
complete-3NF realistic-grid **execution** gate: the realistic `Np=82, Nq=12`
W1 build has not been run to completion (its cost is a linear extrapolation
~19.4 h per low-order J-parity sector, divisible across machines but not yet
executed), and the higher-W1-order / J / parity / paired-2NF+3NF convergence
ladders remain open.  The low-energy Ay path has a machine-readable, fail-closed
acceptance audit and a separate feasibility report; both reject interpretation of
the existing diagnostic curves.  Converged radial/angular/J ladders and the
paired result must precede the final comparison.  No existing Ay curve is yet
admissible as a complete-N2LO result.

## GLM-5.2 review disposition

An OpenCode `paratera/GLM-5.2` read-only architecture review was used as an
adversarial second opinion.  Its channel-field map and recommendation to keep a
slow reference distinct from the scalable model were accepted.  The following
claims were rejected:

- it equated the legacy `1/(8pi^3)` factor with the two-coordinate
  `(2pi)^-6` Fourier normalization; these are not equal and the accepted
  normalization is fixed by the contact and Golak gates above;
- it asserted without derivation that the complete `c1/c3` terms require only
  one angular integral and `c4` only two.  Hebeler Eq. (6) leaves three
  nontrivial integrals for a general local spin-dependent interaction after
  analytic factorization;
- it described `tau1.(tau2 x tau3)` as real.  In the chosen product basis its
  off-diagonal recoupling is purely imaginary, for example
  `-2 sqrt(3) i`; the complete `c4` matrix element becomes real only after the
  spin and isospin factors combine.

No rejected GLM statement has been incorporated into production code.
