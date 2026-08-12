# Complete chiral N2LO 3NF: implementation and validation status

**Audit date:** 2026-08-12  
**Branch:** `fix/3nf-physics-contract`  
**Audited commit:** `9d402b04bf986dda3d0ec043296b0d3dc9a9544e`

## Scope and publication gate

The target is a convention-locked, independently validated implementation of the
complete local chiral N2LO three-nucleon force (3NF), including the
`c1`, `c3`, `c4`, `cD`, and `cE` operators, its partial-wave projection, wave-packet
(WP) matrix elements, and its use in the three-body scattering equation.  A
low-energy neutron-deuteron analyzing-power calculation is an end-to-end
benchmark, not a substitute for operator-level validation.

**Current disposition: not publication-ready.**  The current production-selectable
model is deliberately named `chiral_N2LO_c1c3cDcE_approx`.  It is a useful
fail-closed development implementation, but it is not the complete N2LO 3NF and
must not be cited as such.

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
| C++ tests | `ctest --test-dir build --output-on-failure` | **Verified result:** 10/10 passed. |
| Direct 3NF oracle test | `build/tests/test_3nf_operator_oracle` | **Verified result:** 26 passed, 0 failed. |
| Faddeev ordering discriminator | `build/tests/test_faddeev_operator_order` | **Superseded audit-start result:** the former test only proved conformity to the repository's incorrect assumed ordering.  It has been replaced by a primary-equation discriminator described below. |
| 3NF matrix-element tests | `build/tests/test_3nf_matrix_elements` | **Verified result:** 432 passed, 0 failed; the complete `c4` requirement remains an expected failure/skip. |
| Python regressions | `python3 -m unittest tests/test_190mev_data_pipeline.py tests/test_2nf_miller_baseline.py tests/test_3nf_matrix_elements.py tests/test_3nf_physics.py tests/test_3nf_regression.py tests/test_coupling_coefficients.py` | **Verified result:** 25 passed, 5 skipped.  The 40 checked 2NF amplitudes reproduce the stored baseline with `max|delta| = 0`. |
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
   independent Pauli enumeration now use this unit coefficient.  The complete
   Fourier/PW normalization remains open.

2. **Literature fact:** Appendix A of the same paper contains the full partial-wave
   expressions.  In particular, the `cD` term is not restricted to the present
   spectator-S-wave rank-zero piece, and the contact normalization/recoupling in
   Eq. (A-4) still has to be traced into the Tic-tac normalization convention.

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

5. **Verified citation defect:** Witala et al., *Phys. Rev. C* **77**, 034004
   ([arXiv:0801.0367](https://arxiv.org/abs/0801.0367),
   [DOI](https://doi.org/10.1103/PhysRevC.77.034004)) do not establish the
   repository's asserted 3NF AGS kernel.  Their Eq. (3) is a basis-coupling
   formula, and their scattering equation is the 2NF-only Eq. (50).  The present
   citation in `docs/three_nf_equation_contract.md` is therefore invalid.

6. **Literature fact and corrected implementation:** Deltuva, *Phys. Rev. C*
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
| `cE` contact | Spin scalar, pair-isospin eigenvalue, diagonal channel selection, regulator, unit spectator-component coefficient | Ordered-pair counting is independently verified and regression-locked; the complete Fourier/PW normalization and phases still need a primary-source benchmark. |
| `cD` one-pion-contact | Pair-contact/spectator-S-wave rank-zero subset | Incomplete.  Required rank-two and higher-orbital structures are absent. |
| `c1`, `c3` two-pion exchange | Diagonal rank-zero azimuthal/monopole approximation | Incomplete.  Off-diagonal and higher-rank angular-momentum couplings are absent. |
| `c4` two-pion exchange | Constructor rejects nonzero `c4` | Missing by design; production cannot represent the complete N2LO force. |
| Regulator | Squared nonlocal Gaussian associated with Epelbaum Eq. (3.19) | Present; convention and cutoff pairing must remain explicit in every benchmark. |
| WP cache | Four-dimensional radial-bin quadrature with model/coupling/grid/truncation fields in the key | Present, but the key does not yet encode all physical constants or a stable operator-definition version.  Cache/direct parity is only tested for a one-point cell. |
| Scattering insertion | Code builds `W1*C` and then its sparse left-permuted `P*W1*C`, including the complete intermediate-channel contraction | Corrected to the primary-source kernel `(1+P)W1`.  The noncommuting test now derives the reference matrix directly from Deltuva Eq. (7a). |

The implementation also exposes `w1_scale`.  It is correctly marked as a
diagnostic fault-injection control and must remain exactly one in physical runs.

## Independent-oracle status

`tools/3nf_oracle/angular_oracle.py` is presently a diagnostic oracle, not the
required independent unprojected oracle:

- its `cE` coefficient and Fourier factor duplicate production conventions;
- its `cD` calculation covers only the rank-zero S-wave subset;
- its advertised "full" `c1/c3` integral fixes both incoming Jacobi directions
  and integrates only three angular variables, whereas the primary projection
  requires five after rotational reduction;
- it contains no `c4` operator and no complete magnetic-quantum-number summation.

Consequently, existing oracle passes do not validate the full operator.  The
replacement oracle must start from full Jacobi vectors and explicit spin/isospin
states, implement all five LEC terms independently of production PWD code, and
reproduce both direct angular projections and published fixed-momentum matrix
elements.

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
2. Lock the Jacobi spectator, momentum-transfer, permutation, state-normalization,
   Fourier-transform, and spherical-harmonic phase conventions.
3. Trace the remaining `cE` Fourier/state-normalization factors and complete
   contact PWD to Epelbaum Eq. (A-4); the ordered-pair coefficient itself is
   resolved and exact-golden tested.
4. Implement a genuinely independent five-angle/full-vector magnetic-substate
   oracle and reproduce Golak Table 2 before trusting production PWD values.
5. Select and implement an exact full PWD strategy for `c1`, `c3`, `c4`, and
   `cD`; translate any Hebeler factorization into Tic-tac's spectator convention.
6. Version the W1 cache by the complete operator definition and all physical
   constants, then demonstrate cache-on/cache-off equality at quadrature order
   at least two and convergence under higher quadrature.
7. Demonstrate stable two-body binding, Padé honesty, WP/rank/J convergence,
   symmetry/Hermiticity/permutation checks, zero-LEC and 2NF-only limits, and a
   reproducible low-energy nd Ay comparison with uncertainty/convergence tables.

## Acceptance decision

The repository currently satisfies useful infrastructure and regression gates,
but fails the decisive completeness and independent-validation gates.  The
elastic scattering-equation ordering is now convention-locked and corrected.
The next safe implementation milestone is the independent full-vector oracle;
only after it agrees with published fixed-momentum values should the production
PWD and cache format be replaced.
