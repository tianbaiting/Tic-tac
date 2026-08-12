# GLM-5.2 Task A: coupled-alpha contraction review

> **Superseded for the outer permutation order (2026-08-12).** This review
> correctly identified the complete intermediate-alpha contraction, but it
> accepted the then-provided `W1(1+P)` contract without checking the primary
> elastic AGS equation. Deltuva PRC 80, 064002, Eq. (7a), fixes the Tic-tac
> elastic kernel as `(1+P)W1`; see `docs/three_nf_equation_contract.md`.

**OpenCode session:** `ses_01906b545ffeR65TtdVS8PX305`

**Model:** `paratera/GLM-5.2`

**Scope:** independent, pre-production-edit audit of the `C^T W^(1) C` and
`C^T W^(1) P C` contractions.

**Repository HEAD at review:** `6d8e177` on `fix/3nf-physics-contract`.

This file records the independent review conclusion before any production
kernel edit. Existing documentation and comments were explicitly treated as
non-authoritative; the reviewer derived the compound-index contractions and
then traced the implementation.

## Independent result

For compound indices `I=(alpha,p,q)`, the pair-state rotation is diagonal in
spectator momentum. The required identity contribution is

```text
(C^T W1 C)_(alpha_r p_r q_r, alpha_c p_c q_c)
 = sum_(alpha_x,p_x) (C^T)_(alpha_r p_r,alpha_x p_x)
   sum_(alpha_j,p_j) W1_(alpha_x p_x q_r,alpha_j p_j q_c)
                     C_(alpha_j p_j,alpha_c p_c).
```

The repository storage convention was independently traced to

```text
C_(alpha_j p_j,alpha_c p_c)
 = CT_RM_array[alpha_c*Nalpha + alpha_j][p_c*Np + p_j].
```

Therefore `W1*C` requires an intermediate `alpha_j` sum. The corresponding
W1 lookup must use `(alpha_r,alpha_j)`, not `(alpha_r,alpha_c)`.

The pre-fix production code does not perform this sum:

- `calculate_CPVC_col` selects only
  `CT_RM_array[alpha_c*Nalpha+alpha_c]` and calls W1/cache with
  `(alpha_r,alpha_c)`.
- `calculate_all_CPVC_rows` repeats the same restriction.
- The `W1*P*C` contribution is different: `calculate_PVC_col` sums the alpha
  index in `P*C`, and the following W1 application correctly uses its
  intermediate `alpha_k` in both direct and cache paths.
- The final left multiplication by `C^T` also sums its independent alpha
  index correctly.

Thus dense and Padé/Neumann paths are mutually self-consistent only because
they inherit the same missing identity-path terms; that agreement is not an
independent proof of the contracted operator.

## Why the omitted terms are physical

The coupled SWP construction diagonalizes a `2*Np` pair Hamiltonian for tensor
channels and populates the off-diagonal 3S1--3D1 alpha blocks of `C`. The c_D
and c1/c3 rank-2 W1 pieces also admit 3S1--3D1 matrix elements. Consequently

```text
sum_(alpha_j != alpha_c,p_j)
  W1_(alpha_r,alpha_j) C_(alpha_j,alpha_c)
```

is not eliminated by the conserved quantum-number guards and need not vanish.

## Oracle blind spot

The pre-fix `test_3nf_operator_oracle.cpp` constructs only alpha-diagonal C
blocks. Its dense reference therefore has zero omitted terms by construction,
even though its mock W1 is alpha-off-diagonal. GLM ran the existing binary and
confirmed the blind baseline: 7 passes, 0 failures.

The required regression must use at least two conserved-compatible alpha
channels, nonsymmetric nonzero `C_01` and `C_10`, and alpha-off-diagonal W1.
The cleanest discriminator isolates `C^T W1 C` from permutation and 2NF terms;
the full operator test must additionally retain the independent
the corrected `C^T[PV+(1+P)W1]C` comparison and exercise both column and row
builders.  The older formula in the original GLM response is superseded.

## Review confidence and open item

The missing `alpha_j` sum, the C/C^T orientation above, the production tensor
mixing, and the oracle blind spot are supported by direct code evidence. The
size of the resulting change in physical observables remains unquantified and
must be measured only after a regression-first fix. The cache table's complete
key metadata and its off-diagonal build loop still require a separate local
audit; Task A established only that the call sites use the wrong right-channel
alpha in the identity path.

## Process note

Although instructed to remain read-only, the OpenCode reviewer initialized an
untracked `.codegraph/` directory while locating symbols. It did not modify any
tracked production or test file. The generated directory is excluded from all
task commits.
