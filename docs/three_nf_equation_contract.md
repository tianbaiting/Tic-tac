# Three-Nucleon-Force Equation Contract for Tic-tac

**Status:** authoritative reference. Every C++ object in the 3NF / Faddeev path
must be traceable to a symbol defined here. **No empirical scaling factor,
post-hoc Hermitian averaging, or hand-flipped sign is permitted to make a curve
agree with a reference.** When a measurement disagrees with this contract, the
*derivation* in this document is the court of appeal, not the reference curve.

This document is Phase 0 of the fix branch `fix/3nf-physics-contract`. It is the
prerequisite for every later phase: the finite-dimensional operator oracle
(Phase 0 acceptance), the c₄ / c_E fixes (Phases 1–2), the independent angular
oracle (Phase 3), and the kernel unification (Phase 5) all read their *expected*
operator algebra from §6 below.

---

## 1. Conventions and notation

| Symbol | Meaning |
|--------|---------|
| `p`, `q` | Jacobi momenta in the spectator-1 frame: `p` = pair-relative (particles 2,3), `q` = spectator (particle 1). Units: fm⁻¹. |
| `G0 = (E − H0 + iε)⁻¹` | free 3-body resolvent |
| `t` | fully-dressed two-body t-matrix in the (2,3) pair subspace (a function of `p,p'` at fixed spectator `q`) |
| `P = P₁₂₃ + P₁₃₂` | sum of cyclic particle permutations. For fully **antisymmetric** pair states `P₁₂₃ = P₁₃₂`, so `P = 2·P₁₂₃`. |
| `V^(1)` | the spectator-1 component of the 3NF: the operator is symmetric under 2↔3 and the full force is `V₄=V₄^(1)+V₄^(2)+V₄^(3)`. |
| `W^(1)` | the bare spectator-1 3NF component used by the code; `W^(1)=V₄^(1)=u` in the references below. |
| `φ`, `Φ` | asymptotic deuteron-spectator channel state, antisymmetric in the pair only |
| `C` | SWP→WP basis-rotation matrix (the `C_WP` arrays) |
| `G` | **channel resolvent**: the SWP-diagonal resolvent obtained after rotating through `C`. `G ≠ G0`. |

Literature lock (re-audited 2026-08-12):

- E2002 — Epelbaum et al., PRC **66** (2002) 064001. Operator forms, LECs,
  regulator eq. (3.19), contact PWD App. A.
- G2010 — Golak et al., EPJ A **43** (2010) 241. General 5D PWD algorithm.
- H2015 — Hebeler et al., PRC **91** (2015) 044001. LS-coupled factorised PWD.
- W2011 — Witała et al., PRC **83** (2011) 044001, Eqs. (49)--(50).
  Faddeev breakup component and elastic operator with a 3NF.
- D2009 — Deltuva, PRC **80** (2009) 064002, Eqs. (1)--(10), especially
  the symmetrised elastic AGS Eq. (7a).
- M2022 — Miller et al., PRC **106** (2022) 024001, Eqs. (1) and (10).
  The 2NF channel-resolvent/WPCD equation implemented by Tic-tac.

The former citation to Witała et al., PRC **77**, 034004 (2008), Eq. (3), was
wrong: that equation is a basis-coupling relation in a relativity/Ay paper and
does not define a 3NF AGS kernel.

---

## 2. The continuous Faddeev equation with 3NF

The starting point is W2011 Eq. (49) (equivalently E2002 Eq. (3.15)): the
spectator-1 Faddeev breakup component `T` of the full
3-body scattering operator, driven by both the 2-body t-matrix and the 3NF:

```
T|φ⟩ = t·P|φ⟩
      + (1 + t·G0)·V^(1)·(1 + P)|φ⟩              … (i) 3NF Born
      + t·P·G0·T|φ⟩                              … (ii) 2NF kernel
      + (1 + t·G0)·V^(1)·(1 + P)·G0·T|φ⟩        … (iii) 3NF kernel
```

The four pieces have distinct roles:

1. `t·P|φ⟩` — 2NF Born (one application of the dressed pair t, with spectator
   permutation `P`).
2. `(1 + t·G0)·V^(1)·(1 + P)|φ⟩` — **3NF Born**, direct. The `(1+P)` antisymmetrises
   the *ket*; the `(1+t·G0)` dresses the 3NF by the pair t-matrix.
3. `t·P·G0·T|φ⟩` — 2NF rescattering kernel.
4. `(1 + t·G0)·V^(1)·(1 + P)·G0·T|φ⟩` — **3NF rescattering kernel**. Same dressing
   `(1+t·G0)` as (ii); the extra `G0·T` is one more rescattering.

The unknown operator is `T`, the spectator-1 Faddeev breakup component.  The
breakup operator is `U0=(1+P)T`.  The elastic operator is different and is
given by W2011 Eq. (50); one must not transfer the right-hand `(1+P)` in the
breakup-component equation directly into the channel-resolvent AGS kernel.

---

## 3. Reduction of the symmetrised AGS equation to WPCD form

D2009 derives the multichannel AGS equations from the full resolvent and then
specialises to identical nucleons. With `u=W^(1)=V4^(1)`, its Eq. (7a) is

```
U = P·G0^-1 + (1+P)·u
  + P·t·G0·U
  + (1+P)·u·G0·(1+t·G0)·U .
```

This `U` is the symmetrised elastic transition operator. It is the operator
used in M2022 Eq. (1) and represented by Tic-tac's `U_array`; it is not the
breakup component `T` in §2.

Let `G=G0+G0·t·G0=(E-H0-V+i0)^-1`. The Lippmann--Schwinger equation gives
the two exact identities

```
t·G0 = V·G,
G0·(1+t·G0) = G.
```

Therefore the two iteration terms in D2009 Eq. (7a) reduce without reordering
any noncommuting factors:

```
P·t·G0·U                         = P·V·G·U,
(1+P)·u·G0·(1+t·G0)·U           = (1+P)·u·G·U.
```

For an incoming channel eigenstate `|phi>` of `H0+V`,
`G0^-1|phi>=V|phi>`. Hence the driving term has the same on-shell action as
`P·V+(1+P)·u`. The channel-resolvent/WPCD equation is consequently

```
U = K + K·G·U,
K = P·V + (1+P)·W^(1).                         (contract)
```

This reduces to M2022 Eq. (1) when `W^(1)=0`.

### 3.1 Operator ordering

The order is fixed by D2009 Eq. (7a): `(1+P)` is on the **left** of the bare
spectator component. `W^(1)·(1+P)` occurs in the Faddeev breakup-component
equation of §2, but that equation solves a different operator. Moving its
right symmetriser into the elastic channel-resolvent equation is invalid.

The distinction is observable before external-state symmetrisation whenever
`[P,W^(1)] != 0`. `tests/cpp/test_faddeev_operator_order.cpp` evaluates the
D2009 kernel with Hermitian noncommuting finite matrices, checks `P^2=P+2`,
verifies both resolvent identities above, accepts
`[P V+(1+P)W^(1)]G`, and rejects `[P V+W^(1)(1+P)]G`.

### 3.2 Dressing and scope

`W^(1)` is the bare spectator component. No `(1+tG0)` is inserted into
`W1_element`: the pair dressing is already carried by the channel resolvent in
`K G U`. Conversely, `G` contains only the pair interaction; the 3NF enters
the driving matrix `K` and every subsequent Neumann/Padé iteration through
that matrix.

The breakup amplitude still requires the separate W2011 relation
`U0=(1+P)T`. The present Tic-tac breakup-output path has not yet been audited
against D2009 Eq. (7b)/(9b); no publication claim is made for 3NF breakup
amplitudes until that audit is complete.

---

## 4. WPCD / SWP discretisation of the AGS equation

### 4.1 Bases

Tic-tac uses three nested bases for the (α, p, q) discretisation:

1. **WP basis** (wave packets): momentum-space bins `[p_i, p_i+Δp_i]` and
   `[q_j, q_j+Δq_j]` with midpoint quadrature. State label
   `|α, i, j⟩_WP`. Dimension `Nα·Np_WP·Nq_WP` (= `dense_dim` in the code).
2. **SWP basis** (shifted wave packets): the WP basis diagonalising the
   pair Hamiltonian `H_pair = H0_pair + V_pair`. The transformation is `C`
   (`C_WP_unco_array`, `C_WP_coup_array`). SWP eigenvalues are the
   `e_SWP_*_array` used by `make_resolvent`.
3. **Channel basis**: the (α) partial-wave label `α = (L₂N, S₂N, J₂N, T₂N,
   L₁N, 2j₁N; 2J₃N, 2T₃N, P₃N)`. Conservation of `2J₃N, 2T₃N, P₃N` block-diagonalises
   every operator.

### 4.2 The WP-basis AGS kernel

Rotate the AGS kernel `K = P·V + (1+P)·W^(1)` into the WP basis and apply the
SWP rotation `C` on the pair line:

```
A_α'α  ≡  (C^T · K · C)_{α'α}
        =  (C^T · P · V · C)_{α'α}                  … 2NF driving  (PV term)
         +  (C^T · W^(1) · C)_{α'α}                 … 3NF identity  (W1·C term)
         +  (C^T · P · W^(1) · C)_{α'α}              … 3NF permutation (P·W1·C term)
```

This is exactly the three contributions the code assembles inside
`calculate_CPVC_col` / `calculate_all_CPVC_rows`:

| Math object | C++ array / function |
|-------------|----------------------|
| `C` (pair SWP rotation) | `C_WP_unco_array`, `C_WP_coup_array`; viewed row-major as `CT_RM_array[α·Nα+α']` |
| `V` (2NF pair potential) | `V_WP_unco_array`, `V_WP_coup_array`; combined with `C` into `VC_CM_array[α·Nα+α']` |
| `P = P₁₂₃ + P₁₃₂` (cyclic permutation, sparse) | `P123_sparse_{val,row,col}_array` |
| `W^(1)` (bare spectator-1 3NF matrix element) | `three_nucleon_force_model::W1_element(α_r, α_c, p_r, q_r, p_c, q_c, pw_states)` |
| `PVC = P·V·C` column | `calculate_PVC_col` → `PVC_col[]` |
| `(1+P)·W^(1)·C` column | `add_one_plus_P_W1_C_col`: first contract `W^(1)·C`, then add the identity column and its sparse left-permuted `P·W^(1)·C` column |
| `A = C^T·(PVC + W1·C + P·W1·C)` column | `calculate_CPVC_col` → `col_array[]`; row form `calculate_all_CPVC_rows` → `row_arrays[]` |

The two code paths (`calculate_CPVC_col` and `calculate_all_CPVC_rows`) compute
**the same operator `A`** — one column-at-a-time (for the sparse precalculated
kernel and the dense solver) and one row-at-a-time (for the on-shell driving
term and the Neumann iteration). Both call the single
`add_one_plus_P_W1_C_col` helper and must agree to machine precision
(operator-level oracle, §6).

### 4.3 The Neumann / Padé iteration

Once `A = C^T·K·C` is formed, the AGS equation `U = K + K·G·U` becomes the
matrix iteration

```
U_WP = A + A·G_WP·U_WP ,
```

with `G_WP` the SWP-diagonal channel resolvent (`G_array` in the code, built by
`make_resolvent` from `e_SWP_*_array`). Tic-tac expands this as the Neumann
series

```
a_0 = A ,   a_{n+1} = a_n · G_WP · A ,    U_WP ≈ Σ_{n=0}^{N} a_n ,
```

then Padé-resums `[N/N]` of the `a_n` to recover the on-shell elastic
`U_array`. The 3NF enters **only** through `A` (i.e. through `K`); the
resolvent `G_WP` is *independent* of the 3NF LECs `c_D, c_E, c_1, c_3, c_4`
(the channel resolvent is the *2-body* pair resolvent, not a 3-body object).

### 4.4 How `(1+P)·W^(1)·C` is assembled

The expensive contraction is performed once to form a dense `W^(1)·C` column.
The identity contribution is added directly. The same column is then used as
the right vector of the stored CSC permutation matrix, giving
`P·W^(1)·C`. This order is both cheaper than recomputing a second 3NF
contraction and faithful to D2009 Eq. (7a).

---

## 5. Operator ordering: the four questions answered

| Question (from task §Phase 0) | Answer |
|------------------------------|--------|
| What is the unknown operator? | The symmetrised elastic AGS `U` of D2009 Eq. (7a) and M2022 Eq. (1), represented by `U_array`. It is not the Faddeev breakup component `T`. |
| Why is the 2NF driving `P·V` (not `P·t`)? | On the incoming channel, `P G0^-1|phi>=P V|phi>`, while `tG0=VG` converts the iteration term to `P V G U`. |
| What is the 3NF driving, `V^(1)·(1+P)` or `(1+P)·V^(1)`? | `(1+P)·V^(1)`, fixed without reordering by D2009 Eq. (7a). The right-ordered form belongs to the distinct breakup-component equation. |
| Is `(1+t·G0)` fully absorbed? | Yes: `G0(1+tG0)=G`. The code evaluates the bare `V^(1)` and the channel resolvent supplies the pair dressing. |
| Does the direct `(1+P)·V^(1)` term appear once in the elastic amplitude? | Yes, as the 3NF part of the `n=0` driving matrix `A=C^T K C`; later occurrences are rescattering iterations. |
| Does breakup `(1+P)·T` include the 3NF? | The primary equations say yes, but Tic-tac's 3NF breakup-output mapping is not yet independently audited; it remains outside the present publication claim. |

---

## 6. Finite-dimensional operator-level oracle (Phase 0 acceptance)

Before any physics tuning, the algebra of §4.2 must hold on a small
artificial matrix. The oracle lives in
`tests/cpp/test_3nf_operator_oracle.cpp` and is deliberately independent of
the production 3NF recoupling / kernels: it builds *random* but
symmetry-respecting `P, V, W1, G1, C` and checks that **every code path** that
claims to compute `A = C^T·(P·V + (1+P)·W^(1))·C` agrees with the mathematical
expression, and that the resulting Neumann / Padé iteration reproduces
`U = (1 − A·G)⁻¹·A`.

### 6.1 Setup

Finite dimensions: `Nα = 2`, `Np = Nq = 2` (so `DIM = 8`). The mock operators
are **deliberately non-symmetric** so that transposition / ordering bugs are
not masked:

- `P` — a generic sparse permutation in the production-contraction oracle; the
  separate ordering test uses the physical sum with `P^2=P+2`.
- `V` — random pair potential, `V ≠ V^T` (asymmetric test).
- `W1` — random spectator-1 3NF, `W1 ≠ W1^T` (asymmetric test).
- `C` — random SWP rotation, `C ≠ C^T` (so `C` and `C^T` are distinguished).
- `G1` — diagonal channel resolvent, `G1[i] = 1/(E − e_i + iε)`.

A single shared seed makes the test reproducible.

### 6.2 Equalities checked (all to `1e-10`)

Let `M_math = C^T·(P·V + (1+P)·W1)·C` (dense, assembled directly from the
definitions). Then:

1. **`calculate_CPVC_col` path** — for each column `j`,
   `col_array[j] == M_math[:, j]`. This is the sparse column builder used by the
   dense solver and the precalculated-kernel path.
2. **`calculate_all_CPVC_rows` path** — for each on-shell row `i`,
   `row_array[i, :] == M_math[i, :]`. This is the row builder used by the
   Neumann driving term.
3. **Dense solver path** — `faddeev_dense_solver` with `K = M_math` and
   `G = G1` produces `U_dense = (1 − M_math·G1)⁻¹·M_math`, compared to the
   analytic inverse.
4. **Padé / Neumann path** — the `a_n = A·(G·A)^n` series, Padé-resummed,
   matches `U_dense` to within the Padé tolerance (looser, `1e-6`, because
   Padé is an approximation).
5. **Operator-ordering guard** — the wrong-order kernel
   `M_wrong = C^T·(P·V + W1·(1+P))·C` is *not* equal to `M_math`, and the
   difference `||M_math − M_wrong||_max > 0` (sanity, not a tolerance).
6. **Transpose guard** — with asymmetric `W1`, the code path
   `calculate_CPVC_col` must compute `W1·C` (not `C·W1` and not `W1^T·C`); a
   deliberate swap inside the test confirms the expected value matches
   `W1·C` and not the alternatives.

This oracle is the **gate** for Phase 3 (independent angular oracle): until
every check in §6.2 passes, no claim about partial-wave matrix elements is
admissible.

---

## 7. Unit contract

The contract fixes the units of every object. No `(2π)³` factor may be inserted
empirically; if a `(2π)³` appears, it must be derivable from the Fourier
convention below.

| Object | Unit | Source |
|--------|------|--------|
| `p, q, p', q'` (Jacobi) | fm⁻¹ | definition |
| `m_π, f_π, Λ_χ, Λ_3NF` | fm⁻¹ | constants.h converted via `÷ħc` |
| `c_1, c_3, c_4` | fm (= GeV⁻¹ × `ħc/1000`) | constructor `m_c1 = c1·ħc/1000` |
| `c_D, c_E` | dimensionless | constructor `m_c_D = c_D` |
| `V^(1)(p',q'; p,q)` (bare 3NF) | fm⁵ | from `(gA/2fπ)²·fm² × 1/fπ² × ...`; closed-form for c_E |
| `p_WP, q_WP` (stored WP boundaries) | MeV | consistent with 2NF `V_WP` |
| `G_array` (channel resolvent) | MeV⁻¹ | from `(E − e_SWP)⁻¹` |
| `W1_WP` (WP-basis 3NF, after bin averaging) | MeV | `W1_raw[fm⁵] × (p_r·q_r·p_c·q_c)·√(dp_r·dq_r·dp_c·dq_c) × 1/ħc⁵`, where the four momenta (MeV) + four `√bin` (MeV^{1/2} each) give MeV⁶, times `1/ħc⁵` (MeV⁻⁵·fm⁻⁵ → but `W1_raw` is fm⁵, so fm⁵×MeV⁶×(1/ħc⁵) = MeV⁶×MeV⁻⁵ = MeV). |

**Fourier convention.** Tic-tac applies `1/(2π)³` for each independent relative
coordinate.  The 2NF has one such coordinate; a 3NF Jacobi state has `p` and
`q`, hence a raw 3NF partial-wave projection is multiplied by `(2π)^-6`.
For the contact, Epelbaum Eq. (A-4) supplies the raw S-wave angular factor
`(4π)²`, so the normalized spectator coefficient is
`(4π)²(2π)^-6 = 1/(4π⁴)`.  The independent five-angle projector reproduces
both this contact factor and Golak Table 2 before the Tic-tac normalization is
applied.  The WP cache integrates only the radial measure and contains no
hidden `(2π)^3` cancellation.  The old single factor `1/(8π³)` made `cE` too
large by `π/2`; W1 cache schema v5 prevents reuse of those blocks.

The incomplete `cD` and `c1/c3` rank-zero development kernels retain a clearly
named legacy normalization until they are replaced by the exact full angular
projection.  Their numerical normalization is not part of the locked complete
N2LO contract.

**Regulator.** Squared-Gaussian per E2002 eq. (3.19):
`f_R(p,q;Λ) = exp(−((4p²+3q²)/(4Λ²))²)`, applied symmetrically to bra and ket.
No alternative regulator form is admissible without amending this contract.

---

## 8. What this contract forbids (carried into the audit report)

The following are **forbidden** as physics-fix tools on this branch (they may
appear only as labelled debug fault-injection, gated to `w1_scale=1.0` in
production):

1. `w1_scale ≠ 1.0` in any production run. The parameter is retained as a
   debug knob; non-unity values must print a strong warning and are
   contract-violating.
2. Global empirical multiplicative factors applied to `W^(1)` to match a
   reference (the historical `~2.4–3.0×` residual is explicitly disallowed).
3. Hand-flipped signs in the kernel to fix a magnitude/sign discrepancy.
   E2002 Eq. (2.10) has `+1/2 E`; any alternative LEC convention must be
   translated explicitly and documented at the input boundary.
4. Post-hoc Hermitian symmetrisation of `W^(1)` (e.g. averaging
   `W1(r,c) + W1(c,r)^T`); the manifestly Hermitian kernel must be derived from
   the 5D angular integral (Phase 3 oracle), not imposed after the fact.
5. Azimuthal-angle averaging dressed up as "true partial-wave projection" in
   comments or docs; the current c_D / c_1 / c_3 path is a restricted rank-zero
   approximation and is labelled as such (it is *not* the full G2010 5D PWD).
6. Advertising the model as `chiral_N2LO` (full) while `c_4` is unimplemented;
   the honest name is `chiral_N2LO_c1c3cDcE_approx` until c_4 is independently
   verified (Phase 1).

---

## 9. Acceptance of Phase 0

Phase 0 is accepted when:

- [x] This document exists and every C++ array in §4.2 maps to a math object.
- [x] `tests/cpp/test_3nf_operator_oracle.cpp` exists and passes its production
      contraction, row/column, cache, dense, Neumann, Padé, and ordering checks
      of §6.2 on the finite-dimensional random system.
- [x] The algebra oracle is independent of production 3NF recoupling/kernels:
      it builds its own `P, V, W1, G1, C` and compares against direct dense
      assembly. This acceptance is limited to contraction algebra and does not
      validate physical `W1_element` values.
- [x] `tests/cpp/test_faddeev_operator_order.cpp` independently reduces D2009
      Eq. (7a) and rejects the former reversed order for Hermitian,
      noncommuting matrices.

These checks accept the elastic AGS/WPCD algebra only. The full-vector physical
oracle and complete N2LO partial-wave projection remain later publication
gates.
