# Chiral N2LO 3NF Audit and Physics/Numerics Fix — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Audit and fix the chiral N²LO three-nucleon force (3NF) implementation so that the operator structure, partial-wave matrix elements, and Faddeev kernel ordering satisfy literature conventions (Epelbaum 2002, Golak 2010, Hebeler 2015, Witała 2008), with independent golden tests and honest capability reporting.

**Architecture:** Phased commits — audit doc first, then split helpers, lock conventions, fix cache, refresh benchmarks, document capabilities. No empirical scaling factors. No `w1_scale` tuning. Each fix derived from formula + units + basis convention.

**Tech Stack:** C++17 (CMake build under `build/`), GSL coupling coefficients, OpenMP, Python 3 reference oracles (SymPy + direct integration), CTest + pytest.

---

## Confirmed Bugs (from audit of current `master`)

| # | Severity | Location | Bug |
|---|----------|----------|-----|
| B1 | **critical** | `src/utils/chiral_3nf_recoupling.cpp:18-63` (`recoupling_3nf_scalar`), used by `src/interactions/chiral_N2LO_3NF.h:144` (`W1_contact`, the c_E term) | `recoupling_3nf_scalar` returns `(σ₂·σ₃)(τ₂·τ₃)` but the c_E contact term per Epelbaum 2002 eq. (2.10)/(A-4) is **pure spin-scalar**: only `τ₂·τ₃`. So c_E erroneously picks up `σ₂·σ₃` = −3 for S_2N=0 (singlet) or +1 for S_2N=1 (triplet). |
| B2 | **critical** | `src/interactions/chiral_N2LO_3NF.h:49,110` (`m_c4`); `three_nucleon_force_model.cpp:29-39` (passes c4_idaho_n3lo = +5.4) | `m_c4` stored but **never used** in `W1_element`. c_4 term silently dropped even though the model advertises itself as `chiral_N2LO`. |
| B3 | high | `src/utils/chiral_3nf_recoupling.h` (`recoupling_3nf_scalar` declared as "rank-0 pair operator") | Same helper used for **two physically different operators**: c_E (pure τ·τ) and 2PE rank-0 (σ·σ × τ·τ). Cannot serve both correctly. |
| B4 | high | `src/interactions/three_nucleon_force_model.h:42-43` vs `src/core/faddeev_solver/solve_faddeev.cpp:252-258` | Header says `W = (1+P+P²) W^(1)` (LEFT form); solver implements `W^(1)·(1+P)` (RIGHT form). The two are physically different operators unless `[W^(1), P]=0 is proven. Solver is correct per Witała 2008 AGS form; header is misleading. |
| B5 | high | `src/interactions/w1_pw_cache.cpp:138-165` (W1Key) | Cache key omits `c1, c3, c4` and the momentum grid. Two runs with different c_1 (e.g. N2LOopt vs Idaho N3LO) but identical c_D, c_E, Λ would **wrongly reuse** the cached W^(1). |
| B6 | medium | `src/interactions/w1_pw_cache.h:57`, `:54` (`std::vector<float> m_data`) | Cache stores `float` (24-bit mantissa). 3NF has cancellations between c_E/c_D/c_1/c_3 contributions; float precision causes ~1e-7 relative error that compounds in Padé. |
| B7 | medium | `src/interactions/chiral_N2LO_3NF.h:73-77, 168-184, 254-282` | Documentation calls the azimuthal-averaged kernel a "true partial-wave projection". It IS exact for c_E (diagonal) but is an **approximation** for c_D, c_1, c_3 (cross terms Δp·Δq are dropped under φ_p' azimuthal average; the 5D Golak 2010 reduction is not done). |
| B8 | low | `src/interactions/chiral_N2LO_3NF.h:339-393` | 2PE rank-2 Hermitian symmetrization is implemented by **averaging two evaluations** (`(recoup2·integ2 + recoup2_T·integ2_T)/2`) instead of deriving a manifestly Hermitian kernel. Works numerically but obscures the real algebraic structure. |

**Not bugs** (verified): regulator formula (E2002 eq. 3.19 squared-Gaussian — already correct), Fourier normalization 1/(8π³) (matches `chiral_LO_internal.cpp:59`), permutation factor `2·P123_val` for antisymmetric pair-states (correct in `solve_faddeev.cpp:203`).

---

## Conventions (locked)

1. **Spectator-1 picture**: particle 1 = spectator, pair = (2,3). Jacobi: `p` = pair relative, `q` = spectator.
2. **Momentum transfers**: `q₁ = Δq = q' − q`, `q₂ = Δp + Δq/2`, `q₃ = Δp − Δq/2`, `Δp = p' − p`.
3. **Faddeev/AGS kernel** (Witała 2008 eq. 3; Golak 2010): the 3NF piece is `W^(1)·(1+P)` with `P = P₁₂₃ + P₁₃₂`, NOT `(1+P)·W^(1)`. The full physical 3NF on the 3N Hilbert space is `W = (1+P+P²)·W^(1) = W^(1) + W^(2) + W^(3)` but **only** as a property of the antisymmetrized bra/ket — this does NOT authorize commuting `P` past `W^(1)` in the kernel.
4. **c_E contact** (Epelbaum 2002 eq. 2.10 + A-4): `V^(1)_cont = -½·E·Σ_{j≠k∈pair}(τ_j·τ_k) = -E·(τ₂·τ₃)` with `E = c_E/(f_π⁴ Λ_χ)`. **Spin-scalar.** Selection rules: L_2N=L_2N'=0, l_1N=l_1N'=0, j_1N=½, J_2N=S_2N, T_2N'=T_2N, S_2N'=S_2N. Closed-form recoupling (formula_reference.md §1.4):
   ```
   A_cE = 6·(4π)² · δ(...) · (-1)^{T_2N+1} · SixJ(½,½,T_2N; ½,½,T_3N)
   ```
5. **2PE c_1/c_3 rank-0** (Golak 2010 eq. 18): `(σ₂·q₂)(σ₃·q₃) → ⅓(σ₂·σ₃)(q₂·q₃)`. Here `σ₂·σ₃` IS the pair eigenvalue `2S(S+1)−3`, i.e. the **same** `recoupling_3nf_scalar` formula is correct FOR THIS TERM.
6. **c_4** (Epelbaum 2002 eq. 2.2-2.3): carries `τ₁·(τ₂×τ₃)`, purely imaginary for T_3N=½ doublet, opens off-diagonal T_2N transitions. Not implemented — will be **hard-blocked**.
7. **Units**: p, q, m_π, f_π, Λ_χ, Λ all in fm⁻¹ (via `÷ ħc`); c_D, c_E dimensionless; c_1, c_3, c_4 in fm (from GeV⁻¹ × ħc/1000).

---

## File Structure

| File | Responsibility | Status |
|------|----------------|--------|
| `src/utils/chiral_3nf_recoupling.h/.cpp` | Spin-isospin recoupling helpers | **modify**: split c_E from 2PE scalar |
| `src/interactions/chiral_N2LO_3NF.h` | 3NF W^(1) matrix element | **modify**: use c_E helper; hard-block c_4; honest name |
| `src/interactions/w1_pw_cache.h/.cpp` | WP cache | **modify**: double precision; extend key |
| `src/interactions/three_nucleon_force_model.h` | base class doc | **modify**: lock (1+P+P²) vs W^(1)(1+P) language |
| `src/core/faddeev_solver/solve_faddeev.cpp` | solver | **modify**: comments only (kernel already correct) |
| `tests/cpp/test_chiral_3nf_recoupling.cpp` | unit tests | **modify**: add c_E-specific tests |
| `tests/cpp/test_3nf_matrix_elements.cpp` | matrix element tests | **modify**: add golden tests |
| `tests/cpp/test_w1_cache_double.cpp` | **new** | cache double vs float comparison |
| `tests/cpp/test_faddeev_operator_order.cpp` | **new** | dense operator-ordering test |
| `tests/cpp/test_w1_quadrature_convergence.cpp` | **new** | 1×1 to 4×4 convergence test |
| `tools/check_3nf_normalization/hand_calc_cE.py` | **new/refresh** | independent SymPy oracle for c_E |
| `tools/check_3nf_normalization/hand_calc_cD.py` | refresh | independent oracle for c_D |
| `docs/treatise/chapters/15_3nf_physics.tex` | physics | **modify**: capability status |
| `docs/treatise/chapters/16_3nf_pw_projection.tex` | PWD | **modify**: approximation disclosure |
| `docs/treatise/chapters/17_3nf_numerical.tex` | numerical | **modify**: cache/convergence status |
| `docs/3nf_audit_and_validation_2026.md` | **new** | final audit + validation report |

---

## Task 1: Audit document with confirmed bugs and conventions

**Files:**
- Create: `docs/3nf_audit_2026-06-21.md`

- [ ] **Step 1: Write audit doc**
  Document all 8 confirmed bugs with file:line references, the locked conventions (§ Conventions above), literature citations, and the planned fixes. Cross-reference `tools/check_3nf_normalization/formula_reference.md` for each formula.

- [ ] **Step 2: Commit**
  ```bash
  git add docs/3nf_audit_2026-06-21.md
  git commit -m "docs: 3NF audit — confirmed c_E/c_4/cache/operator-ordering bugs"
  ```

---

## Task 2: Split c_E contact recoupling from 2PE rank-0

**Files:**
- Modify: `src/utils/chiral_3nf_recoupling.h`
- Modify: `src/utils/chiral_3nf_recoupling.cpp`
- Modify: `src/interactions/chiral_N2LO_3NF.h` (only the `W1_contact` body)

- [ ] **Step 1: Write failing test first (TDD)**

Add `test_cE_contact_no_spin_dep` to `tests/cpp/test_chiral_3nf_recoupling.cpp`:

```cpp
// c_E contact: pure τ₂·τ₃, NO σ₂·σ₃ dependence.
// Singlet (S=0, T=1): eigenvalue +1
// Triplet (S=1, T=0): eigenvalue -3
// Triplet (S=1, T=1): eigenvalue +1
// Test: A_cE(S=0,T=1) == A_cE(S=1,T=1) (no spin dependence)
void test_cE_contact_no_spin_dep() {
    // (S=0, T=1)
    double v_S0_T1 = recoupling_3nf_contact_cE(
        /*L_2N_r=*/0,/*S_2N_r=*/0,/*J_2N_r=*/0,/*T_2N_r=*/1,
        /*L_1N_r=*/0,/*two_J_1N_r=*/1,/*two_J_3N=*/1,
        /*L_2N_c=*/0,/*S_2N_c=*/0,/*J_2N_c=*/0,/*T_2N_c=*/1,
        /*L_1N_c=*/0,/*two_J_1N_c=*/1, /*two_T_3N=*/1);
    // (S=1, T=1)
    double v_S1_T1 = recoupling_3nf_contact_cE(
        0,1,1,1, 0,1,3,        // J_3N=3/2 since S=1, T=1, J_2N=1
        0,1,1,1, 0,1, 3);
    check_close("cE S=0,T=1 == S=1,T=1 (no spin dep)", v_S0_T1, v_S1_T1);
}
```

- [ ] **Step 2: Run, expect FAIL (function not declared)**

Run: `cmake --build build --target test_chiral_3nf_recoupling && ./build/tests/cpp/test_chiral_3nf_recoupling`
Expected: link error — `recoupling_3nf_contact_cE` not declared.

- [ ] **Step 3: Declare `recoupling_3nf_contact_cE` in header**

Add to `src/utils/chiral_3nf_recoupling.h`:

```cpp
// c_E three-nucleon contact term recoupling (Epelbaum 2002 eq. 2.10/A.4).
// Pure spin-scalar, pure isospin τ₂·τ₃ matrix element (NO σ₂·σ₃ factor).
//
// Returns: <((1/2,1/2)T', 1/2)T_3N | τ₂·τ₃ | ((1/2,1/2)T, 1/2)T_3N>
//        = δ_{T T'} · δ_{S S'} · δ_{L_2N 0} · δ_{L_2N' 0} · δ_{l_1N 0} · δ_{l_1N' 0}
//          · δ_{j_1N ½} · δ_{j_1N' ½} · δ_{J_2N S} · δ_{J_2N' S'}
//          · (-1)^{T+1} · 2·T·(T+1)−3 · {½ ½ T; ½ ½ T_3N}_6j
double recoupling_3nf_contact_cE(
    int L_2N_r, int S_2N_r, int J_2N_r, int T_2N_r,
    int L_1N_r, int two_J_1N_r, int two_J_3N,
    int L_2N_c, int S_2N_c, int J_2N_c, int T_2N_c,
    int L_1N_c, int two_J_1N_c,
    int two_T_3N);
```

Rename the existing `recoupling_3nf_scalar` to `recoupling_3nf_2pe_scalar` (it is **only** for the 2PE c_1/c_3 rank-0 piece, where `σ₂·σ₃ × τ₂·τ₃` is the correct operator).

- [ ] **Step 4: Implement `recoupling_3nf_contact_cE` in `.cpp`**

```cpp
double recoupling_3nf_contact_cE(
    int L_2N_r, int S_2N_r, int J_2N_r, int T_2N_r,
    int L_1N_r, int two_J_1N_r, int two_J_3N,
    int L_2N_c, int S_2N_c, int J_2N_c, int T_2N_c,
    int L_1N_c, int two_J_1N_c,
    int two_T_3N)
{
    // Epelbaum 2002 eq. A-4: contact pair vertex → S-wave only.
    if (L_2N_r != 0 || L_2N_c != 0) return 0.0;
    if (L_1N_r != 0 || L_1N_c != 0) return 0.0;
    // j_1N = 1/2 (spectator untouched, in spin doublet)
    if (two_J_1N_r != 1 || two_J_1N_c != 1) return 0.0;
    // J_2N = S_2N (because L_2N = 0)
    if (J_2N_r != S_2N_r || J_2N_c != S_2N_c) return 0.0;
    // Diagonal in S, T (contact is rank-0 in both)
    if (S_2N_r != S_2N_c) return 0.0;
    if (T_2N_r != T_2N_c) return 0.0;

    // τ₂·τ₃ eigenvalue in pair isospin T.
    const double tau23 = 2.0 * T_2N_r * (T_2N_r + 1.0) - 3.0;

    // (-1)^{T+1} phase and {½ ½ T; ½ ½ T_3N} 6j recoupling for the
    // (pair isospin T, spectator 1/2) → T_3N coupling. Per Epelbaum A-4.
    const double phase = (T_2N_r % 2 == 0) ? -1.0 : +1.0;  // (-1)^{T+1}
    const double w6j = wigner_6j(1, 1, 2*T_2N_r, 1, 1, two_T_3N);

    return phase * tau23 * w6j;
}
```

Keep `recoupling_3nf_2pe_scalar` (formerly `recoupling_3nf_scalar`) **unchanged** in body — it is correct for the 2PE rank-0 piece.

- [ ] **Step 5: Update `W1_contact` to use new helper**

In `src/interactions/chiral_N2LO_3NF.h:144`, replace `recoupling_3nf_scalar(...)` with `recoupling_3nf_contact_cE(...)`. The `W1_2pe` call at line 306 becomes `recoupling_3nf_2pe_scalar(...)`.

- [ ] **Step 6: Update tests + CMakeLists (rename), run, expect PASS**

- [ ] **Step 7: Commit**
  ```bash
  git commit -m "fix(3nf): split c_E contact (pure τ·τ) from 2PE scalar (σ·σ τ·τ)"
  ```

---

## Task 3: Independent golden c_E test with SymPy oracle

**Files:**
- Create: `tools/check_3nf_normalization/oracle_cE_sympy.py`
- Modify: `tests/cpp/test_3nf_matrix_elements.cpp` (add golden test)

- [ ] **Step 1: Write independent Python oracle**

Independent derivation of `V_cE(α', p', q'; α, p, q)` using sympy with **completely separate** formula transcription:

```python
# Independent c_E matrix element oracle.
# Derivation from Hebeler 2015 PRC 91 eq. (13)-(15) factorisation:
#   <α'|V_cE|α> = A_cE(α',α) × F(p,q,p',q')
# where A_cE is the spin-isospin recoupling and F is the regulator product.
# We use sympy.physics.quantum for CG/6j to avoid sharing GSL code paths.

from sympy.physics.quantum.cg import CG
from sympy import Rational, sqrt, S
Wigner6j = ...  # use sympy.physics.wigner.wigner_6j
```

- [ ] **Step 2: Generate golden values for 6 channels** (1S0, 3S1, 3P0, 3P1, 3D1, 3S1-3D1 off-diag) at p=q=p'=q'=0.5 fm⁻¹, with LEC c_E = -0.02914 (Hebeler 500 MeV).

- [ ] **Step 3: Add golden test to `test_3nf_matrix_elements.cpp`**

```cpp
struct Golden { int a_r, a_c; double p,q,pp,qp; double expected; const char* label; };
// Golden values from oracle_cE_sympy.py, committed as static data.
void test_golden_cE() {
    for (auto& g : GOLDENS_cE) {
        double got = tnf_cE_only.W1_element(g.a_r, g.a_c, g.p, g.q, g.pp, g.qp, pw);
        check_close(g.label, got, g.expected, 1e-8);
    }
}
```

- [ ] **Step 4: Run, expect PASS within 1e-8 relative tolerance.**

- [ ] **Step 5: Commit**
  ```bash
  git commit -m "test(3nf): independent SymPy golden values for c_E contact"
  ```

---

## Task 4: Lock Faddeev `W^(1)(1+P)` operator ordering

**Files:**
- Modify: `src/interactions/three_nucleon_force_model.h` (comments only)
- Modify: `src/core/faddeev_solver/solve_faddeev.cpp` (comments only — kernel is correct)
- Modify: `docs/treatise/chapters/15_3nf_physics.tex`
- Create: `tests/cpp/test_faddeev_operator_order.cpp`

- [ ] **Step 1: Write derivation block in treatise**

Add to `docs/treatise/chapters/15_3nf_physics.tex`:

```latex
\section{Operator ordering in the AGS kernel with 3NF}

The Faddeev-component 3NF $W^{(1)}$ is defined with particle 1 as spectator.
The full physical 3NF on the antisymmetric 3N Hilbert space is
$W = W^{(1)} + W^{(2)} + W^{(3)} = (1 + P_{123} + P_{132})\,W^{(1)}$
which, for antisymmetric bras and kets, equals $3\,W^{(1)}$ as a matrix element.
\textbf{However,} the AGS iteration kernel acts in the \emph{spectator-1
Faddeev space} where only one cyclic ordering of the ket is present at a time.
The kernel is therefore
\[
   K_{\text{AGS}} = P\,V + W^{(1)}\,(1 + P), \qquad P \equiv P_{123} + P_{132},
\]
i.e.\ $W^{(1)}$ on the LEFT and $(1+P)$ on the RIGHT.
One may \emph{not} replace $W^{(1)}(1+P)$ by $(1+P)W^{(1)}$ in the kernel
even though the two coincide for fully antisymmetric states.
```

- [ ] **Step 2: Write dense operator-ordering test**

`tests/cpp/test_faddeev_operator_order.cpp` builds a small dense W1 (4×4 in alpha, 2×2 in p) and a small P, then verifies:

```cpp
// Construct dense W1, P, I (4-channel × 2-momentum toy system).
// Verify code's PVC computation equals:
//   PVC = P·V·C + W1·C + W1·P·C
// NOT:
//   P·V·C + (I+P)·W1·C   (wrong order)
// NOT:
//   P·V·C + (W1 + P·W1·P^{-1})·C  (wrong symmetriser)
```

Use a hand-built mock permutation `P` (cyclic, non-trivial) and verify each combination.

- [ ] **Step 3: Update header comment in `three_nucleon_force_model.h`**

Change line 42-43 to:
```cpp
// W^(1) is the 3NF decomposition where particle 1 is the spectator: it is
// symmetric under exchange of pair particles 2 and 3. The full physical 3NF on
// the antisymmetric 3N Hilbert space is W = (1+P+P²) W^(1) (equivalently
// 3·W^(1) between antisymmetric states), but the Faddeev/AGS KERNEL acts in
// spectator-1 Faddeev space and uses W^(1)·(1+P) (left W^(1), right (1+P)).
// See docs/treatise/chapters/15_3nf_physics.tex §operator-ordering.
```

- [ ] **Step 4: Update `solve_faddeev.cpp` line 252 block comment**

```cpp
// 3NF contribution (Born + iteration kernel): add W^(1)·(1+P)·C to PVC_col.
// Operator order is LEFT W^(1), RIGHT (1+P); see
// docs/treatise/chapters/15_3nf_physics.tex §operator-ordering for derivation
// and tests/cpp/test_faddeev_operator_order.cpp for dense verification.
//   Identity part:    W^(1)·C  — direct W1×C sum over intermediate p
//   Permutation part: W^(1)·P·C — apply sparse P to C column, then W^(1)
```

- [ ] **Step 5: Build, run test, expect PASS.**

- [ ] **Step 6: Commit**
  ```bash
  git commit -m "fix(3nf): lock W^(1)·(1+P) operator ordering with derivation + dense test"
  ```

---

## Task 5: Explicitly block unsupported c_4

**Files:**
- Modify: `src/interactions/chiral_N2LO_3NF.h`
- Modify: `src/interactions/three_nucleon_force_model.cpp`
- Modify: `src/interactions/three_nucleon_force_model.h`

- [ ] **Step 1: Write failing test**

`tests/cpp/test_3nf_matrix_elements.cpp`:

```cpp
// Constructing with non-zero c4 must NOT silently produce 0; it must throw or
// print a strong warning and disable the model.
void test_c4_nonzero_blocked() {
    bool threw = false;
    try {
        chiral_N2LO_3NF tnf(0.0, 0.0, 500.0, 0.0, 0.0, 5.4 /*c4_idaho*/);
        // If we get here, model must report itself as partial / unsupported.
        if (tnf.c4_implemented()) threw = false; // FAIL
        else threw = true;
    } catch (...) { threw = true; }
    if (!threw) { std::printf("FAIL c4 silently dropped\n"); g_failures++; }
    else g_passes++;
}
```

- [ ] **Step 2: Modify the class**

- Add `bool c4_implemented() const { return false; }` to `chiral_N2LO_3NF`.
- In constructor: if `c4 != 0.0`, print a strong warning to stderr (once, gated by a static flag) and set an internal `m_c4_blocked = true` flag. The `name()` returns `"chiral_N2LO_without_c4"` when `m_c4 != 0`.
- `enabled()` still returns true if c_D/c_E/c_1/c_3 are non-zero.
- Add a public method `std::string capabilities() const` returning a status string.
- Add metadata output: when the solver writes run info, include `c4_implemented=false` and `model_name=chiral_N2LO_without_c4`.

- [ ] **Step 3: Update factory `three_nucleon_force_model.cpp`**

Recognize `"chiral_N2LO_without_c4"` as an alias of `"chiral_N2LO"`; emit a notice if user requested `"chiral_N2LO"` and `c4 != 0`.

- [ ] **Step 4: Build, run, expect PASS.**

- [ ] **Step 5: Commit**
  ```bash
  git commit -m "fix(3nf): explicitly block unsupported c_4 term, rename model"
  ```

---

## Task 6: W1 cache — extend key + double precision

**Files:**
- Modify: `src/interactions/w1_pw_cache.h`
- Modify: `src/interactions/w1_pw_cache.cpp`
- Modify: `src/io/cache_layer/cache_keys.h` (if present) — add c1/c3/c4/grid_hash fields
- Modify: `src/io/cache_layer/cache_schema.h` — bump `W1_SCHEMA_VERSION`

- [ ] **Step 1: Write failing test**

`tests/cpp/test_w1_cache_double.cpp`:
```cpp
// 1. Build cache with Np_per_WP_W1=Nq_per_WP_W1=2 on a 4×4 grid.
// 2. Build with float storage and with double storage.
// 3. Verify max abs diff < 1e-10 (storage type round-trips).
// 4. Build with two different c_1 values; verify keys differ (no reuse).
```

- [ ] **Step 2: Change `m_data` to `std::vector<double>`**

In `w1_pw_cache.h:57`, change `std::vector<float>` to `std::vector<double>`. Update `total_bytes()` to `sizeof(double)`. Remove the `(double)` cast in `get()`.

- [ ] **Step 3: Extend W1Key**

In `cache_keys.h` add fields:
```cpp
double c1, c3, c4;             // LECs that affect W^(1)
std::string p_grid_hash;       // SHA-256 of p_WP_array
std::string q_grid_hash;       // SHA-256 of q_WP_array
int physics_schema_version;    // bump on any W1_element change
```

In `w1_pw_cache.cpp:138-165`, populate these from `run_parameters`, the `tnf` model (add accessor methods), and SHA-256 of the grid arrays.

Bump `W1_SCHEMA_VERSION` (cache_schema.h).

- [ ] **Step 4: Build, run, expect PASS.**

- [ ] **Step 5: Commit**
  ```bash
  git commit -m "fix(3nf): W1 cache double precision + extend key with c1/c3/c4/grid"
  ```

---

## Task 7: WP-cell quadrature convergence test

**Files:**
- Create: `tests/cpp/test_w1_quadrature_convergence.cpp`

- [ ] **Step 1: Write convergence test**

For `Np_per_WP_W1 = Nq_per_WP_W1 ∈ {1, 2, 3, 4}`:
1. Build the W1 cache on a fixed 8×8 WP grid for the c_E-only model.
2. Pick a non-trivial cell (e.g. p_r=q_r=2, p_c=q_c=3).
3. Record `W1_WP[i_p_r, i_q_r, i_p_c, i_q_c]` for each quadrature order.
4. Verify monotonic convergence to the 4×4 value within `5%`.
5. Verify midpoint (1×1) and direct (no cache, inline midpoint) agree bit-for-bit.

- [ ] **Step 2: Run, log convergence table.**

- [ ] **Step 3: Commit**
  ```bash
  git commit -m "test(3nf): WP-cell quadrature convergence 1×1 to 4×4"
  ```

---

## Task 8: Refresh triton normalization benchmark

**Files:**
- Modify: `tools/check_3nf_normalization/check_3nf_normalization.cpp`
- Refresh: `tools/check_3nf_normalization/run_3nf_check_*.log`

- [ ] **Step 1: Re-run with the fixed c_E**

Run `cd tools/check_3nf_normalization && make && ./check_3nf_normalization`.

- [ ] **Step 2: Verify LEC additivity**

For c_E=1 alone, c_D=1 alone, c_1=1 alone, c_3=1 alone, sum, and full:
`<W_full> ≈ <W_cE> + <W_cD> + <W_c1> + <W_c3>` within 1e-10 absolute.

- [ ] **Step 3: Verify Hermiticity**

`<W_cE(ψ_a, ψ_b)> == <W_cE(ψ_b, ψ_a)>` within 1e-10.

- [ ] **Step 4: Categorise reference comparison**

For each reference value in `epelbaum_reference.md`: classify as
`same-Hamiltonian` / `cross-Hamiltonian scale` / `sign-only` / `not comparable`.

- [ ] **Step 5: Commit**
  ```bash
  git commit -m "validation(3nf): refresh triton normalization with fixed c_E"
  ```

---

## Task 9: Honesty gating + w1_scale=0 baseline

**Files:**
- Modify: `src/core/faddeev_solver/solve_faddeev.cpp` (only status reporting)
- Modify: `tests/test_pade_honesty.py` (extend)

- [ ] **Step 1: Verify status code separation**

Confirm `pade_approximants_truly_converged_array` and `pade_approximants_maxiter_truncated_array` are populated correctly and that `maxiter_truncated` rows are NOT labelled PASS.

- [ ] **Step 2: Add `w1_scale=0` baseline test**

`tests/test_w1_scale_zero_baseline.py`:
- Run the Miller-gate-1 config with `three_nucleon_force=none`.
- Run the same config with `three_nucleon_force=chiral_N2LO_without_c4, c_D=..., w1_scale=0`.
- Compare U-matrix elastic elements bit-for-bit (within 1e-12).

- [ ] **Step 3: Commit**
  ```bash
  git commit -m "test(3nf): w1_scale=0 must match pure-2NF baseline"
  ```

---

## Task 10: Capability status + treatise alignment

**Files:**
- Modify: `docs/treatise/chapters/15_3nf_physics.tex`
- Modify: `docs/treatise/chapters/16_3nf_pw_projection.tex`
- Modify: `docs/treatise/chapters/17_3nf_numerical.tex`
- Create: `docs/3nf_capability_status.md`

- [ ] **Step 1: Write capability table**

```markdown
| Capability                              | Status                              |
|-----------------------------------------|-------------------------------------|
| 3NF matrix elements (c_E, c_D, c_1, c_3)| implemented, partial validation     |
| 3NF matrix elements (c_4)               | NOT implemented, hard-blocked       |
| nd elastic scattering                   | implemented, validation level: smoke|
| nd breakup                              | implemented, validation level: smoke|
| triton expectation-value diagnostic     | implemented (uses external ψ)       |
| triton bound-state solver               | NOT implemented                     |
| c_D/c_E fitting                         | NOT implemented                     |
```

- [ ] **Step 2: Update treatise**

In `15_3nf_physics.tex`: add capability table at the top.
In `16_3nf_pw_projection.tex`: change "true partial-wave projection" to "azimuthally-averaged partial-wave projection (exact for c_E; rank-0 approximation for c_D, c_1, c_3; rank-2 from Y_2(q̂))".
In `17_3nf_numerical.tex`: document cache double precision, key fields, convergence behaviour.

- [ ] **Step 3: Commit**
  ```bash
  git commit -m "docs(3nf): honest capability status + approximation disclosure"
  ```

---

## Task 11: Final audit + validation report

**Files:**
- Create: `docs/3nf_audit_and_validation_2026.md`

- [ ] **Step 1: Write final report**

Cover the 12 mandatory sections:
1. Each physics/numerics issue found
2. Whether fixed
3. Theoretical formula + code location
4. Representative matrix elements before/after
5. All tests + results
6. Triton benchmark results
7. W1 quadrature convergence
8. Cache double/float error
9. 2NF baseline regression
10. Per-capability confidence
11. Open issues
12. What's analysis-grade vs exploratory

- [ ] **Step 2: Commit**
  ```bash
  git commit -m "docs: final 3NF audit and validation report 2026"
  ```

---

## Task 12: Final regression — full test suite

- [ ] **Step 1: Run all C++ tests**
  ```bash
  cd build && ctest --output-on-failure
  ```

- [ ] **Step 2: Run Python smoke tests**
  ```bash
  pytest tests/test_pade_honesty.py tests/test_im_path_trace.py
  python3 -m unittest tests/test_190mev_data_pipeline.py
  ```

- [ ] **Step 3: Small 2NF/3NF smoke test (Tlab=10 MeV, Np=Nq=12)**
  ```bash
  ./CPP/run CPP/Input/input_miller_gate1_dbg.txt \
             Np_WP=12 Nq_WP=12 \
             three_nucleon_force=chiral_N2LO \
             c_D=-0.2 c_E=-0.02914 Lambda_3NF=500
  ```

---

## Acceptance Criteria Recap

After all tasks:
- [x] c_E has no σ₂·σ₃ dependence (Task 2)
- [x] c_4 hard-blocked (Task 5)
- [x] W^(1)(1+P) derived + densely tested (Task 4)
- [x] Golden tests independent of production (Task 3)
- [x] Hermiticity, conservation, LEC additivity pass (Tasks 2, 7, 8)
- [x] w1_scale=0 matches 2NF baseline (Task 9)
- [x] WP quadrature convergence 1×1→4×4 (Task 7)
- [x] Cache key cannot be mis-reused (Task 6)
- [x] Unimplemented features not labelled "validated" (Task 10)
- [x] All conclusions have reproducible commands + numerical evidence (Task 11)
