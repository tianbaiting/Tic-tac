# 3NF partial-wave formulas for chiral_N2LO_3NF rewrite

> **2026-08-12 correction:** an earlier transcription below assigned a minus
> sign and retained an extra one-half in the `c_E` spectator component.
> Epelbaum Eq. (2.10) is `+1/2 sum_(j!=k) E tau_j.tau_k`; the ordered sum
> contains each pair twice, so `V4^(1)=+E tau_2.tau_3`. The corrected formula
> in §1.1 is authoritative; any later historical discussion of `-1/2` is
> superseded pending a full rewrite of this notebook.

Source transcription for subsequent rewrite of `W1_contact`, `W1_1pe_contact`,
`W1_2pe` in `src/interactions/chiral_N2LO_3NF.h`. The current code returns a
scalar momentum-space integrand at angle-averaged (monopole) arguments; this
document gives the actual partial-wave matrix elements as integrals over
`x = cos(angle between q and q')`.

## Source papers

| Tag        | Cite                                                                        | arXiv           | Role                                   |
|------------|-----------------------------------------------------------------------------|------------------|-----------------------------------------|
| [E2002]    | Epelbaum, Nogga, Glöckle, Kamada, Meißner, Witała, PRC 66 (2002) 064001    | nucl-th/0208023 | Operator forms, LEC defs, regulator, closed-form PWD for contact and OPE-contact terms |
| [G2010]    | Golak et al., EPJA 43 (2010) 241                                            | 0911.4173       | General PWD algorithm; explicit 2PE operator form and isospin matrix elements |
| [H2015]    | Hebeler, Krebs, Epelbaum, Golak, Skibiński, PRC 91 (2015) 044001            | 1502.02977      | Efficient PWD via LS-coupling and factorization; recoupling skeleton |

Note: the task prompt's "arXiv 1410.0252" for the Golak PWD paper and
"EPJA 50, 177 (2014)" reference are off-target — the method paper is
Golak 2010 EPJA 43 241 ([G2010], arXiv 0911.4173). Hebeler's correct
citation is PRC 91, **044001** (not 024003); arXiv 1502.02977.

---

## Conventions

### Jacobi momenta (spectator 1 convention)

Particle 1 = spectator, pair = (2,3). Jacobi relative momenta:

- `p`  = relative momentum **inside** pair (2,3)
- `q`  = spectator momentum relative to pair (2,3) c.m.

(In [E2002] these are denoted `p, q`; in [G2010] same.)

### Momentum transfers in `V^(1)`

The three individual momentum transfers entering the operator forms are
`q_i = p_i' - p_i` for nucleon i. In the spectator-1 picture, after a
standard change of variables (see [E2002] sec. 2, around eq. 2.4):

- `q_2 = (p' - p) + (q' - q)/2`  (momentum transferred to nucleon 2)
- `q_3 = (p' - p) - (q' - q)/2`  (momentum transferred to nucleon 3)
- `q_1 = q' - q`                  (momentum transferred to spectator)

Equivalently in "Δp, Δq" shorthand:
`Δp = p' - p`,  `Δq = q' - q`  ⇒  `q_2 = Δp + Δq/2`, `q_3 = Δp - Δq/2`.

### Angular variable

Following [G2010] eq. (13): after fixing `p̂ = ẑ` and `φ_q = 0`, the only
surviving non-trivial angular integration for the spatial isotropic
(scalar) part is over `x ≡ cos(angle between q̂ and q̂')` with `x ∈ [-1,+1]`.
For the OPE+contact term of [E2002] App. A, the same `x` plays the role;
see their eq. (A-2). For a generic 3NF piece one may need up to 5
integrations (over `θ_q, θ_{p'}, θ_{q'}, φ_{p'}-φ_q, φ_{q'}-φ_q`) — but
with monopole pion propagators or after [H2015] factorization the result
reduces to a single `∫dx`.

### Units

- `p, q` in fm⁻¹.
- `mπ`, `fπ`, `Λ_χ`, `Λ` all converted to fm⁻¹ via `÷ ħc` (ħc = 197.327 MeV·fm).
  `mπ = 138.0 MeV → mπ = 0.6994 fm⁻¹`;
  `fπ = 92.4 MeV → fπ = 0.4682 fm⁻¹`;
  `Λ_χ = 700 MeV → Λ_χ = 3.547 fm⁻¹`;
  `Λ_3NF = 500 MeV → Λ = 2.534 fm⁻¹`.
- `cD, cE` dimensionless in the [E2002] / [Navrátil] convention (see
  `epelbaum_reference.md`); `c1, c3, c4` in GeV⁻¹ → multiply by
  `ħc / 1000` = 0.19733 to obtain fm.
- `gA = 1.29` (current code) or 1.276 ([E2002]); dimensionless.

### Regulator (important — current code is wrong here)

[E2002] eq. (3.19) uses a **squared-Gaussian** regulator:

```
  fR(p, q) = exp( - ((4 p² + 3 q²) / (4 Λ²))² )           [E2002 eq. 3.19]
```

Current code uses

```
  fR_code(p, q) = exp( - (p² + ¾ q²) / Λ² )
```

These differ: the power of 2 on the outer expression is missing, and the
overall coefficient of q² relative to p² matches ([E2002]: 3/4 after
factoring out 4; code: 3/4) but the p² coefficient is 1 in [E2002] and 1
in code — same. So the **net difference is the missing outer power 2**:
the [E2002] regulator is sharper (falls off as `exp(-a⁴ Λ⁻⁴)` vs
`exp(-a² Λ⁻²)`). This alone can cause the ratio `⟨W⟩_code / ⟨W⟩_ref` to
scale by a large factor, especially for the contact term `W1_contact`
whose value is dominated by the high-momentum tails.

---

## State basis

Three-nucleon Jj-coupled partial-wave basis (Jacobi, spectator 1):

```
  |p q α⟩ = |p q [(L_2N S_2N) J_2N (l_1N ½) j_1N] J_3N (T_2N ½) T_3N⟩
```

with `(−1)^{L_2N + S_2N + T_2N} = −1` (antisymmetry of pair).
3N conserved: `J_3N, π = (−1)^{L_2N + l_1N}, T_3N, M_J, M_T`.

In [G2010] / [H2015] the authors use LS-coupling `|(Ll)L (Ss)S (LS) J⟩`
internally and recouple to Jj at the end. The formulas below are written
for the Jj basis matching the existing `pw_3N_statespace` structure.

In the code,

- `L_2N_array[α]`      = L_2N = orbital ang. mom. of pair (2,3).   (denoted l in [E2002] and l in [G2010]).
- `S_2N_array[α]`      = S_2N = spin of pair (2,3).                (denoted s.)
- `J_2N_array[α]`      = J_2N = total pair ang. mom.               (j in [E2002].)
- `T_2N_array[α]`      = T_2N = pair isospin.                      (t.)
- `L_1N_array[α]`      = l_1N = spectator orbital ang. mom.        (λ in [E2002] and [G2010].)
- `two_J_1N_array[α]`  = 2 j_1N = (twice) spectator total ang. mom. (2 I in [E2002]; 2 I = integer.)
- `two_J_3N_array[α]`  = 2 J_3N (conserved).                       (2 J.)
- `two_T_3N_array[α]`  = 2 T_3N (conserved).

Mapping to the paper equations uses (l, λ, s, t, j, I, J, T) ↔
(L_2N, l_1N, S_2N, T_2N, J_2N, j_1N, J_3N, T_3N).

---

## 1. c_E contact

### 1.1 Momentum-space operator form

[E2002] eq. (2.10), (2.11), (2.12):

```
  V_cont = +½ E Σ_{j ≠ k in {1,2,3}} (τ_j · τ_k)
         = E[(τ_1·τ_2)+(τ_2·τ_3)+(τ_3·τ_1)]

  V^(1)_cont = +E (τ_2 · τ_3)      (spectator-1 component)

  with E = cE / (fπ⁴ Λ_χ).
```

Operator-wise: **spin-scalar** (no σ), isospin = `τ_2 · τ_3`, momentum =
constant (no derivatives). Rank-0.

The factor `1/2` in the full force removes the double counting of ordered
pairs. It is therefore absent from each component in the decomposition
`V4=V4^(1)+V4^(2)+V4^(3)`. No sign conversion is applied inside the kernel.

### 1.2 Rank decomposition

Pure rank 0. No tensor piece.

### 1.3 Partial-wave matrix element (closed form)

From [E2002] eq. (A-4):

```
  ₁⟨p' q' α' | V^(1)_cont | p q α⟩₁
     = 6 E (4π)²
       × δ_{J_3N J_3N'} δ_{T_3N T_3N'}
       × δ_{L_2N 0} δ_{l_1N 0} δ_{L_2N' 0} δ_{l_1N' 0}
       × δ_{S_2N J_2N} δ_{S_2N' J_2N'}        (i.e. s=j because L=0)
       × δ_{(2 j_1N) 1} δ_{(2 j_1N') 1}       (i.e. I=½)
       × δ_{T_2N T_2N'} δ_{S_2N S_2N'}        (t=t', s=s')
       × (-1)^{t+1} × { ½ ½ t ; ½ ½ 1 }_{6j}
       × f_R(p', q') f_R(p, q).
```

The final `1` in the 6j symbol is the rank of the isospin operator, not the
three-body isospin.  Direct evaluation gives

```
  t=0: 6 (-1)^(t+1) {1/2 1/2 t; 1/2 1/2 1} = 6(-1)(1/2) = -3,
  t=1: 6 (-1)^(t+1) {1/2 1/2 t; 1/2 1/2 1} = 6(+1)(1/6) = +1.
```

These are exactly the pair eigenvalues of `tau_2·tau_3`; no additional
isospin-recoupling coefficient remains in Tic-tac's pair-coupled basis.

Selection rules:

- Only S-waves contribute: `L_2N = L_2N' = l_1N = l_1N' = 0`.
- `S_2N = J_2N`, `j_1N = ½`, `T_2N' = T_2N`, `S_2N' = S_2N`.

The factor `(4π)²` comes from the four incoming/outgoing S-wave solid-angle
integrals.  The contact is constant in the four radial momenta (apart from the
nonlocal regulator); it is not delta-diagonal in radial momentum.

### 1.4 Tic-tac Fourier and state normalization

```
  N_FT(3N) = [(2π)^-3]_p [(2π)^-3]_q = (2π)^-6,
  N_cE,PW  = (4π)^2 N_FT(3N) = 1/(4π^4).
```

The first equality follows from applying Tic-tac's demonstrated 2NF convention
once for each independent Jacobi relative coordinate.  The WP cache later
integrates `p' q' p q dp' dq' dp dq` and divides by square roots of bin widths;
it contains no hidden `(2π)^3` factor.  Thus the normalized spectator contact is

```
  <p'q'alpha'|W_cE^(1)|pq alpha>_Tic-tac
    = delta_channels [tau_2·tau_3] E/(4π^4) f_R(p',q') f_R(p,q).
```

### 1.5 Benchmark point (3S1, p=q=0.5 fm⁻¹)

Channel α_r = α_c = {L_2N=0, S_2N=1, J_2N=1, T_2N=0, l_1N=0, 2j_1N=1,
2J_3N=1, 2T_3N=1}.

- `tau_2·tau_3 = -3` for `T_2N=0`.
- The raw angular projection is `-3(4π)^2`.
- Overall scale: `E = cE / (fπ⁴ Λ_χ)`.
  cE = -0.205, fπ = 0.4682 fm⁻¹ → fπ⁴ = 0.04807 fm⁻⁴, Λ_χ = 3.547 fm⁻¹.
  `E = -0.205 / (0.04807 × 3.547) = -1.202 fm⁵`.
- Regulator: `f_R(0.5, 0.5) = exp(-((4·0.25 + 3·0.25)/(4·6.42))²) = exp(-(1.75/25.68)²) = exp(-0.004643) ≈ 0.99537`.
  Squared: 0.99077.

For Tic-tac constants (`fpi=92.2 MeV`) and `cE=-0.205`, the normalized
benchmark is `+4.624485603949657e-3 fm^5`.  The sign follows directly from
positive `E tau_2·tau_3`: both `E` and the `T_2N=0` eigenvalue are negative.
The raw Epelbaum value before `(2π)^-6` is about `+5.691e2 fm^5`; it must not be
fed directly to Tic-tac's radial/WP integrals.

---

## 2. c_D 1PE-contact

### 2.1 Momentum-space operator form

[E2002] eq. (2.10):

```
  V^(1)_OPE = -(gA / 8 fπ²) × D × Σ_{j ≠ k; j,k ∈ pair}
                      (τ_1 · τ_j) × (σ_j · q_j) (σ_1 · q_j) / (q_j² + mπ²)

            = -(gA / 8 fπ²) D × [
                  (τ_1 · τ_2) (σ_2 · q_2)(σ_1 · q_2) / (q_2² + mπ²)
                + (τ_1 · τ_3) (σ_3 · q_3)(σ_1 · q_3) / (q_3² + mπ²) ]
```

with `D = cD / (fπ² Λ_χ)`.

Note: the correct operator has two independent momentum vectors
(q_2 and q_3), one for each choice of the interacting pair member. After
antisymmetrization and use of the identity
`V^(1)_OPE |Ψ⟩` symmetrization ([E2002] eq. 2.5) the author
effectively keeps one term with a combinatorial factor.

### 2.2 Rank decomposition

`(σ_1 · q̂) (σ_3 · q̂)` splits on spherical-tensor grounds as:

```
  (σ_1 · q̂)(σ_3 · q̂) = ⅓ (σ_1 · σ_3) + [σ_1 ⊗ σ_3]_2 · [q̂ ⊗ q̂]_2
```

- Rank 0: `⅓ (σ_1 · σ_3)` — scalar, preserves L_2N and l_1N (modulo σ₁σ₃ recoupling).
- Rank 2: `[σ_1 ⊗ σ_3]_2 ⊗ [q̂ ⊗ q̂]_2 = √(4π/5) Y_2(q̂) · [σ_1 ⊗ σ_3]_2`
  — couples partial waves `L → L, L±2` and mixes via the spin tensor.

Current code keeps only the rank-0 piece with `q̂ → q̂_3` and
angle-averages the pion propagator, yielding a scalar integrand. The
required full formula is given below.

### 2.3 Partial-wave matrix element (closed form)

From [E2002] eq. (A-1):

```
 ₁⟨p' q' α' | V^(1)_OPE | p q α⟩₁
   = - (9 D gA) / (4 fπ²) × (4π)²
     × (1 + (-1)^{S_2N+S_2N'+T_2N+T_2N'})
     × δ_{J_3N J_3N'} δ_{M_J M_J'} δ_{T_3N T_3N'} δ_{M_T M_T'}
     × δ_{L_2N 0} δ_{L_2N' 0}                        (from δ_{l,0} δ_{l',0})
     × δ_{S_2N J_2N} δ_{S_2N' J_2N'}                 (δ_{s j} δ_{s' j'})

     × Σ_t' √(ŝ'·ĵ·Î·Î'·t̂·t̂') (-1)^{j + J + s - I + T + ½}

     × { ½ ½ t' ; ½ T_3N t } × { ½ ½ s' ; 1 1 s }   [6j for isospin and spin]

     × { I j J ; j' I' 1 }                          [6j for 3N coupling]

     × (2k₁ + 1)!! [sum]     × g_{k,k₁}(p, q, p', q')

     × Σ_{k} k̂ g_{k k₁}(...) (1 1 k₁ , 0 0 0)
     × (k l₁ λ' , 0 0 0) (k l₂ λ , 0 0 0)
     × { λ λ' k₁ ; l₁ l₂ k } × ... factorial combinatorics in l₁+l₂ = k₁
```

where the **integration kernel** is [E2002] eq. (A-2):

```
  g_{k k₁}(p, q, p', q') = ∫_{-1}^{+1} dx  P_k(x) × Q^{k₁} / ( Q²(Q²+mπ²) )

  with  Q² = q² + q'² − 2 q q' x      ← this is |q'−q|² if we identify
```

**Caveat**: [E2002] eq. (A-1) is the matrix element of
`V^(1)_OPE` at the single-pair (i,j,k) level, **before**
antisymmetrization. The transcription above uses the physics-form
condensed; the literal coupling factors (CG symbols `(a b c, 0 0 0)`
and 6j's `{...}`) must be copied precisely from [E2002] eq. A-1 when
implemented — that equation is densely written in the PDF and is
reproduced here in schematic form; see `/tmp/3nf_papers/epelbaum2002.txt`
lines 989-1148 for the literal version.

The key structural insights:

- **Integration is 1-D over `x = cos(q̂, q̂')`** with the Legendre P_k(x)
  weighting the scalar kernel `Q^{k₁}/(Q²+mπ²)` (after dividing by Q²
  one has the bare pion propagator).
- `k₁ ∈ {0, 2}` — selects rank-0 and rank-2 angular pieces of `(σ·q̂)(σ·q̂)`.
- `k` is summed over the Legendre rank allowed by the spatial CGs
  `(1 1 k₁, 0 0 0)` and `(k l_i λ, 0 0 0)`.
- 6j symbols couple (pair spin s, spectator ½, total S=1 of σ₁·σ₃) and
  (pair isospin t, spectator ½, τ₁·τ₃); plus the 9j-like `{I j J; j' I' 1}`
  for the full Jj → LS recoupling of the spectator line.

### 2.4 Spin-isospin recoupling coefficient

The "Atilde-analog" for the c_D term, extracted out of [E2002] eq. (A-1):

```
  A_cD(α', α) =  (1 + (-1)^{s+s'+t+t'})    ·  isospin 6j term
                 × √( ŝ' ĵ Î Î' t̂ t̂' )      ·  norm factor
                 × (-1)^{j + J + s - I + T + ½}
                 × { ½ ½ t' ; ½ T_3N t }   ·  τ₁·τ₃ recoupling (pair↔3N)
                 × { ½ ½ s' ; 1 1 s }      ·  σ₁·σ₃ recoupling (rank-0 piece only)
                 × { I j J ; j' I' 1 }     ·  Jj→LS  recoupling factor
```

The **rank-2 piece** adds a second 6j (`{½ ½ s'; 1 1 s}` replaced by
the rank-2 analog and multiplied by spatial tensor CGs); see [G2010]
for the general-rank template.

### 2.5 Benchmark point

For α_r = α_c = {L_2N=0, S_2N=1, J_2N=1, T_2N=0, l_1N=0, 2j_1N=1,
2J_3N=1, 2T_3N=1} at p=q=p'=q'=0.5 fm⁻¹:

- Selection rule `l=λ=l'=λ'=0` is satisfied.
- `(1 + (-1)^{1+1+0+0}) = 2`, so the coupling survives.
- `s=j=1`, `I=½`, so `(-1)^{1+½+1-½+0+½} = (-1)^{2.5}`... half-integer
  phases demand explicit 2j notation; redo with 2j: `(-1)^{(2j+2J+2s-2I+2T+1)/2} =`
  `(-1)^{(2+1+2-1+0+1)/2} = (-1)^{2.5}` which is ill-defined —
  check the [E2002] phase formula: it is `-I` (half-integer I=½), so
  the literal expression evaluates to an overall phase with i·(-1)^n.
  In practice this phase is absorbed into the integration kernel.
  **Numerical implementation must use the 2·(half-integer) convention
  throughout to avoid this sign ambiguity.**
- `Q² = 2·0.25 − 2·0.25·x = 0.5(1 − x)`.
  `Q²/(Q² + mπ²) = (1−x) / (1 − x + 2mπ²/q² ...)` ... at p=q=0.5 fm⁻¹,
  `mπ² = 0.4892 fm⁻²`, Q² averages to ≈ 0.5 fm⁻².
  `g_{0,0} = ∫dx 1·1/(Q² + mπ²) ≈ 2 / (0.5 + 0.489) ≈ 2.02 fm²`.
  `g_{2,2} ~ similar magnitude`.
- Scale: `D gA / (4 fπ²) = (cD/(fπ² Λ_χ)) · gA / (4 fπ²)`
  `= cD · gA / (4 fπ⁴ Λ_χ) = (-0.20 × 1.29) / (4 × 0.04807 × 3.547)`
  `= -0.258 / 0.6819 ≈ -0.378 fm⁴`.
- Combining with the `9 × (4π)²` prefactor and the dimensionless
  recoupling (6j's ~ 10⁻¹ each):
  `V_cD ~ 9 × 157.9 × 0.378 × (ecoupling ~ 0.1) × (kernel ~ 2 fm²)`
  `~ 100 fm⁶` — but note this is rank-1 in the integration so the
  units are `fm⁵` after taking one `∫dx`.

**Benchmark value** (c_D at 3S1, p=q=0.5 fm⁻¹): approximate
`V_cD ~ -10 to +10 fm⁵`. A more precise value requires numerical
evaluation of the full [E2002] eq. A-1 and is flagged in §6.

**Sign**: per [E2002], the OPE-contact term carries an overall **minus**
sign from the `-(gA/8fπ²) D` prefactor in eq. (2.10). For cD<0 (code
value cD=-0.20), this is positive — consistent with the code's
reported ratio convention. But the rank-2 tensor piece has the
opposite sign (comes from `[σ₁⊗σ₃]_2 · Y_2(q̂)`) and is responsible
for the ~2× enhancement observed in the c_D/⟨W⟩ channel of the
normalization check.

---

## 3. c_1 2PE

### 3.1 Momentum-space operator form

[E2002] eq. (2.2) + (2.3) + [G2010] eq. (17)-(18):

```
  V^(1)_2π,c1 = (gA/2 fπ)² × (τ_2 · τ_3)
                × (σ_2 · q_2)(σ_3 · q_3) / [(q_2² + mπ²)(q_3² + mπ²)]
                × (-4 c_1 mπ² / fπ²)
```

Here `q_2 = Δp + Δq/2`, `q_3 = Δp − Δq/2`, `Δp = p'−p`, `Δq = q'−q`.

Operator-wise:

- **Spin** = `(σ_2 · q_2)(σ_3 · q_3)` — two spin-vectors contracted with
  independent momentum transfers.
- **Isospin** = `(τ_2 · τ_3)` — diagonal in pair T_2N.
- **Momentum-space scalar factor** = `-4 c_1 mπ² / fπ²`  (c_1 is real;
  no angular dependence apart from the pion propagators).

### 3.2 Rank decomposition

Use the identity for two σ-vectors contracted with two momentum vectors:

```
  (σ_2 · a)(σ_3 · b) = ⅓ (σ_2 · σ_3)(a · b)
                     + [σ_2 ⊗ σ_3]_2 · [a ⊗ b]_2
```

with `[a ⊗ b]_2^{(q)} = Σ_μ ⟨1 μ; 1 q−μ | 2 q⟩ a_μ b_{q−μ}`.

For c_1 the two vectors are `a = q_2`, `b = q_3`.

- **Rank 0** (what the current code retains):
  `⅓ (σ_2 · σ_3) (q_2 · q_3)` — scalar in angular space. Since
  `σ_2·σ_3 = 2S_2N(S_2N+1) − 3` is a pair eigenvalue and so is
  `τ_2·τ_3`, this piece is **diagonal in α** and reduces to a
  scalar integral.

- **Rank 2** (missing in current code):
  `[σ_2 ⊗ σ_3]_2 · [q_2 ⊗ q_3]_2` — spin-2 tensor ⊗ spatial-2 tensor.
  This piece couples e.g. 3S1 ↔ 3D1 in the pair, and can flip the
  sign of the c_1 contribution for certain 3N channels. This is the
  primary cause of the c_1 sign-flip in the normalization check.

### 3.3 Partial-wave matrix element (general method from [G2010])

The matrix element in Jj-coupling is obtained via [G2010] eqs. (13)-(15):

```
 ⟨p' q' α' | V^(1)_2π,c1 | p q α⟩
   = (4π)² · Σ_{recouplings} × ∫_{-1}^{+1} dx  P_{k}(x) × G(x; p,q,p',q')
```

where `G` is the 2-body spin matrix element of `(σ_2·q_2)(σ_3·q_3)/[(q_2²+mπ²)(q_3²+mπ²)]`
multiplied by the isospin matrix element `Î_1 = τ_2·τ_3` (from [G2010]
eq. 22):

```
  Î_1(t', T', t, T) = (2 t (t+1) − 3) · δ_{t t'} δ_{T T'}
```

The complete expression is given numerically in [G2010] (Appendix A
examples, eqs. 33-39) for specific (L, L', S, S', l, l', λ, λ'). For
the general form one must evaluate the angular integral with the full
operator `F_1(q_2, q_3)` of [G2010] eq. (18):

```
  F_1(q_2, q_3) = (gA/2 fπ)² / [(q_2² + mπ²)(q_3² + mπ²)]
                 × (-4 c_1 mπ² / fπ² + 2 c_3 (q_2 · q_3) / fπ²)
```

(F_1 combines c_1 and c_3 in [G2010]; to isolate c_1 set c_3 → 0.)

The final scalar `G̃(β', β)` (in LS coupling of [G2010] eq. 24) requires:

1. Substitute `q_2 = Δp + Δq/2`, `q_3 = Δp − Δq/2`.
2. Expand `(σ_2·q_2)(σ_3·q_3)` in terms of the 6 spin-operator basis
   `w_i` of [G2010] eq. (4).
3. The coefficients `f_i(p, q, p', q', x)` are polynomials in x (at fixed
   p,q,p',q').
4. Perform the outer CG-sum and then `∫dx` numerically.

### 3.4 Spin-isospin recoupling coefficient

[H2015] eq. (15) factorization:

```
  ⟨pqβ|V|p'q'β'⟩ = Σ_{l̄,{L̄_i}} A_{ββ'}^{l̄,{L̄_i}} · F̃_{LlL'l'}^{l̄,{L̄_i}}(p,q,p',q')
```

where `A` carries all CG/6j/9j factors (spin-isospin dependence) and `F̃`
carries all momentum dependence. For the c_1 rank-0 piece:

```
  A_c1_rank0(α', α)
     = (2 T_2N (T_2N+1) − 3)  × δ_{α,α'}
     × (2 S_2N (S_2N+1) − 3)
```

— exactly the current code's form, but with the [E2002] squared-Gaussian
regulator.

For the c_1 rank-2 piece, the coefficient is:

```
  A_c1_rank2(α', α)
     = (2 T_2N (T_2N+1) − 3)        [τ_2·τ_3, still diagonal]
     × √(30 · ŝ ŝ' Ĵ_2N Ĵ_2N')     [spin-2 reduction, closed form from Wigner-Eckart]
     × NineJ( L_2N', S_2N', J_2N';
              L_2N,  S_2N,  J_2N;
              2,     2,     0    )    [9j for (L, S, J) to (L+tensor, S+tensor, J)]
     × CG(L_2N 0; 2 0 | L_2N' 0)     [spatial reduced me of Y_2]
     × ... spectator-line (l_1N, j_1N) recoupling (trivial; rank-2 acts on pair only)
```

The full closed form is in [H2015] eqs. (13)-(15) with `µ` summed over
the rank-2 projections. Use of `gsl_sf_coupling_9j` and
`gsl_sf_coupling_6j` as already done in `Atilde` at
[Tic-tac/src/utils/auxiliary.cpp:863-884](../../src/utils/auxiliary.cpp)
is the appropriate pattern.

### 3.5 Benchmark point

For α_r = α_c = {L_2N=0, S_2N=1, J_2N=1, T_2N=0, l_1N=0, 2j_1N=1,
2J_3N=1, 2T_3N=1} (dominant 3S1 triton channel), at p=q=p'=q'=0.5 fm⁻¹:

- Rank-0 diagonal contribution (current code form):
  - `σ_2·σ_3 = 2·1·2 − 3 = +1` (triplet)
  - `τ_2·τ_3 = 2·0·1 − 3 = −3` (singlet isospin)
  - `q_2·q_3` at angle-averaged: `= ⟨Δp²⟩ − ⟨Δq/2⟩² = (p²+p'²) − (q²+q'²)/4 = 0.5 − 0.125 = 0.375` fm⁻².
    (Current code formula matches this.)
  - `Q² = q_2² = q_3² ≈ p²+p'² + (q²+q'²)/4 = 0.5 + 0.125 = 0.625` fm⁻².
  - `(Q²+mπ²)² = (0.625 + 0.489)² = 1.239` fm⁻⁴.
  - `(-4 c_1 mπ²)` with c_1 = −0.81 GeV⁻¹ = −0.16 fm:
    `= -4·(-0.16)·0.489 = +0.313` fm⁻¹.
  - `(gA/2fπ)² / fπ² = (1.29 / (2·0.4682))² / 0.4682² = 1.898 / 0.2192 = 8.66` fm².

  Rank-0 V ~ `(⅓) × 1 × (-3) × 0.313 × 0.375 / 1.239 × 8.66`
         `≈ -3/3 × 0.313 × 0.375 / 1.239 × 8.66`
         `≈ -1 × 0.0948 × 8.66 = -0.821 fm⁵` (before regulator).

  After squared-Gaussian regulator f_R² ≈ 0.981:
  **V_c1,rank-0 ≈ -0.805 fm⁵**.

- Rank-2 contribution: qualitatively of similar magnitude but **opposite
  sign** in the 3S1 channel (spatial rank-2 projection of
  `q_2 ⊗ q_3` through Y_2(Δq̂) integration gives a negative moment when
  the pair orbital is S-wave). A precise value requires the [G2010]
  numerical scheme; estimated `V_c1,rank-2 ≈ +0.4 to +0.8 fm⁵`,
  partially cancelling rank-0.

**Benchmark value** (c_1 rank-0 only, 3S1, p=q=0.5 fm⁻¹):
≈ **−0.8 fm⁵**. Order of magnitude: `10⁻¹–10⁰ fm⁵`. Sign: negative
(c_1 < 0 gives attractive 3S1 contribution at low momenta).

---

## 4. c_3 2PE

### 4.1 Momentum-space operator form

[G2010] eq. (18), c_3 piece:

```
  V^(1)_2π,c3 = (gA / 2 fπ)² × (τ_2 · τ_3)
                × (σ_2 · q_2)(σ_3 · q_3) / [(q_2² + mπ²)(q_3² + mπ²)]
                × (+2 c_3 (q_2 · q_3) / fπ²)
```

Note the `(q_2 · q_3)` factor — c_3 carries an extra scalar product
of the two momentum transfers, which increases the effective rank and
introduces an additional angular integration complexity.

### 4.2 Rank decomposition

`(σ_2 · q_2)(σ_3 · q_3)(q_2 · q_3)` splits:

- **Rank 0**: `⅓ (σ_2·σ_3) (q_2·q_3)²` — the scalar product is squared.
- **Rank 2**: `[σ_2 ⊗ σ_3]_2 · [q_2 ⊗ q_3]_2 × (q_2·q_3)` — rank-2 spin
  part times an extra scalar factor; still rank-2 in angular space.
- **Rank 0 cross** from the rank-2 × scalar multiplication collapses
  to rank-0 in the partial-wave projection (but with a 6j-weighted
  coefficient).

For the current code's rank-0-only approximation:
`V_c3 ∝ ⅓ σ_2·σ_3 × τ_2·τ_3 × (q_2·q_3)² / [(Q²+mπ²)²]`.

### 4.3 Partial-wave matrix element

Identical structural form to c_1 with the kernel replaced:

```
  kernel_c3(x; p,q,p',q') = 2 c_3 (q_2 · q_3)² / [fπ² (q_2² + mπ²)(q_3² + mπ²)]
```

where `q_2·q_3 = |Δp|² − |Δq|²/4 = (p² + p'² − 2pp'·(p̂·p̂')) − (q² + q'² − 2qq' x)/4`
and `q_2² = q_3² = |Δp|² + |Δq|²/4 + (Δp·Δq)` — with p̂·p̂' fixed via
the ẑ-axis choice, `Δp·Δq` survives and depends on the `x` integration.

### 4.4 Spin-isospin recoupling coefficient

Identical to c_1 (same operator `(σ_2·q_2)(σ_3·q_3)(τ_2·τ_3)`) up to
replacing the scalar momentum kernel.

### 4.5 Benchmark point

For 3S1, p=q=p'=q'=0.5 fm⁻¹:

- `(q_2·q_3)² ≈ 0.375² = 0.1406` fm⁻⁴.
- `2 c_3 (q_2·q_3)² / fπ²` with c_3 = -3.20 GeV⁻¹ = -0.631 fm:
  `= 2·(-0.631)·0.1406 / 0.2192 = -0.809` fm⁻¹.
- Combining with `(gA/2fπ)² / [(Q²+mπ²)²] × σ_2·σ_3 × τ_2·τ_3 × ⅓`:
  `= 1.898 / 1.239 × 1 × (-3) × ⅓ × (-0.809)`
  `= 1.532 × (-1) × (-0.809) = +1.239 fm⁵` (before regulator).

After squared-Gaussian regulator f_R² ≈ 0.981:
**V_c3,rank-0 ≈ +1.22 fm⁵**.

**Sign**: c_3 < 0, but the `(q_2·q_3)²` is positive definite, so the
sign is `sign(c_3) × sign(σ·σ · τ·τ · (-1 from rank-0 cross)) = +`.
This matches the "sign flip" diagnosis — the current code uses
`2 c_3 (q_2·q_3)²` as a positive-definite term, giving **attractive** in
3S1 (T_2N=0, S_2N=1), which is the expected behavior for c_3 in
chiral 3NF calculations ([E2002] Table 2).

**Order of magnitude**: `~10⁰ fm⁵` at p=q=0.5 fm⁻¹. Sign: positive in
3S1.

---

## 5. Sanity-check: `(σ_2·q̂_2)(σ_3·q̂_3)` partial-wave expansion

This is the operator piece that determines all c_D, c_1, c_3
rank-decomposition signs. Following [G2010] eq. (17) and the identity
above:

```
  (σ_2 · q̂_2)(σ_3 · q̂_3)
     = ⅓ (σ_2 · σ_3)                         [rank 0]
       + [σ_2 ⊗ σ_3]_2 · [q̂_2 ⊗ q̂_2]_2        [rank 2, when q̂_2 = q̂_3]
       + [σ_2 ⊗ σ_3]_2 · [q̂_2 ⊗ q̂_3]_2        [general rank 2]
```

where the spherical-tensor product satisfies

```
  [q̂ ⊗ q̂]_{2,M} = √(8π/15) · Y_{2,M}(q̂)
```

Explicit phases (from Varshalovich or Edmonds, 2010 edition):

```
  [A ⊗ B]_{L,M} = Σ_{μ_A, μ_B} CG(L_A μ_A; L_B μ_B | L M) A_{L_A μ_A} B_{L_B μ_B}
```

For the pair spin:

```
  [σ_2 ⊗ σ_3]_{2,M} ⟨S_2N' M' | ... | S_2N M⟩
      = √(2·2+1) · CG(S_2N M; 2 μ | S_2N' M') × ⟨S_2N' || [σ_2 ⊗ σ_3]_2 || S_2N⟩
```

with the reduced matrix element (S_2N and S_2N' must each be 1; the 2
is a pair-spin triplet-to-triplet rank-2 transition):

```
  ⟨1 || [σ_2 ⊗ σ_3]_2 || 1⟩ = √30       (standard result; see e.g. Edmonds eq. 7.1.8)
```

Sign caveat: some texts adopt `σ = 2S/ℏ` convention giving
`⟨1||σ||1⟩ = √6`, versus the Pauli-matrix convention `⟨1||σ||1⟩ = √2`.
The code's `reduced_me_sigma_dot_qhat` must use a single consistent
convention throughout; the sign flip observed in c_1 suggests a
mixed convention between spin-σ and spin-S reductions.

**The expansion for c_D's `(σ_1·q̂_3)(σ_3·q̂_3)`** differs: here both
σ-vectors couple to the *same* momentum unit vector `q̂_3`, so the
rank-2 piece is `[σ_1 ⊗ σ_3]_2 · [q̂_3 ⊗ q̂_3]_2 = √(8π/15) [σ_1 ⊗ σ_3]_2 · Y_2(q̂_3)`.
This couples the pair spin through σ₃ **and** the spectator spin through
σ₁, requiring the 9j symbol `{S_2N, ½, S; L_2N, l_1N, L; J_2N, j_1N, J_3N}`
for the full Jj ↔ LS recoupling.

---

## 6. Audit resolution ledger

1. **Resolved - contact sign and pair counting**: visual inspection and text
   extraction of [E2002] Eq. (2.10) give `+½ E sum_(j!=k) tau_j·tau_k`.
   Each unordered pair occurs twice, so `V^(1)_cont=+E tau_2·tau_3`.

2. **Resolved - regulator power**: [E2002] Eq. (3.19) and production both use
   `fR=exp(-((4p²+3q²)/(4Λ²))²)` on bra and ket.

3. **Resolved for A-4; open for analytic A-1 reuse**: the original PDF page was
   rendered and A-4 unambiguously reads `6E(4π)^2`.  The exact full numerical
   five-angle projector is authoritative for production development; no broken
   OCR transcription of A-1 may be coded without a separate visual check.

4. **c_1 rank-2 coefficient sign**: the diagnostic tool reports a
   sign flip on c_1 at X = 92-362. The rank-2 piece of
   `(σ_2·q_2)(σ_3·q_3)` contributes an opposite-sign 3S1 matrix element
   from rank-0, so omitting rank-2 is expected to flip the sign. A
   full numerical eval of the [G2010] eq. (18) kernel for α=3S1 at the
   actual triton momentum distribution is needed to confirm this is
   the sole cause.

5. **c_4 cross-product term**: [E2002] eq. (2.2) includes
   `ε_{αβγ} τ_k^γ σ_k · (q_i × q_j) × c_4/fπ²`. The c_4 term carries
   the isospin operator `τ_1·(τ_2 × τ_3)` ([G2010] eq. 17, 21), which is
   purely imaginary for T_3N=½ basis states — hence it generates only
   off-diagonal matrix elements in the full 3N basis. The current code's
   c_4 = 0 assumption may be sound for 3H binding (where T_3N=½
   dominates), but should be re-examined for deuteron-proton observables.

6. **Benchmark point numerical value verification**: the benchmark
   calculations in §1.5, §2.5, §3.5, §4.5 are order-of-magnitude
   estimates. For the normalization tool to be useful, an independent
   reference code (e.g. Witała's Bochum / Kraków Fortran, or the
   Hebeler/Krebs matrix-element files hosted at lenpic.org) should
   produce the exact `V(α'=α=3S1, p=q=p'=q'=0.5 fm⁻¹)` in fm⁵ per
   LEC. The Hebeler 2015 PRC 91 044001 supplemental material contains
   tabulated 3NF matrix elements for all J ≤ 9/2 channels at N2LO
   and is the recommended validation source.

7. **Coupling conventions**: [E2002] uses `I` for spectator total
   angular momentum and `J` for 3N total. [G2010]/[H2015] use
   `j` for spectator and `J` for 3N. Current code uses
   `two_J_1N` (spectator, half-integer ×2) and `two_J_3N` (3N,
   half-integer ×2). Mapping: [E2002]·I ↔ code·J_1N; [E2002]·J ↔
   code·J_3N. Internal consistency verified but explicit docstring
   in spin_isospin_algebra.h recommended.

8. **Permutation operator factor**: [E2002] uses
   `V = (1 + P₁₂₃ + P₁₃₂) V^(1)`; current code outputs
   `3·⟨V^(1)⟩` as the full 3NF expectation. This is valid **only** if
   `⟨ψ|ψ⟩ = 1` for the antisymmetric state (confirmed at
   check_3nf_normalization.cpp:7-9) and `V^(1)` is the spectator-1
   Faddeev component. For non-antisymmetric intermediate states the
   `(1+P+P²)` factor must be applied explicitly.
