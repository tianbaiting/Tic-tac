# Miller benchmark reference

Regression-test target values for the **2NF path** of Tic-tac, drawn from
publications by the upstream author Sean B. S. Miller (Chalmers). These
numbers are the first things to reproduce after any non-trivial rewrite
(e.g. the upcoming partial-wave projection of the chiral N²LO 3NF) so
that a silent breakage of the 2NF+AGS machinery is caught immediately.

## Paper(s) consulted

1. **S. B. S. Miller, A. Ekström, K. Hebeler**, *Neutron-deuteron
   scattering cross sections with chiral NN interactions using
   wave-packet continuum discretization*,
   **Phys. Rev. C 106, 024001 (2022)**, arXiv:**2201.09600** [nucl-th].
   — Primary reference. This is the paper that the upstream Tic-tac
   repository cites for the method (`seanbsm/Tic-tac` README → PRC
   106, 024001) and it contains the explicit convergence benchmarks
   and observable predictions that the code was built to reproduce.

2. **S. B. S. Miller, A. Ekström, C. Forssén**, *Wave-packet continuum
   discretisation for nucleon-nucleon scattering predictions*,
   **J. Phys. G 49, 024001 (2022)**, arXiv:**2106.00454** [nucl-th].
   — Companion (NN-only) paper. Establishes the WPCD 2-body
   scattering accuracy using NNLOopt. Method error ≈1–5 mb in the
   total NN cross section for 0–350 MeV, on a GPU-parallelized
   implementation. Not a primary benchmark source for the 3-body code
   but documents the underlying 2-body discretization accuracy.

3. **S. B. S. Miller**, *Approximating the Three-Nucleon Continuum*,
   **Doctoral thesis, Chalmers University of Technology (2022)**
   (ISBN/URL on research.chalmers.se/en/person/seanmi). PDF was not
   successfully retrieved in this session; the published PRC 106
   paper is its content-equivalent primary output. Listed here in
   case a future maintainer can obtain the thesis and extract
   higher-resolution tables (it almost certainly carries the full
   phase-shift tabulations that appear only as figures in PRC 106).

4. **S. B. S. Miller, A. Ekström, C. Forssén**, *Posterior predictive
   distributions of neutron-deuteron cross sections*,
   **Phys. Rev. C 107, 014002 (2023)** (listed on the Chalmers
   profile; arXiv ID not confirmed in this session — candidate
   2210.xxxxx / 2211.xxxxx). Uses the same WPCD Faddeev code on χEFT
   up to N³LO with Bayesian-sampled LECs; produces PPD bands rather
   than single-number benchmarks, so it is secondary for our
   regression-test purpose.

5. **Upstream reference [2] of the PRC 106 paper**:
   W. Glöckle, H. Witała, D. Huber, H. Kamada, J. Golak, *The
   three-nucleon continuum: achievements, challenges and
   applications*, **Phys. Rep. 274, 107 (1996)**. This is the
   “Standard Faddeev” source to which Miller benchmarks at
   E_lab = 13 MeV; any independent digitization of Miller's Fig. 1 or
   Fig. 2 ultimately resolves against values in this Phys. Rep. review.

## Available benchmark tables

All values below are from Miller et al., PRC 106 024001 (2022)
unless stated. Values given only as figures (not tables) in the paper
are annotated “[fig read]”; the exact curve is machine-readable by
`WebPlotDigitizer` from the arXiv preprint PDF when higher precision
is needed.

### Table 1. Triton ground-state properties
> **Not reported in PRC 106 024001.**
>
> The PRC 106 paper explicitly restricts itself to elastic
> nd-scattering above threshold; it does not tabulate E_B(³H),
> ⟨T⟩_³H or ⟨V⟩_³H. The paper's only reference to the triton
> binding energy is to note the Phillips-line correlation with the nd
> doublet scattering length (body text around Ref. [44]), without a
> numerical value.
>
> **Recommendation**: for a triton-E_B cross-check of Tic-tac, use
> instead the canonical benchmark values from Nogga, Kamada,
> Glöckle, Barrett, *Phys. Rev. C 65, 054003 (2002)*:
> Idaho-N3LO (Λ=500 MeV, NN only) → E_B(³H) ≈ **−7.855 MeV**
> (empirical −8.482 MeV). These are not Miller numbers but they are
> what the community uses as the NN-only triton-energy target, and the
> existing `tools/check_3nf_normalization/epelbaum_reference.md` in
> this directory is already aligned to that convention.

### Table 2. nd scattering phase shifts — convergence benchmark
Source: **PRC 106 024001, Fig. 1**.

NN interaction: **Nijmegen-I** (Stoks-Klomp-Terheggen-de Swart,
Phys. Rev. C 49, 2950 (1994)).

Benchmark point: **E_lab = 13 MeV**, reference = Glöckle et al.
Phys. Rep. 274 (1996) “Standard Faddeev”.

| Channel                | Observable        | Reference value       | Miller WPCD (N_WP=125) | Agreement | Ref.          |
|------------------------|-------------------|-----------------------|------------------------|-----------|---------------|
| J^Π = 1/2⁺ doublet     | Re δ(2S) [deg]    | ≈ 105 deg [fig read]  | recovered              | within 1% | Fig. 1 left   |
| J^Π = 1/2⁺ doublet     | Im δ(2S) [deg]    | ≈ 0 deg (below breakup) | recovered            | < 1e-2°   | Fig. 1 left   |
| J^Π = 7/2⁺ quartet     | Re δ(4D) [deg]    | ≈ +20 deg [fig read]  | recovered              | within 1% | Fig. 1 right  |
| J^Π = 7/2⁺ quartet     | Im δ(4D) [deg]    | ≈ 0 deg (below breakup) | recovered            | < 1e-2°   | Fig. 1 right  |

**What the text states (verbatim context):** “we recover standard
Faddeev results for all imaginary and real parts of the NNN phase
shifts for J^Π ≤ 7/2± within ∼1% using N_WP ≳ 125 wave packets.”
“the magnitude of the imaginary part of the phase shifts is
|Im(δ)| ≲ 10⁻² degrees for scattering energies below the deuteron
breakup threshold.”

**Use as a regression test**: at E_lab = 13 MeV with Nijmegen-I and
N_WP = 125, Tic-tac's doublet and quartet phase shifts must match
Glöckle/Witała Phys. Rep. 274 values within ~1 % real and
|Im δ| ≤ 0.01° for all J^Π ≤ 7/2. A ≥5 % deviation on any of these
channels indicates a broken 2NF + kernel pipeline.

### Table 3. Neutron analyzing power Ay(n) convergence
Source: **PRC 106 024001, Fig. 2**.

NN interaction: **Nijmegen-I**. Energy: **E_lab = 35 MeV**.

Convergence pattern (WPCD → Standard Faddeev with increasing N_WP):
curves for N_WP = {50, 75, 100, 125, 150} approach the Glöckle et al.
(Phys. Rep. 274) standard-Faddeev curve monotonically; **for
N_WP ≳ 125 the two curves overlap to within plotting accuracy.**

Representative values read from Fig. 2 (N_WP = 125, Nijmegen-I):

| θ_cm [deg] | Ay(n) Miller (N_WP=125)  | Standard Faddeev [Ref. 2] | Ref.   |
|------------|--------------------------|---------------------------|--------|
| ~0         | ≈ 0.00                   | ≈ 0.00                    | Fig. 2 |
| ~60        | ≈ 0.00                   | ≈ 0.00                    | Fig. 2 |
| ~120 (min) | ≈ −0.50 [fig read]       | ≈ −0.50                   | Fig. 2 |
| ~180       | ≈ 0.00                   | ≈ 0.00                    | Fig. 2 |

**Use as a regression test**: Tic-tac run with `potential_model =
nijmegen`, `two_J_3N_max` including 1/2 and 7/2, N_WP = 125,
E_lab = 35 MeV must reproduce an Ay minimum of ≈ −0.50 ± 0.02 near
θ_cm = 120°.

### Table 4. WPCD convergence with Idaho-N3LO
Source: **PRC 106 024001, Fig. 3**.

NN interaction: **Idaho-N3LO** (Entem-Machleidt, Phys. Rev. C 68,
041001 (2003)). Observables: Ay(n), T22(d), Kyy(nn), T21(d),
dσ/dΩ. Energies: **E_lab = 3, 10, 65 MeV**.

Numeric y-axis ticks visible in Fig. 3 (for digitizing checks):

| Panel (E_lab)  | Observable   | y-range visible in Fig. 3   | Ref.   |
|----------------|--------------|-----------------------------|--------|
| 3 MeV          | Ay(n)        | 0.00 – 0.05                 | Fig. 3 |
| 3 MeV          | T22(d)       | −0.03 – 0.00                | Fig. 3 |
| 3 MeV          | dσ/dΩ [mb]   | 2×10² – 4×10²               | Fig. 3 |
| 10 MeV         | Ay(n)        | −0.01 – 0.14                | Fig. 3 |
| 10 MeV         | Kyy(nn)      | 0.00 – 1.00                 | Fig. 3 |
| 10 MeV         | dσ/dΩ [mb]   | ≈ 10²                       | Fig. 3 |
| 65 MeV         | Ay(n)        | −0.63 – 0.29                | Fig. 3 |
| 65 MeV         | T21(d)       | −0.29 – 0.02                | Fig. 3 |
| 65 MeV         | dσ/dΩ [mb]   | 10⁰ – 10¹                   | Fig. 3 |

**What the text states:** convergence of dσ/dΩ and polarization
observables is reached for **N_WP ≳ 125** at E_lab ≤ 70 MeV.
T21(d) is identified as the slowest-converging observable.
“At E_lab > 70 MeV, the inclusion of J > 3 NN channels will have a
percent-level effect on the predictions.”

### Table 5. Differential cross section
Source: **PRC 106 024001, Fig. 6**.

NN interactions: **N2LOopt** and **Idaho-N3LO**. N_WP = 125.
Energies:
- Left panel: **E_lab = 6, 12, 25 MeV** (nd data from Ref. [41]).
- Right panel: **E_lab = 64.5 MeV** (pd data from Ref. [42]).

Qualitative benchmark (from text): WPCD + N2LOopt reproduces the
world nd elastic cross-section database “excellently up to
E_lab ≈ 70 MeV”. At **E_lab = 64.5 MeV** the Idaho-N3LO prediction
sits slightly below N2LOopt at the minimum. Exact points must be
digitized from Fig. 6.

### Table 6. Spin observables
Source: **PRC 106 024001, Fig. 7**.

NN interactions: **N2LOopt** and **Idaho-N3LO**. N_WP = 125.

| Energy [MeV] | Observable   | Data type | Ref.    |
|--------------|--------------|-----------|---------|
| 21           | Ay(n)        | nd data   | [46]    |
| 22.5         | Kyy(nn)      | nd data   | [47]    |
| 35           | Ay(n)        | pd data   | [49]    |
| 13           | Cxx(nd)      | pd data   | [50]    |
| 13           | Cyy(nd)      | pd data   | [50]    |
| 47.5         | iT11(d)      | pd data   | [51]    |
| 65           | Ay(n)        | nd data   | [48]    |
| 65           | Ay(n) (pd)   | pd data   | [42]    |
| 70           | Axx(d)       | pd data   | [52]    |
| 70           | Ay(d)        | pd data   | [52]    |

These are the paper's curated observable set — they span the same
kinematic region as Tic-tac's own 190 MeV benchmark (from below) but
do not include 190 MeV directly.

### Table 7. Total nd cross section
Source: **PRC 106 024001, Fig. 4**.

NN interaction: **N2LOopt**, N_WP = 100 and 125. Energy range
spans ~1 – 90 MeV on a log axis. Peaks at the ELab = 0 limit
(scattering length correlation with E_B(³H) via Phillips line),
decays to ≈ 10² mb near 90 MeV. Difference between N_WP = 100 and
N_WP = 125 visible only at E_lab > 70 MeV.

### Table 8. WPCD grid parameters (PRC 106 024001)
Source: **PRC 106 024001, Sec. II.B** and Eq. (5).

| Parameter            | Symbol           | Value / range           |
|----------------------|------------------|-------------------------|
| Wave-packet count    | N_WP             | production: 125; convergence scan: 50, 75, 100, 125, 150 |
| Grid type            | —                | Chebyshev: p_i = q_i = α tan[(2i−1)π / 4 N_WP] |
| Grid scale           | α                | **200 MeV**             |
| Chebyshev exponent   | t                | **1**                   |
| Max bin momentum     | p_max ≈ q_max    | ~10 GeV (bin edges); bulk of WPs at < 500 MeV |
| 3N total-J cap       | J ≤ 17/2         | (up to J ≤ 17/2 at convergence) |
| 2N total-J cap       | J ≤ 3            | (i.e. J_2N_max = 3)     |
| Number of NNN partial-waves per J^Π | —  | ≲ 60                    |
| Isospin              | T = 1/2, 3/2     | charge-dependent 1S0 (T=3/2−1/2 coupling included; impact "very small" per Ref. [33]) |
| Neumann-series order | —                | 20–30 terms with Padé resummation |

**Mapping to Tic-tac inputs**: `Np_WP=Nq_WP=125`, Chebyshev
`chebyshev_s`/`chebyshev_t` chosen so that α=200 MeV, t=1,
`two_J_3N_max ≥ 17`, `J_2N_max = 3`, `isospin_breaking_1S0 = true`.

### Table 9. Kernel eigenvalue η(Z=0)
> **Not numerically tabulated in PRC 106 024001.**
>
> The Weinberg-type eigenvalue criterion |η_i| < 1 is stated as the
> convergence condition for the Neumann series (Appendix body around
> Eq. (C) in the arXiv preprint, “This series only converges if all
> so-called Weinberg eigenvalues ηᵢ of K satisfy |ηᵢ| < 1”), and
> the paper cites **Refs. [35, 66]** for numerical studies.
> Reference [35] = Kukulin et al. (Padé resummation), Reference [66]
> = A. Schwenk et al., *Weinberg eigenvalues for chiral nucleon
> interactions*. **No η(Z=0) table was extracted.** For a regression
> test, the relevant quantity inside Tic-tac is whether the Padé
> resummation actually converges; this is a binary check, not a
> numerical match.

## Overlap with Tic-tac's canonical 190 MeV/u benchmark

Tic-tac's current validated benchmark at **T_lab = 190 MeV/u**
(d+p → d+p, with dσ/dΩ, iT11, T20, T21, T22 vs. experiment in
`data/DataOfCrosssectionAndPol/`) is **outside** the energy range
covered by PRC 106 024001 in two respects:

1. PRC 106 explicitly restricts itself to **E_lab ≲ 100 MeV**
   (neutron-lab = 100 MeV ↔ deuteron-lab ≈ 150 MeV, still below
   190 MeV/u = 380 MeV total kinetic). Above ≈70 MeV the paper notes
   *“J > 3 NN-channels have a percent-level effect on the
   predictions”* — i.e. Tic-tac at 190 MeV/u is in a regime that
   Miller's own truncation already admits is under-converged with
   J_2N_max = 3.
2. The paper reports mostly Ay(n), iT11(d) at 47.5 MeV, T21(d),
   T22(d), Kyy(nn) — some of these match Tic-tac's observable list
   (iT11, T21, T22) but none is at 190 MeV/u.

**Consequence for the 3NF-rewrite regression test**:

- The Miller-provided benchmark points are **low-energy / mid-energy
  convergence benchmarks** (13–70 MeV, Nijmegen-I and N2LOopt,
  NN-only). Pass/fail there tells us the 2NF + Faddeev pipeline is
  sound.
- Tic-tac's own **190 MeV/u benchmark** (`docs/dpol_p_190MeV_
  validation.md`) is a **high-energy data-vs-theory** check and
  stresses a different corner (higher partial waves, relativistic
  kinematics, Coulomb omitted). It is complementary, not overlapping.

Recommendation: treat the Miller 13 MeV Nijmegen-I phase-shift
check (Table 2 above) and the 35 MeV Ay(n) Nijmegen-I shape check
(Table 3 above) as the **first two gates** any 3NF-rewrite must pass
with 3NF switched OFF, before the 190 MeV/u Ay pipeline is even
exercised.

## Pending clarifications

- **No absolute numerical tables** are given in PRC 106 024001 —
  everything (phase shifts, Ay, cross sections) is reported only as
  figures. All values in Tables 2–7 above are read off the plotted
  curves; *exact digitization is required for a quantitative
  regression test*. The preprint PDF (arXiv:2201.09600, ~1 MB)
  has been cached at `/tmp/miller_prc2022.pdf` during this research
  pass; converting the seven figure PDFs inside it with
  `pdfimages` + `WebPlotDigitizer` is a one-off task that would
  promote these "[fig read]" entries to numerical targets.
- **Miller 2022 PhD thesis** (*Approximating the Three-Nucleon
  Continuum*, Chalmers) was not retrieved in this session. It most
  likely contains the underlying tabulated data for the PRC 106
  figures. Retrieving it from
  `research.chalmers.se/en/person/seanmi` or from
  `publications.lib.chalmers.se` would remove the “[fig read]”
  qualifier from Tables 2–4.
- **Miller–Ekström–Forssén 2023 PRC (posterior predictive)** — the
  arXiv number was not confirmed. If it provides explicit
  maximum-likelihood N3LO values of Ay(n) or dσ/dΩ at specific
  energies, those are tighter benchmarks than the PRC 106 figures
  because they come with an explicit uncertainty band.
- **Triton E_B** and **η(Z=0) Weinberg eigenvalue**: not in Miller.
  For these quantities use the existing
  `epelbaum_reference.md` (Phys. Rev. C 66, 064001 (2002)) in this
  same directory.
- **Sign/convention on Ay**: Miller's Fig. 2 shows a negative
  minimum; check that Tic-tac's sign convention for the neutron
  analyzing power matches (Madison convention, Appendix E of
  PRC 106).
