# Authoritative N2LO three-nucleon-force conventions

**Status:** locked for the exact-oracle and production-PWD work after
`7ab45a0` (2026-08-13).  A row marked *open* may not be filled by numerical
tuning.  This document and `three_nf_equation_contract.md` supersede older
normalization statements in historical audit notes.

## Convention table

| Quantity | Tic-tac convention | Evidence |
|---|---|---|
| Spectator component | Particle 1 is the spectator label; particles 2 and 3 form the antisymmetric pair. `W1` denotes the normalized PWD of `V^(1)`, and `V_3N=V^(1)+V^(2)+V^(3)`. | Epelbaum 2002 Eqs. (2.2), (2.10), (3.14); production interface and pair-exchange oracle. |
| Jacobi momenta | In the 3N center-of-mass frame, `k1=q`, `k2=p-q/2`, `k3=-p-q/2`; hence `p=(k2-k3)/2` and `q=k1`. | Golak 2010 Eq. (26); explicit conservation test. |
| Transfers | `Qi=ki'-ki`, so `Q1=q'-q`, `Q2=p'-q'/2-(p-q/2)`, `Q3=-p'-q'/2-(-p-q/2)` and `Q1+Q2+Q3=0`. | Golak 2010 Eqs. (17), (26); `test_jacobi_transfers_conserve_momentum`. |
| Plane-wave radial state | `<p'|p>=delta(p'-p)/(p'p)`, and analogously for `q`. The 3N state is the tensor product of the two radial states and a normalized spin-angular channel. | Miller 2022 Eq. (3) and App. A Eq. (A-4); WP construction and cache measure. |
| Spherical harmonics and CG phases | Unit-normalized `Y_lm` with the Condon-Shortley phase; CG coefficients use the same convention as SymPy/GSL. | Miller 2022 App. A; independent CG/6j/9j tests. |
| Pair antisymmetry | Allowed channels obey `(-1)^(l_pair+s_pair+t_pair)=-1`. | Channel builder and independent LS-oracle constructor. |
| Spin/isospin operators | `sigma` and `tau` are Pauli matrices, not spin/isospin generators divided by two. Thus `tau2.tau3=(-3,+1)` for pair isospin `(0,1)`. | Explicit 8-state Pauli enumeration and Epelbaum Eq. (A-4). |
| Raw operator units | With momenta and masses in fm^-1, `c1,c3,c4` in fm, and `cD,cE` dimensionless, the unnormalized momentum-space operator and raw PWD have units fm^5. | Dimensional audit; Golak Table 2 reports fm^5. |
| Fourier normalization | Apply one `(2*pi)^-3` factor for each independent Jacobi relative coordinate: `N_FT=(2*pi)^-6`. | Tic-tac 2NF convention applied to `p` and `q`; WP cache has no later `(2*pi)` factor. |
| Contact angular normalization | A constant S-wave projection gives `(4*pi)^2`; consequently the normalized `cE` PWD coefficient is `(4*pi)^2(2*pi)^-6=1/(4*pi^4)`. | Epelbaum 2002 Eq. (A-4), rendered PDF; independent five-angle projector. |
| Regulator | `f_R(p,q)=exp(-[(p^2+3q^2/4)/Lambda^2]^2)` on both bra and ket. | Epelbaum 2002 Eq. (3.19). |
| Constants in production | `hbarc=197.327 MeV fm`, `gA=1.289`, `fpi=92.2 MeV`, `mpi=138.039 MeV`, `Lambda_chi=700 MeV`. Published-number tests may override these explicitly. | `include/constants.h`; constructor. |
| LEC units/signs | `c1,c3,c4[GeV^-1] -> c_i*hbarc/1000 [fm]`; `D=cD/(fpi^2 Lambda_chi)` and `E=cE/(fpi^4 Lambda_chi)`. The spectator contact is `+E tau2.tau3`. | Epelbaum 2002 Eqs. (2.10), (2.12). |
| LEC-set compatibility | A published `(cD,cE)` pair is admissible only with its stated regulator, cutoff, pion/decay constants, and `+E tau.tau` sign convention. The final benchmark set is *open* and must be documented before use. | Scientific acceptance requirement; no set selected yet. |
| AGS insertion | Tic-tac's elastic kernel is `K=P v+(1+P)W1`; the permutation is on the left of the spectator component. | Deltuva 2009 Eq. (7a), reduced with `tG0=vG` and `G0(1+tG0)=G`; noncommuting matrix discriminator. |
| WP cell measure | Integrate `p' q' p q dp' dq' dp dq`, divide by the square root of the four bin widths, and convert fm^5 to MeV with `hbarc^-5`. | `W1_PW_cache::build`; cache/direct operator test. |

## Spectator-1 N2LO operator before PWD

With `D_i=Qi^2+mpi^2`, the raw full-vector component is

```text
V_2pi^(1) = gA^2/(4 fpi^4) (sigma2.Q2)(sigma3.Q3)/(D_2 D_3)
             * {[-4 c1 mpi^2 + 2 c3 (Q2.Q3)] (tau2.tau3)
                + c4 [tau1.(tau2 x tau3)] [sigma1.(Q2 x Q3)]},

V_D^(1)   = -gA cD/(8 fpi^4 Lambda_chi) 1/D_1
             * {(tau1.tau2)(sigma1.Q1)(sigma2.Q1)
                +(tau1.tau3)(sigma1.Q1)(sigma3.Q1)},

V_E^(1)   = +cE/(fpi^4 Lambda_chi) (tau2.tau3).
```

Multiply every term by `f_R(p',q')f_R(p,q)`.  These expressions are raw:
the angular projector acts first, followed by `(2*pi)^-6` for Tic-tac.  No
empirical scale, post-hoc Hermitian averaging, or `Ay`-dependent sign choice is
allowed.

Two phase discriminators are locked:

- `<t23'=0|tau1.(tau2 x tau3)|t23=1> = -2 sqrt(3) i` for total isospin
  `T=1/2`, `M_T=+1/2` in the stated Condon-Shortley basis;
- the reverse element is its complex conjugate.  Combined with the spin and
  momentum cross product, the complete `c4` matrix element is Hermitian.

## Evidence boundaries

- **Verified:** all five full-vector terms, momentum conservation, pair
  symmetry, Hermiticity, the `c4` phase, Golak Table 2 `G(1,1)` and `G(2,1)`,
  the generic five-angle projector, and the absolute `cE` normalization.
- **Not yet verified:** the LS-to-Jj production recoupling for general channels,
  exact production PWD values for `c1,c3,c4,cD`, and their WP convergence.
- **Forbidden inference:** passing legacy rank-zero tests does not validate the
  omitted tensor, off-diagonal, or higher-partial-wave sectors.
