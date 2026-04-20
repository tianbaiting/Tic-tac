# Per-LEC ⟨W⟩_3H reference values

**Primary source**: Epelbaum, Nogga, Glöckle, Kamada, Meißner, Witała,
Phys. Rev. C 66, 064001 (2002). arXiv:nucl-th/0208023.
Table 2 (³H column, Λ=500 MeV NNLO row).

**Backup / related**: Witała, Golak, Skibiński, Kamada, Nogga, Epelbaum,
Phys. Rev. C 77, 034004 (2008) — uses a distinct LEC convention
(c_D/c_E defined via eq. (3.3) of Navrátil, Few-Body Syst. 41, 117 (2007),
with the 4πc_D/F_π²Λ_χ redefinition). The LEC values used in our code
(c_D=−0.20, c_E=−0.205) follow the modern Navrátil/Witała convention;
Epelbaum's 2002 values (c_D=+3.6, c_E=+0.37 at Λ=500) are in the older
convention. Values are linearly rescaled below. The rescaling is valid
because W^(1) is strictly linear in each LEC separately.

## Our code LECs and Epelbaum reference LECs

| LEC  | Our code value (Witała / Navrátil conv.) | Epelbaum 2002 value (original conv.) |
|------|------------------------------------------|--------------------------------------|
| c_D  | −0.20                                    | +3.6      (eq. 3.20)                 |
| c_E  | −0.205                                   | +0.37     (eq. 3.20)                 |
| c_1  | −0.81 GeV⁻¹                              | −0.81 GeV⁻¹ (inherited from 2NF [10]) |
| c_3  | −3.20 GeV⁻¹                              | ≈ −3.40 GeV⁻¹ (Epelbaum EGM500 2NF)  |
| Λ_3NF| 500 MeV                                  | 500 MeV                              |

Sign/magnitude difference in c_D, c_E between 2002 Epelbaum and
Navrátil/Witała conventions arises from the D, E definitions:
Epelbaum (eq. 2.11–2.12): D = c_D / (f_π² Λ_χ), E = c_E / (f_π⁴ Λ_χ),
Λ_χ = 700 MeV. Navrátil: D = c_D / (f_π² Λ_χ) with 4πc_D and flipped sign
convention for the contact vertex. The overall W^(1) is linear in the
physical coupling, so a proportional rescale of the LEC value rescales
⟨W⟩ by the same factor.

## Epelbaum's raw numbers (Table 2, ³H column, NNLO Λ=500 MeV)

| Epelbaum channel | ⟨W⟩_ref at Epelbaum LECs [MeV] |
|------------------|-------------------------------|
| c-terms (c_1+c_3 combined) | −0.39    |
| D-term (c_D only)          | +0.81    |
| E-term (c_E only)          | −0.74    |
| all (full 3NF)             | −0.32    |

## Rescaled to our code's LEC values

Linear rescaling: ⟨W⟩_ours = ⟨W⟩_Epelbaum × (LEC_ours / LEC_Epelbaum).

| Channel       | LEC value  | Rescaling factor            | ⟨W⟩_ref [MeV] | Source / method                |
|---------------|------------|-----------------------------|---------------|-------------------------------|
| c_E_only      | c_E=−0.205 | (−0.205 / 0.37) = −0.5541   | +0.410        | Epelbaum Table 2 × rescale    |
| c_D_only      | c_D=−0.20  | (−0.20 / 3.6)   = −0.05556  | −0.045        | Epelbaum Table 2 × rescale    |
| c_1_only      | c_1=−0.81  | c-terms split ~20%          | −0.073        | ESTIMATE — Epelbaum paper gives only combined c-terms (see note) |
| c_3_only      | c_3=−3.20  | c-terms split ~80%, × 0.94  | −0.294        | ESTIMATE — Epelbaum paper gives only combined c-terms (see note) |
| c_1+c_3 (comb)| both       | × (3.20/3.40) ≈ 0.94        | −0.367        | Epelbaum Table 2 × (our c_3 / Epelbaum c_3) — approximation      |
| full_Witala   | all four   | sum of four ≈ −0.002        | ≈ 0.00        | Derived sum — note near-zero from cancellations at our LECs; NOT equal to Epelbaum's −0.32 MeV, which was at Epelbaum's own fitted c_D, c_E |

## Notes on the c_1 vs c_3 split (ESTIMATE)

Epelbaum's Table 2 reports only the **combined c-terms** (−0.39 MeV at
Λ=500 MeV, their c_1=−0.81, c_3≈−3.40 GeV⁻¹). A per-LEC decomposition
is not directly published. We estimate the split using the typical
literature observation that the c_3 two-pion-exchange term dominates
over the c_1 M_π² term by roughly 4:1. Hence:

- c_3 share ≈ 80% of combined c-terms
- c_1 share ≈ 20% of combined c-terms

These are NOT direct paper values. The more robust comparison below is
the **combined c_1+c_3 ratio**, computed from the sum of our c_1_only
and c_3_only rows against Epelbaum's −0.39 × (our c_3 / Epelbaum c_3)
≈ −0.367 MeV.

## Sign convention

In Epelbaum's table ⟨W⟩ is the ³H expectation value of the 3NF in the
Hamiltonian H = T + V_NN + V_3NF. A negative ⟨W⟩ means the 3NF is net
attractive (increases binding). The total ΔE_3NF ≈ −0.32 MeV at
Epelbaum's LECs is consistent with this.

Our code's reported "3·⟨W⟩" should directly correspond to Epelbaum's
⟨W⟩_ref if the normalization is correct. A code that is too strong
by a factor X would report 3·⟨W⟩_code = X · ⟨W⟩_ref.

## Final reference table (feed to check_3nf_normalization.cpp)

| Channel       | ⟨W⟩_ref [MeV]   | Type                    |
|---------------|-----------------|--------------------------|
| c_E_only      | +0.410          | transcribed + rescaled  |
| c_D_only      | −0.045          | transcribed + rescaled  |
| c_1_only      | −0.073          | ESTIMATE (c_1/c_3 split)|
| c_3_only      | −0.294          | ESTIMATE (c_1/c_3 split)|
| full_Witala   | −0.002          | sum of four (dominated by near-cancellation at our LECs) |

**Robust aggregate**: c_1+c_3 combined = −0.367 MeV (transcribed + rescaled).
This is the most defensible comparison point because it does not rely
on the c_1/c_3 split estimate.

Replace the ESTIMATE rows with direct transcription if Epelbaum Table 2
is republished at finer granularity, or if one digs into Epelbaum's
raw H-matrix code to extract the per-LEC split.
