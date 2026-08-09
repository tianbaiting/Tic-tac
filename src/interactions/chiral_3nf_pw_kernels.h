#ifndef CHIRAL_3NF_PW_KERNELS_H
#define CHIRAL_3NF_PW_KERNELS_H

#include <cmath>

// ============================================================================
// [EN] Scalar momentum-space kernels for the partial-wave decomposition of the
//      chiral N2LO three-nucleon force.
//
//      These helpers factor out the pure *momentum-space* pieces (pion
//      propagators, LEC prefactors, regulator) so that the outer partial-wave
//      machinery (recoupling coefficients, spin/isospin matrix elements,
//      x-integration) can stay clean. All functions are side-effect free,
//      header-only, and inline.
//
//      Formula sources (transcribed in
//      tools/check_3nf_normalization/formula_reference.md):
//        - [E2002] Epelbaum et al., PRC 66 (2002) 064001 — eqs. (2.10), (3.19),
//                  App. A eq. (A-2).
//        - [G2010] Golak et al., EPJA 43 (2010) 241 — eqs. (17), (18), (22).
//        - [H2015] Hebeler et al., PRC 91 (2015) 044001 — eqs. (13)–(15).
//
//      Conventions (see §5 of formula_reference.md):
//        - Jacobi momenta `p`, `q` in fm⁻¹ (spectator-1 convention).
//        - After fixing p̂ = ẑ and φ_q = 0, the only non-trivial angular
//          variable in the scalar partial-wave projection is
//          x ≡ cos(q̂, q̂'). Azimuthal averages are absorbed implicitly.
//        - Momentum transfers:
//            Δp = p' − p,    Δq = q' − q
//            q₁ = Δq                       (transferred to spectator 1)
//            q₂ = Δp + Δq/2                (to pair-particle 2)
//            q₃ = Δp − Δq/2                (to pair-particle 3)
//          With p̂ = ẑ and the azimuthal average over the p'-direction
//          orthogonal to the x-integration, cross-terms Δp·Δq drop out at
//          the scalar (rank-0) level, leaving:
//            |Δq|² = q² + q'² − 2 q q' x       ← [E2002] eq. (A-2)
//            ⟨|Δp|²⟩ = p² + p'²
//            ⟨q₂²⟩ = ⟨q₃²⟩ = ⟨|Δp|²⟩ + |Δq|²/4
//            ⟨q₂·q₃⟩ = ⟨|Δp|²⟩ − |Δq|²/4
//
//      Unit conventions (see §Units of formula_reference.md):
//        - momenta (p, q, p', q', m_π, Λ): fm⁻¹.
//        - fπ in fm (converted from MeV via fπ/ħc).
//        - Λ_χ in fm (converted from MeV via Λ_χ/ħc).
//        - c_E dimensionless; c_D dimensionless; c_1, c_3 in fm (converted
//          from GeV⁻¹ via ×ħc/1000 = 0.19733).
//
// [CN] 手征 N2LO 三体核力分波分解的标量动量空间核函数。
//      这些辅助函数将纯动量空间部分（π 传播子、LEC 前因子、正规化因子）
//      分离出来，使得外层分波机制（重耦合系数、自旋/同位旋矩阵元、x 积分）
//      可以保持清晰。所有函数无副作用，只含头文件，全部为 inline。
// ============================================================================

namespace chiral_3nf {

// -----------------------------------------------------------------------------
// Partial-wave Fourier normalization constant, 1 / (2π)³ = 1 / (8π³).
//
// This mirrors the 2NF convention in src/interactions/chiral_LO_internal.cpp:59
// where every two-body matrix element is multiplied by 1/(8π³) as part of the
// Tic-tac partial-wave basis Fourier convention. Diagnosed via the Epelbaum
// target comparison: the code/ref magnitude ratio was ~246× without this
// factor, essentially equal to (2π)³ = 248.05. Applied to every non-trivial
// momentum-space kernel below (kernel_contact, kernel_1pe_contact,
// kernel_2pe_c1c3). *Not* applied to the dimensionless regulator.
// -----------------------------------------------------------------------------
constexpr double fourier_norm_3nf = 1.0 / (8.0 * M_PI * M_PI * M_PI);

// -----------------------------------------------------------------------------
// Squared-Gaussian regulator (per [E2002] eq. 3.19):
//
//     f_R(p, q; Λ) = exp( - ((4 p² + 3 q²) / (4 Λ²))² )
//
// Input (all in fm⁻¹):
//   p           — pair relative momentum
//   q           — spectator momentum
//   Lambda      — 3NF cut-off Λ
//
// Returns: dimensionless regulator ∈ (0, 1].
//
// Note: the *squared* outer argument is essential; the earlier plain-Gaussian
// form exp(-(p²+¾q²)/Λ²) is too soft and was the cause of the ~10² overshoot
// in the c_E normalization diagnostic (see §Regulator of formula_reference.md).
// -----------------------------------------------------------------------------
inline double regulator_gauss(double p, double q, double Lambda) noexcept
{
    const double num    = 4.0 * p * p + 3.0 * q * q;
    const double denom4 = 4.0 * Lambda * Lambda;
    const double a      = num / denom4;
    return std::exp(-a * a);
}

// -----------------------------------------------------------------------------
// 3N contact scalar LEC factor (Navrátil/Witała convention):
//
//     V^(1)_cont prefactor = +½ c_E / (f_π⁴ Λ_χ)
//
// Sign convention note: [E2002] eq. 2.10 writes V^(1)_cont = -½ E Σ(τ_j·τ_k),
// but the code LEC values are in the Navrátil/Witała convention (matching the
// chiral_N2LO 2NF LECs that ship with Tic-tac and the Witała triton
// benchmarks). This convention absorbs the overall sign into the LEC
// definition — see epelbaum_reference.md lines 26-30. Diagnostic confirms:
// with -½ the magnitude matches Epelbaum Table 2 but the sign is flipped;
// with +½ the sign matches. See formula_reference.md §1.5 for the sign
// analysis and the commit fixing this (post-Task-5 triangulation).
//
// Input:
//   c_E         — dimensionless LEC (Navrátil/Witała convention)
//   fpi4_fm     — f_π⁴ in fm⁻⁴ (i.e. (fπ/ħc)⁴ with fπ in MeV)
//   Lambda_chi  — chiral breaking scale Λ_χ in fm⁻¹
//
// Returns: scalar in fm⁵ (momentum-independent LEC factor multiplying
//          (τ₂·τ₃) × f_R(p',q') f_R(p,q) × fourier_norm_3nf in the full
//          matrix element).
// -----------------------------------------------------------------------------
inline double kernel_contact(double c_E, double fpi4_fm, double Lambda_chi) noexcept
{
    // 1/(8π³) mirrors chiral_LO_internal.cpp:59 Fourier convention.
    return fourier_norm_3nf * (+0.5 * c_E / (fpi4_fm * Lambda_chi));
}

// -----------------------------------------------------------------------------
// 1PE-contact (c_D) scalar pion propagator at fixed x (per [E2002] eq. A-2):
//
//     kernel(x) = 1 / ( Q²(x) + m_π² )
//
//  with  Q²(x) = q² + q'² − 2 q q' x  =  |Δq|²
//
// This is the bare pion propagator appearing in §2 of formula_reference.md
// for the 1PE part of the c_D 1PE-contact term. In [E2002] eq. A-2 the
// full integrand is (Q^{k₁} / [Q² (Q² + m_π²)]); this function provides the
// (Q² + m_π²)⁻¹ factor, leaving the Q^{k₁}/Q² spatial-multipole reduction
// to the caller (Task 3 will wrap it in the Legendre expansion).
//
// Input (all in fm⁻¹):
//   p, q        — bra-side Jacobi momenta (unused at this scalar level; kept
//                 in the signature so later tasks can refine the kernel once
//                 p-direction effects enter via the rank-2 pieces)
//   pp, qp      — ket-side Jacobi momenta (p', q')
//   x           — cos(q̂, q̂') ∈ [-1, +1]
//   m_pi_fm     — pion mass in fm⁻¹
//
// Returns: scalar in fm² (units of 1/momentum²).
// -----------------------------------------------------------------------------
inline double kernel_1pe_contact(double /*p*/, double q,
                                 double /*pp*/, double qp,
                                 double x, double m_pi_fm) noexcept
{
    const double Q2  = q * q + qp * qp - 2.0 * q * qp * x;  // |Δq|²
    const double mp2 = m_pi_fm * m_pi_fm;
    // 1/(8π³) mirrors chiral_LO_internal.cpp:59 Fourier convention.
    return fourier_norm_3nf * (1.0 / (Q2 + mp2));
}

// -----------------------------------------------------------------------------
// 2PE c₁/c₃ scalar (rank-0) kernel at fixed x (per §3 and §4 of
// formula_reference.md; [G2010] eq. 18):
//
//     kernel(x) = [(-4 c₁ m_π² / fπ²) + (2 c₃ (q₂·q₃) / fπ²)]
//                  × (q₂·q₃) / [ (q₂² + m_π²) (q₃² + m_π²) ]
//
// Here we keep only the *rank-0* angular decomposition:
//
//     (σ₂·q₂)(σ₃·q₃)  →  ⅓ (σ₂·σ₃) (q₂·q₃)        [rank-0, §3.2 / §4.2]
//
// The σ₂·σ₃ pair eigenvalue, τ₂·τ₃ isospin factor, the outer (gA/2fπ)²,
// the ⅓ rank-0 coefficient, and the regulator product f_R(p',q')·f_R(p,q)
// are *not* included here — they live in the outer recoupling helper
// (`recoupling_3nf_scalar`) called by Tasks 3–5.
//
// The x-dependent scalar momentum quantities, consistent with the [E2002]
// eq. (A-2) azimuthal reduction (p̂ = ẑ, φ_q = 0, averaging over φ_{p'}):
//
//     dp2   = |Δp|²   ≈  p² + p'²           (azimuthal-averaged)
//     dq2   = |Δq|²   =  q² + q'² − 2qq'x
//     q2sq  = q₂²    =  dp2 + dq2/4
//     q3sq  = q₃²    =  dp2 + dq2/4         (degenerate with q₂² at scalar level)
//     q2q3  = q₂·q₃  =  dp2 − dq2/4
//
// The `p·p'` cross term drops out under the φ_{p'} azimuthal average for the
// rank-0 piece; the surviving cross term Δp·Δq is likewise zero-average.
// The rank-2 tensor piece (which does *not* vanish under azimuthal average
// once recoupled) is handled separately in Task 5 via
// `recoupling_3nf_rank2`.
//
// Input (all in fm⁻¹ except c_1, c_3):
//   p, q        — bra-side Jacobi momenta
//   pp, qp      — ket-side Jacobi momenta (p', q')
//   x           — cos(q̂, q̂') ∈ [-1, +1]
//   m_pi_fm     — pion mass in fm⁻¹
//   c1_fm       — c_1 LEC in fm  (= c_1[GeV⁻¹] × ħc[MeV·fm] / 1000)
//   c3_fm       — c_3 LEC in fm  (= c_3[GeV⁻¹] × ħc[MeV·fm] / 1000)
//   fpi_fm      — f_π in fm      (= fπ[MeV] / ħc[MeV·fm])
//
// Returns: scalar in fm⁵.
//
// Note: the explicit fπ dependence in the two c_i terms comes from [G2010]
// eq. (18): F_1 = (gA/2fπ)² / [(q₂²+mπ²)(q₃²+mπ²)] × [-4c₁mπ²/fπ² + 2c₃ q₂·q₃/fπ²].
// The (gA/2fπ)² global prefactor is *not* included here (pulled out by
// the caller, uniform across all pieces of V^(1)_2π).
// -----------------------------------------------------------------------------
inline double kernel_2pe_c1c3(double p, double q,
                              double pp, double qp,
                              double x,
                              double m_pi_fm,
                              double c1_fm, double c3_fm,
                              double fpi_fm) noexcept
{
    // Azimuthal-averaged scalar momentum pieces (see header comment).
    const double dp2  = p * p + pp * pp;                 // ⟨|Δp|²⟩
    const double dq2  = q * q + qp * qp - 2.0 * q * qp * x;  // |Δq|²(x)
    const double q2sq = dp2 + 0.25 * dq2;                // ⟨q₂²⟩ = ⟨q₃²⟩
    const double q2q3 = dp2 - 0.25 * dq2;                // ⟨q₂·q₃⟩

    const double mp2    = m_pi_fm * m_pi_fm;
    const double fpi_sq = fpi_fm * fpi_fm;

    // Pion-propagator product: 1 / [(q₂² + m_π²)(q₃² + m_π²)].
    const double prop = 1.0 / ((q2sq + mp2) * (q2sq + mp2));

    // Bracketed LEC combination from [G2010] eq. (18):
    //     [-4 c₁ m_π² + 2 c₃ (q₂·q₃)] / fπ²
    const double lec_bracket = (-4.0 * c1_fm * mp2 + 2.0 * c3_fm * q2q3) / fpi_sq;

    // Scalar (rank-0) spin-reduction carries a (q₂·q₃) factor from the
    // ⅓(σ₂·σ₃)(q₂·q₃) identity — see §3.2.
    // 1/(8π³) mirrors chiral_LO_internal.cpp:59 Fourier convention.
    return fourier_norm_3nf * (lec_bracket * q2q3 * prop);
}

}  // namespace chiral_3nf

#endif  // CHIRAL_3NF_PW_KERNELS_H
