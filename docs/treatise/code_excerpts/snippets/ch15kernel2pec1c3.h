// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：155..229
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
