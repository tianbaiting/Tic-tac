// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：230..313
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================

// -----------------------------------------------------------------------------
// 2PE c₁/c₃ rank-2 spatial kernel at fixed x (companion to kernel_2pe_c1c3):
//
//     kernel_rank2(x) = lec_bracket × 2 × pp² / (sqrt(30) × (q₂² + m_π²)²)
//
// This provides the SPATIAL RANK-2 part of the partial-wave matrix element for
// the operator [σ₂⊗σ₃]₂·[q₂⊗q₃]₂ × (tau₂·tau₃) × lec_bracket/fπ² × propagators.
//
// Derivation:
//   The rank-2 spatial tensor [q₂⊗q₃]₂ projected onto pair partial waves
//   L_2N → L_2N+2 (or L_2N → L_2N-2) produces a radial factor via the pair
//   orbital angular integral. Taking p̂_bra = ẑ (reference axis) and
//   integrating over the ket pair direction p̂':
//
//     ⟨L_r,0||[q₂⊗q₃]₂||L_c,0⟩ = (2/√30) × p_c² × δ_{L_r,L_c+2}  +  ...
//
//   where p_c = |p'| (ket pair momentum magnitude). This factor (2/√30) is
//   designed to cancel with the √30 reduced matrix element of [σ₂⊗σ₃]₂
//   that is included in `recoupling_3nf_rank2`, so the full combined factor
//   from recoupling × kernel is:
//
//     (τ₂·τ₃) × √30 × hat × 9j × CG × [lec × 2 × p_c² / (√30 × prop²)]
//     = (τ₂·τ₃) × 2 × hat × 9j × CG × lec × p_c² / prop²
//
//   The factor 2 comes from the pair orbital normalization in the angular
//   integral (2π from the azimuthal average, divided by the 4π pair solid
//   angle, gives the 2 via the spherical harmonic overlap).
//
// Momentum arguments (all in fm⁻¹):
//   p   = bra pair momentum (p_r)          [used for dp2 only]
//   q   = bra spectator momentum (q_r)
//   pp  = ket pair momentum (p_c = p')     [the D-wave factor; pp² is the rank-2 spatial factor]
//   qp  = ket spectator momentum (q_c = q')
//   x   = cos(q̂·q̂') ∈ [-1, +1]
//   m_pi_fm  = pion mass in fm⁻¹
//   c1_fm    = c₁ LEC in fm
//   c3_fm    = c₃ LEC in fm
//   fpi_fm   = f_π in fm
//
// Returns: spatial rank-2 radial kernel in fm⁵.
//
// Note: the (gA/2fπ)² global prefactor, the regulator product f_R(p,q)f_R(p',q'),
// the ⅓ rank-0 coefficient (which does NOT apply here for rank-2), and the
// τ₂·τ₃ isospin factor are NOT included — they live in the outer W1_2pe caller.
// Only the lec_bracket = (-4c₁mπ²+2c₃q₂q₃)/fπ² factor IS included here (same
// as in kernel_2pe_c1c3) because it is the same LEC combination for both ranks.
// -----------------------------------------------------------------------------
inline double kernel_2pe_c1c3_rank2(double p, double q,
                                     double pp, double qp,
                                     double x,
                                     double m_pi_fm,
                                     double c1_fm, double c3_fm,
                                     double fpi_fm) noexcept
{
    // Azimuthal-averaged scalar momentum pieces (same as rank-0 kernel).
    const double dp2  = p * p + pp * pp;                       // ⟨|Δp|²⟩
    const double dq2  = q * q + qp * qp - 2.0 * q * qp * x;  // |Δq|²(x)
    const double q2sq = dp2 + 0.25 * dq2;                      // ⟨q₂²⟩ = ⟨q₃²⟩
    const double q2q3 = dp2 - 0.25 * dq2;                      // ⟨q₂·q₃⟩

    const double mp2    = m_pi_fm * m_pi_fm;
    const double fpi_sq = fpi_fm * fpi_fm;

    // Pion-propagator product: 1 / [(q₂² + m_π²)(q₃² + m_π²)].
    const double prop = 1.0 / ((q2sq + mp2) * (q2sq + mp2));

    // LEC combination from [G2010] eq. (18) — identical to rank-0:
    //     [-4 c₁ m_π² + 2 c₃ (q₂·q₃)] / fπ²
    const double lec_bracket = (-4.0 * c1_fm * mp2 + 2.0 * c3_fm * q2q3) / fpi_sq;

    // Spatial rank-2 radial factor:
    //   pp² from the pair orbital angular integral ⟨L_r||[q₂⊗q₃]₂||L_c⟩
    //   Factor 2/√30: the 2 comes from the pair orbital normalization
    //   (from the 2π/4π azimuthal integral × Y₂ overlap); the 1/√30 cancels
    //   with the √30 reduced matrix element in recoupling_3nf_rank2.
    const double rank2_spatial = 2.0 / std::sqrt(30.0) * (pp * pp);

    return fourier_norm_3nf * (lec_bracket * rank2_spatial * prop);
}

}  // namespace chiral_3nf

#endif  // CHIRAL_3NF_PW_KERNELS_H
