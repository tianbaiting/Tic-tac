// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：278..309
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
