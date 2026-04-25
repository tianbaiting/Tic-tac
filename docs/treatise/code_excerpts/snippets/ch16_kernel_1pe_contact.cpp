// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：145..153
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
inline double kernel_1pe_contact(double /*p*/, double q,
                                 double /*pp*/, double qp,
                                 double x, double m_pi_fm) noexcept
{
    const double Q2  = q * q + qp * qp - 2.0 * q * qp * x;  // |Δq|²
    const double mp2 = m_pi_fm * m_pi_fm;
    // 1/(8π³) mirrors chiral_LO_internal.cpp:59 Fourier convention.
    return fourier_norm_3nf * (1.0 / (Q2 + mp2));
}
