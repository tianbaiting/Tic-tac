// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：116..120
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
inline double kernel_contact(double c_E, double fpi4_fm, double Lambda_chi) noexcept
{
    // 1/(8π³) mirrors chiral_LO_internal.cpp:59 Fourier convention.
    return fourier_norm_3nf * (+0.5 * c_E / (fpi4_fm * Lambda_chi));
}
