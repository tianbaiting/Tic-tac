// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：85..91
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
inline double regulator_gauss(double p, double q, double Lambda) noexcept
{
    const double num    = 4.0 * p * p + 3.0 * q * q;
    const double denom4 = 4.0 * Lambda * Lambda;
    const double a      = num / denom4;
    return std::exp(-a * a);
}
