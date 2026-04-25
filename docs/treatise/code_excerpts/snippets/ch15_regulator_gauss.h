// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：70..92
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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

