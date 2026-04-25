// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：54..121
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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

