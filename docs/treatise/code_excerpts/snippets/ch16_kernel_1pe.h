// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：100..160
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
