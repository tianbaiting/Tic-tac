// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_3nf_pw_kernels.h
// 行号区段：122..153
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
