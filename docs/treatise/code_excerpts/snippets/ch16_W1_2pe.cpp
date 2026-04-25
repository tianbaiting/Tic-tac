// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_N2LO_3NF.h
// 行号区段：295..394
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	double W1_2pe(int alpha_r, int alpha_c,
				  double p_r, double q_r,
				  double p_c, double q_c,
				  const pw_3N_statespace& pw_states) const
	{
		if (m_c1 == 0.0 && m_c3 == 0.0) return 0.0;

		// [EN] Rank-0 recoupling: (σ₂·σ₃)(τ₂·τ₃) pair eigenvalues, diagonal in all
		// pair and spectator quantum numbers. Returns 0 for off-diagonal α pairs.
		// Note: recoupling_3nf_scalar returns sigma*sigma × tau*tau WITHOUT the ⅓
		// rank-0 coefficient — we apply 1/3 explicitly in the overall coefficient below.
		const double recoup0 = recoupling_3nf_scalar(
			pw_states.L_2N_array[alpha_r], pw_states.S_2N_array[alpha_r],
			pw_states.J_2N_array[alpha_r], pw_states.T_2N_array[alpha_r],
			pw_states.L_1N_array[alpha_r], pw_states.two_J_1N_array[alpha_r],
			pw_states.two_J_3N_array[alpha_r],
			pw_states.L_2N_array[alpha_c], pw_states.S_2N_array[alpha_c],
			pw_states.J_2N_array[alpha_c], pw_states.T_2N_array[alpha_c],
			pw_states.L_1N_array[alpha_c], pw_states.two_J_1N_array[alpha_c],
			pw_states.two_T_3N_array[alpha_r]);

		// [EN] Rank-2 recoupling: (τ₂·τ₃) × √30 × hat × 9j × CG(L_c,0;2,0|L_r,0).
		// Non-zero for L_r - L_c = ±2 with S_r=S_c=1 (triplet pairs only).
		// Enables 3S1 ↔ 3D1 coupling and the associated sign flip.
		const double recoup2 = recoupling_3nf_rank2(
			pw_states.L_2N_array[alpha_r], pw_states.S_2N_array[alpha_r],
			pw_states.J_2N_array[alpha_r], pw_states.T_2N_array[alpha_r],
			pw_states.L_1N_array[alpha_r], pw_states.two_J_1N_array[alpha_r],
			pw_states.two_J_3N_array[alpha_r],
			pw_states.L_2N_array[alpha_c], pw_states.S_2N_array[alpha_c],
			pw_states.J_2N_array[alpha_c], pw_states.T_2N_array[alpha_c],
			pw_states.L_1N_array[alpha_c], pw_states.two_J_1N_array[alpha_c],
			pw_states.two_T_3N_array[alpha_r]);

		if (recoup0 == 0.0 && recoup2 == 0.0) return 0.0;

		// Squared-Gaussian regulator per E2002 eq. (3.19).
		const double f_bra = chiral_3nf::regulator_gauss(p_r, q_r, m_Lambda);
		const double f_ket = chiral_3nf::regulator_gauss(p_c, q_c, m_Lambda);

		// fπ in fm⁻¹ (needed for (gA/2fπ)² = gA²/(4fπ²) prefactor).
		// Computed from m_fpi4 = fπ⁴ stored in the class.
		const double fpi_fm = std::sqrt(std::sqrt(m_fpi4));

		// [EN] Rank-2 Hermitian symmetrization: compute the transposed recoupling
		// recoup2_T = recoupling with bra and ket swapped (alpha_c→alpha_r convention).
		// The spatial kernel kernel_2pe_c1c3_rank2 has a pp² = p_c² dependence that
		// makes W1(alpha_r,alpha_c,p_r,q_r,p_c,q_c) ≠ W1(alpha_c,alpha_r,p_c,q_c,p_r,q_r)
		// when only one side is evaluated. Explicit symmetrization:
		//   rank2_contribution = ½ [recoup2(r,c) × integ2(p_r,q_r,p_c,q_c)
		//                          + recoup2(c,r) × integ2(p_c,q_c,p_r,q_r)]
		// ensures the full W1 satisfies W1(r,c,p_r,q_r,p_c,q_c) = W1(c,r,p_c,q_c,p_r,q_r).
		const double recoup2_T = recoupling_3nf_rank2(
			pw_states.L_2N_array[alpha_c], pw_states.S_2N_array[alpha_c],
			pw_states.J_2N_array[alpha_c], pw_states.T_2N_array[alpha_c],
			pw_states.L_1N_array[alpha_c], pw_states.two_J_1N_array[alpha_c],
			pw_states.two_J_3N_array[alpha_c],
			pw_states.L_2N_array[alpha_r], pw_states.S_2N_array[alpha_r],
			pw_states.J_2N_array[alpha_r], pw_states.T_2N_array[alpha_r],
			pw_states.L_1N_array[alpha_r], pw_states.two_J_1N_array[alpha_r],
			pw_states.two_T_3N_array[alpha_c]);

		const bool need_rank2 = (recoup2 != 0.0 || recoup2_T != 0.0);

		// Gauss-Legendre x-quadrature for pion propagator kernels.
		// x = cos(q̂·q̂') ∈ [-1, +1]; Q²(x) = q² + q'² - 2qq'x = |Δq|²(x).
		double integ0 = 0.0, integ2 = 0.0, integ2_T = 0.0;
		for (int ix = 0; ix < N_GL; ++ix) {
			const double x  = m_gl_x[ix];
			const double wx = m_gl_w[ix];
			if (recoup0 != 0.0) {
				integ0 += wx * chiral_3nf::kernel_2pe_c1c3(
					p_r, q_r, p_c, q_c, x, m_mpi_fm, m_c1, m_c3, fpi_fm);
			}
			if (need_rank2) {
				// Direct: kernel with pp = p_c (ket pair momentum)
				const double k2 = chiral_3nf::kernel_2pe_c1c3_rank2(
					p_r, q_r, p_c, q_c, x, m_mpi_fm, m_c1, m_c3, fpi_fm);
				// Transposed: swap (p_r,q_r) ↔ (p_c,q_c) so pp = p_r (bra pair momentum)
				const double k2_T = chiral_3nf::kernel_2pe_c1c3_rank2(
					p_c, q_c, p_r, q_r, x, m_mpi_fm, m_c1, m_c3, fpi_fm);
				if (recoup2   != 0.0) integ2   += wx * k2;
				if (recoup2_T != 0.0) integ2_T += wx * k2_T;
			}
		}

		// Overall coefficient: (gA/2fπ)² = gA²/(4fπ²).
		// The fπ² here matches [G2010] eq. (18): (gA/2fπ)² prefactor in F₁.
		// Note: this is fπ² NOT fπ⁴ (the kernels already carry one fπ² inside lec_bracket).
		// Units: gA² dimensionless, fπ² in fm⁻², so coeff has units fm².
		const double coeff = m_gA * m_gA / (4.0 * fpi_fm * fpi_fm);

		// Rank-0 carries a ⅓ from the (σ₂·q₂)(σ₃·q₃) = ⅓(σ₂·σ₃)(q₂·q₃) decomposition.
		// Rank-2 does NOT carry ⅓ — it is the full [σ₂⊗σ₃]₂·[q₂⊗q₃]₂ piece.
		// The ½ symmetrization factor ensures Hermitian symmetry: each side contributes
		// half of the combined recoup2×integ2 + recoup2_T×integ2_T.
		return coeff * f_bra * f_ket
		       * ((1.0/3.0) * recoup0 * integ0
		          + 0.5 * (recoup2 * integ2 + recoup2_T * integ2_T));
	}
