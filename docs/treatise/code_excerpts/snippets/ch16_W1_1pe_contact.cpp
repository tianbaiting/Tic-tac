// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_N2LO_3NF.h
// 行号区段：195..251
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	double W1_1pe_contact(int alpha_r, int alpha_c,
						  double p_r, double q_r,
						  double p_c, double q_c,
						  const pw_3N_statespace& pw_states) const
	{
		if (m_c_D == 0.0) return 0.0;

		// Rank-0 recoupling: (1/3)(σ₁·σ₃)(τ₁·τ₃) in 3N Jj basis.
		// Selection rules: L_2N=L_2N'=0, L_1N=L_1N' (enforced inside helper).
		const double recoup0 = recoupling_3nf_1pe_ct_scalar(
			pw_states.L_2N_array[alpha_r], pw_states.S_2N_array[alpha_r],
			pw_states.J_2N_array[alpha_r], pw_states.T_2N_array[alpha_r],
			pw_states.L_1N_array[alpha_r], pw_states.two_J_1N_array[alpha_r],
			pw_states.two_J_3N_array[alpha_r],
			pw_states.L_2N_array[alpha_c], pw_states.S_2N_array[alpha_c],
			pw_states.J_2N_array[alpha_c], pw_states.T_2N_array[alpha_c],
			pw_states.L_1N_array[alpha_c], pw_states.two_J_1N_array[alpha_c],
			pw_states.two_T_3N_array[alpha_r]);

		// Rank-2 recoupling: [σ₁⊗σ₃]₂·Y₂(q̂) coefficient.
		// For l_1N=0 channels (dominant triton configs): returns 0 by CG selection rule.
		const double recoup2 = recoupling_3nf_1pe_ct_rank2(
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

		// Gauss-Legendre x-integration for the 1PE pion propagator.
		// Q²(x) = q² + q'² − 2qq'x  (|Δq|²  per E2002 eq. A-2)
		// rank-0 weight: P_0(x) = 1
		// rank-2 weight: P_2(x) = (3x²−1)/2
		double integ0 = 0.0, integ2 = 0.0;
		for (int ix = 0; ix < N_GL; ++ix) {
			const double x  = m_gl_x[ix];
			const double wx = m_gl_w[ix];
			const double k  = chiral_3nf::kernel_1pe_contact(p_r, q_r, p_c, q_c,
			                                                  x, m_mpi_fm);
			if (recoup0 != 0.0) integ0 += wx * k;               // P_0 = 1
			if (recoup2 != 0.0) integ2 += wx * k * (1.5*x*x - 0.5); // P_2(x)
		}

		// Overall coefficient: −g_A c_D / (8 f_π⁴ Λ_χ) × 2
		// (factor of 2 from summing j=2 and j=3 in the operator).
		const double coeff = -m_gA * m_c_D / (8.0 * m_fpi4 * m_Lambda_chi) * 2.0;

		return coeff * f_bra * f_ket * (recoup0 * integ0 + recoup2 * integ2);
	}
