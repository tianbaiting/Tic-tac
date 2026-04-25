// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_N2LO_3NF.h
// 行号区段：134..163
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	double W1_contact(int alpha_r, int alpha_c,
					  double p_r, double q_r,
					  double p_c, double q_c,
					  const pw_3N_statespace& pw_states) const
	{
		// Short-circuit when the LEC is off (pure 2NF runs still build this object).
		if (m_c_E == 0.0) return 0.0;

		// Rank-0 recoupling: pair scalar (sigma2.sigma3)(tau2.tau3) with full
		// Kronecker selection on pair and spectator quantum numbers.
		const double recoup = recoupling_3nf_scalar(
			pw_states.L_2N_array[alpha_r], pw_states.S_2N_array[alpha_r],
			pw_states.J_2N_array[alpha_r], pw_states.T_2N_array[alpha_r],
			pw_states.L_1N_array[alpha_r], pw_states.two_J_1N_array[alpha_r],
			pw_states.two_J_3N_array[alpha_r],
			pw_states.L_2N_array[alpha_c], pw_states.S_2N_array[alpha_c],
			pw_states.J_2N_array[alpha_c], pw_states.T_2N_array[alpha_c],
			pw_states.L_1N_array[alpha_c], pw_states.two_J_1N_array[alpha_c],
			pw_states.two_T_3N_array[alpha_r]);
		if (recoup == 0.0) return 0.0;

		// Squared-Gaussian regulator per E2002 eq. (3.19), applied to bra and ket.
		const double f_bra = chiral_3nf::regulator_gauss(p_r, q_r, m_Lambda);
		const double f_ket = chiral_3nf::regulator_gauss(p_c, q_c, m_Lambda);

		// LEC kernel: -½ c_E / (fπ⁴ Λ_χ) (E2002 eq. 2.10, spectator-1 component).
		const double lec = chiral_3nf::kernel_contact(m_c_E, m_fpi4, m_Lambda_chi);

		return recoup * lec * f_bra * f_ket;
	}
