// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_N2LO_3NF.h
// 行号区段：83..98
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	double W1_element(int alpha_r, int alpha_c,
					  double p_r, double q_r,
					  double p_c, double q_c,
					  const pw_3N_statespace& pw_states) const override
	{
		// 3N conserved quantum numbers
		if (pw_states.two_J_3N_array[alpha_r] != pw_states.two_J_3N_array[alpha_c]) return 0.0;
		if (pw_states.two_T_3N_array[alpha_r] != pw_states.two_T_3N_array[alpha_c]) return 0.0;
		if (pw_states.P_3N_array[alpha_r]     != pw_states.P_3N_array[alpha_c])     return 0.0;

		double result = 0.0;
		result += W1_contact(alpha_r, alpha_c, p_r, q_r, p_c, q_c, pw_states);
		result += W1_1pe_contact(alpha_r, alpha_c, p_r, q_r, p_c, q_c, pw_states);
		result += W1_2pe(alpha_r, alpha_c, p_r, q_r, p_c, q_c, pw_states);
		return result;
	}
