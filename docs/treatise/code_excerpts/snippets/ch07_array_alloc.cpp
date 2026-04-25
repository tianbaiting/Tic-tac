// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_pw_symm_states.cpp
// 行号区段：248..284
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	/* Append the state space size as well - this allows for simpler indexing */
	chn_idx_temp.push_back(Nalpha_temp);
	
	/* Write number of states found to input integer */
	pw_states.Nalpha = Nalpha_temp;     

	/* Write number of 3N channels found to input integer */
	int N_chn_3N = N_chn_3N_temp;
	
	/* Allocate arrays to input array pointers */       
	pw_states.J_2N_max			= J_2N_max;          
	pw_states.L_2N_array		= new int [pw_states.Nalpha];
	pw_states.S_2N_array		= new int [pw_states.Nalpha];
	pw_states.J_2N_array		= new int [pw_states.Nalpha];
	pw_states.T_2N_array		= new int [pw_states.Nalpha];
	pw_states.L_1N_array		= new int [pw_states.Nalpha];
	pw_states.two_J_1N_array	= new int [pw_states.Nalpha];
	pw_states.two_J_3N_array	= new int [pw_states.Nalpha];
	pw_states.two_T_3N_array	= new int [pw_states.Nalpha];
	pw_states.P_3N_array		= new int [pw_states.Nalpha];
	pw_states.chn_3N_idx_array  = new int [N_chn_3N+1];
	pw_states.N_chn_3N			= N_chn_3N;

	//*chn_3N_idx_array_ptr = new int [N_chn_3N+1];

	/* Write temporary vector contents to newly allocated arrays */
	std::copy( L_2N_temp.begin(), L_2N_temp.end(), pw_states.L_2N_array );
	std::copy( S_2N_temp.begin(), S_2N_temp.end(), pw_states.S_2N_array );
	std::copy( J_2N_temp.begin(), J_2N_temp.end(), pw_states.J_2N_array );
	std::copy( T_2N_temp.begin(), T_2N_temp.end(), pw_states.T_2N_array );
	std::copy( L_1N_temp.begin(), L_1N_temp.end(), pw_states.L_1N_array );
	std::copy( two_J_1N_temp.begin(), two_J_1N_temp.end(), pw_states.two_J_1N_array );
	std::copy( two_J_3N_temp.begin(), two_J_3N_temp.end(), pw_states.two_J_3N_array );
	std::copy( two_T_3N_temp.begin(), two_T_3N_temp.end(), pw_states.two_T_3N_array );
	std::copy( P_3N_temp.begin(), P_3N_temp.end(), pw_states.P_3N_array );
	std::copy( chn_idx_temp.begin(), chn_idx_temp.end(), pw_states.chn_3N_idx_array );
}
