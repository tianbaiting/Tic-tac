// ===============================================================
// 抽取自仓库 [current]: src/core/resolvent/make_resolvent.cpp
// 行号区段：133..232
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void calculate_resolvent_array_in_SWP_basis(cdouble* G_array,
											double   E,
											swp_statespace swp_states,
											pw_3N_statespace pw_states,
											run_params run_parameters){

	// [EN] The full AGS/Faddeev iteration only needs the diagonal entries of G in the SWP basis, so this routine
	// fills one complex number per (alpha, q-bin, p-bin). The later Neumann/Padé solver multiplies by these entries
	// as a cheap elementwise step between rescattering iterations. / [CN] 完整的 AGS/Faddeev 迭代在 SWP 基中只需要 G 的对角元，因此这里为每个 (alpha, q-bin, p-bin) 填充一个复数；后续 Neumann/Padé 求解器在每次再散射迭代之间只需做一次廉价的逐元素乘法。
	
	int 	Np_WP			 = swp_states.Np_WP;
	int     Nq_WP			 = swp_states.Nq_WP;
	double* e_SWP_unco_array = swp_states.e_SWP_unco_array;
	double* e_SWP_coup_array = swp_states.e_SWP_coup_array;
	double* q_WP_array		 = swp_states.q_WP_array;

	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* two_T_3N_array = pw_states.two_T_3N_array;

	/* Loop over states along resolvent diagonal */
	for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
	
		int L = L_2N_array[idx_alpha];
		int S = S_2N_array[idx_alpha];
		int J = J_2N_array[idx_alpha];
		int T = T_2N_array[idx_alpha];

		int two_T_3N = two_T_3N_array[idx_alpha];

		// [EN] The channel resolvent is diagonal in the SWP basis, but each alpha state must still select the correct
		// interacting p-packet branch before the analytic cell average can be evaluated. / [CN] 通道分辨算符在 SWP 基中
		// 虽然是对角的，但每个 alpha 态仍需先选中正确的相互作用 p 波包分支，才能计算解析的单元平均。
		double* e_SWP_array_ptr = select_swp_energy_branch(L,
															   S,
															   J,
															   T,
															   two_T_3N,
															   Np_WP,
															   e_SWP_unco_array,
															   e_SWP_coup_array,
															   run_parameters);

		//printf("%d %d %d %d \n", L, S, J, T);

		/* p-momentum index loop */
		for (int idx_p_bin=0; idx_p_bin<Np_WP; idx_p_bin++){

			/* Upper and lower boundaries of current p-bin (expressed in energy) */
			double e_bin_lower = e_SWP_array_ptr[idx_p_bin    ];
			double e_bin_upper = e_SWP_array_ptr[idx_p_bin + 1];

			//printf("%.4f\n", e_bin_lower);

			/* Bound state check, given by p_bin_lower if bound state exists */
			bool bound_state_exists = false;
			double Eb = 0;
			if (e_bin_lower<0){
				Eb = e_bin_lower;
				bound_state_exists = true;
			}

			/* q-momentum index loop */            
			for (int idx_q_bin=0; idx_q_bin<Nq_WP; idx_q_bin++){

				/* Upper and lower boundaries of current q-bin */
				double q_bin_lower = q_WP_array[idx_q_bin    ];
				double q_bin_upper = q_WP_array[idx_q_bin + 1];

				cdouble R = {0, 0};
				cdouble Q = {0, 0};
				if (bound_state_exists){    // Calculate bound-continuum (BC) resolvent part R
					R = resolvent_bound_continuum(E, Eb,
												  q_bin_upper,
												  q_bin_lower);
				}
				else{                       // Calculate continuum-continuum (CC) resolvent part Q
					Q = resolvent_continuum_continuum(E, Eb,
													  q_bin_upper,
													  q_bin_lower,
													  e_bin_upper,
													  e_bin_lower);
					//std::cout << Q << std::endl;
				}

				/* Use identical indexing as used in permutation matrix */
				int G_idx = idx_alpha*Nq_WP*Np_WP + idx_q_bin*Np_WP + idx_p_bin;
				G_array[G_idx] = R + Q;

				//if (e_bin_lower<0){
				//	std::cout << R << " " << Q << std::endl;
				//}
				//std::cout << R << " " << Q << std::endl;
			}
		}
	}
}
