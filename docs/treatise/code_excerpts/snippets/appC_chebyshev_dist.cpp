// ===============================================================
// 抽取自仓库 [current]: CPP/make_wp_states.cpp
// 行号区段：26..75
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void make_chebyshev_distribution(int N_WP,
								 double* boundary_array,
								 double scale,
								 double	sparseness_degree){
	double tan_term = 0;
	double boundary = 0;
	
	for (int i=1; i<N_WP+1; i++){
		tan_term = tan( (2*i-1)*M_PI/(4*N_WP) ); 

		boundary = scale*pow(tan_term, sparseness_degree);
		
		/* Momentum distribution */
		boundary_array[i] = boundary;
		/* Energy distribution */
		//boundary_array[i] = com_energy_to_com_q_momentum(boundary);
	}
	boundary_array[0] = 0.0;
}

void make_p_bin_grid(fwp_statespace& fwp_states, run_params run_parameters){
	if (run_parameters.p_grid_type=="chebyshev"){
		double scale 			 = run_parameters.chebyshev_s;
		double sparseness_degree = run_parameters.chebyshev_t;
		make_chebyshev_distribution(fwp_states.Np_WP, fwp_states.p_WP_array,
									scale,
									sparseness_degree);
	}
	else if (run_parameters.p_grid_type=="custom"){
		read_WP_boundaries_from_txt(fwp_states.p_WP_array, fwp_states.Np_WP, run_parameters.p_grid_filename);
	}
	else{
		raise_error("Unknown p-momentum gridtype specified.");
	}
}
void make_q_bin_grid(fwp_statespace& fwp_states, run_params run_parameters){
	if (run_parameters.q_grid_type=="chebyshev"){
		double scale 			 = run_parameters.chebyshev_s;
		double sparseness_degree = run_parameters.chebyshev_t;
		make_chebyshev_distribution(fwp_states.Nq_WP, fwp_states.q_WP_array,
									scale,
									sparseness_degree);
	}
	else if (run_parameters.q_grid_type=="custom"){
		read_WP_boundaries_from_txt(fwp_states.q_WP_array, fwp_states.Nq_WP, run_parameters.q_grid_filename);
	}
	else{
		raise_error("Unknown q-momentum gridtype specified.");
	}
}
