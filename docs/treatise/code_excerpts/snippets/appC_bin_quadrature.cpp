// ===============================================================
// 抽取自仓库 [current]: CPP/make_wp_states.cpp
// 行号区段：77..105
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void make_p_bin_quadrature_grids(fwp_statespace& fwp_states){
	int 	Np_WP		= fwp_states.Np_WP;
	double* p_WP_array	= fwp_states.p_WP_array;
	int 	Np_per_WP	= fwp_states.Np_per_WP;
	double* p_array		= fwp_states.p_array;
	double* wp_array	= fwp_states.wp_array;

	for (int idx_bin=0; idx_bin<Np_WP; idx_bin++){
		double bin_lower_bound = p_WP_array[idx_bin];
		double bin_upper_bound = p_WP_array[idx_bin+1];

		double*  p_array_ptr =  &p_array[idx_bin*Np_per_WP];
		double* wp_array_ptr = &wp_array[idx_bin*Np_per_WP];
		for (int idx_p=0; idx_p<Np_per_WP; idx_p++){
			gauss(p_array_ptr, wp_array_ptr, Np_per_WP);
			updateRange_a_b(p_array_ptr, wp_array_ptr, bin_lower_bound, bin_upper_bound, Np_per_WP);
		}
	}
}
void make_q_bin_quadrature_grids(fwp_statespace& fwp_states){
	int 	Nq_WP		= fwp_states.Nq_WP;
	double* q_WP_array	= fwp_states.q_WP_array;
	int 	Nq_per_WP	= fwp_states.Nq_per_WP;
	double* q_array		= fwp_states.q_array;
	double* wq_array	= fwp_states.wq_array;

	for (int idx_bin=0; idx_bin<Nq_WP; idx_bin++){
		double bin_lower_bound = q_WP_array[idx_bin];
		double bin_upper_bound = q_WP_array[idx_bin+1];
