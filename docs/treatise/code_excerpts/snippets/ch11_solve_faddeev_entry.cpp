// ===============================================================
// 抽取自仓库 [current]: src/core/faddeev_solver/solve_faddeev.cpp
// 行号区段：1825..1850
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void solve_faddeev_equations(cdouble*  U_array,
							 cdouble*  U_BU_array,
							 cdouble*  G_array,
							 double*   P123_sparse_val_array,
							 int*      P123_sparse_row_array,
							 size_t*   P123_sparse_col_array_csc,
							 size_t    P123_sparse_dim,
							 double*   V_WP_unco_array,
							 double*   V_WP_coup_array,
							 swp_statespace swp_states,
							 fwp_statespace fwp_states,
							 channel_os_indexing chn_os_indexing,
							 pw_3N_statespace pw_states,
							 std::string file_identification,
					         run_params run_parameters,
							 const three_nucleon_force_model* tnf){
	
	/* Make local pointers & variables for on-shell channel-indexing */
	int*   q_com_idx_array		= chn_os_indexing.q_com_idx_array;
	int*   deuteron_idx_array	= chn_os_indexing.deuteron_idx_array;
	size_t num_q_com			= (size_t) chn_os_indexing.num_T_lab;
	size_t num_deuteron_states	= (size_t) chn_os_indexing.num_deuteron_states;
	/* Make local pointers & variables for SWP-statespace */
	size_t    Nq_WP					= (size_t) swp_states.Nq_WP;
	size_t    Np_WP					= (size_t) swp_states.Np_WP;
	double*   C_WP_unco_array		= swp_states.C_SWP_unco_array;
