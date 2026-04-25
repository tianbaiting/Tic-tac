// ===============================================================
// 抽取自仓库 [origin]: CPP/solve_faddeev.h
// 行号区段：25..38
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
							 channel_os_indexing chn_os_indexing,
							 pw_3N_statespace pw_states,
							 std::string file_identification,
					         run_params run_parameters);
