// ===============================================================
// 抽取自仓库 [current]: src/core/faddeev_solver/solve_faddeev.h
// 行号区段：25..65
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
class three_nucleon_force_model;

// [EN] Bundles all data the column-computation hot path needs for the 3NF contribution so that
// existing function signatures stay manageable. When tnf->enabled()==false (null object) the entire
// 3NF branch is skipped via a single test. / [CN] 把列计算热路径所需的全部 3NF 数据打包到一个结构体中，
// 避免已有函数签名膨胀。当 tnf->enabled()==false（null 对象）时，整个 3NF 分支通过一次判断跳过。
struct tnf_kernel_context {
	const three_nucleon_force_model* tnf;
	const pw_3N_statespace*          pw_states;
	const double*                    p_WP_array;   // WP boundaries, size Np_WP+1
	const double*                    q_WP_array;   // WP boundaries, size Nq_WP+1
	double**                         CT_RM_array;  // C^T row-major = C column-major, [Nalpha*Nalpha]
	double                           w1_scale;     // Overall scale factor applied to W^(1) output (diagnostic knob)
};

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
							 const three_nucleon_force_model* tnf = nullptr);

/* Create array of pointers to C^T matrices for product (C^T)PVC in row-major format */
void create_CT_row_maj_3N_pointer_array(double** CT_RM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										int  	 num_2N_unco_states,
										int  	 num_2N_coup_states,
										size_t   Np_WP,
										pw_3N_statespace pw_states,
					         			run_params run_parameters);
