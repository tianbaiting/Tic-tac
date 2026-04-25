// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_permutation_matrix.h
// 行号区段：42..100
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
const int max_size_per_TFC_file = 1; // 1 GB per thread, per file

// [EN] Build the geometric overlap mask that discards packet pairs forbidden by the Jacobi-coordinate map behind the
// permutation operator before any expensive quadrature is attempted. / [CN] 在做任何昂贵求积之前，先构造几何重叠掩码，
// 把那些被置换算符背后的 Jacobi 坐标变换直接禁止的波包配对剔除掉。
void calculate_WP_overlap(bool* pq_WP_overlap_array,
						  int   Np_WP, double *p_array_WP_bounds,
						  int   Nq_WP, double *q_array_WP_bounds,
						  int   Nx,    double* x_array,   double* wx_array,
						  int   Nphi,  double* phi_array, double* wphi_array);

// [EN] Construct one J^pi block of the finite-dimensional permutation matrix P123. This is the expensive, but
// energy-independent, packet-space representation reused by all later Faddeev solves on the same basis. / [CN]
// 构造一个 J^pi 通道块上的有限维置换矩阵 P123；这一步代价高，但在固定基底下与能量无关，因此后续所有 Faddeev 求解都可复用。
void calculate_permutation_elements_for_3N_channel(double** P123_val_dense_array,
												   int*		max_TFC_array,
												   bool     use_dense_format,
												   bool     production_run,
												   int      Np_WP, double *p_array_WP_bounds,
												   int      Nq_WP, double *q_array_WP_bounds,
												   int      Nx, double* x_array, double* wx_array,
												   int      Nphi,
												   bool*    pq_WP_overlap_array,
												   int      J_2N_max,
												   pw_3N_statespace pw_states,
												   run_params run_parameters,
												   std::string P123_folder);

// [EN] Orchestrate the one-time P123 build or recovery flow, then merge the thread-local fragments into the global
// sparse COO object used by the solver. / [CN] 组织 P123 的一次性构造或恢复流程，然后把线程局部片段合并成求解器使用的
// 全局稀疏 COO 对象。
void calculate_permutation_matrices_for_all_3N_channels(double** P123_sparse_val_array,
														int**    P123_sparse_row_array,
														int**    P123_sparse_col_array,
														size_t&  P123_sparse_dim_array,
														bool     production_run,
														int      Np_WP, double *p_array_WP_bounds,
														int      Nq_WP, double *q_array_WP_bounds,
														int      Nx, double* x_array, double* wx_array,
														int      Nphi,
														bool*    pq_WP_overlap_array,
														int      J_2N_max,
														pw_3N_statespace pw_states,
														run_params run_parameters,
														std::string P123_folder);

// [EN] Public entry point used by the main workflow: either generate/store the sparse permutation matrix for the
// current packet basis or load it back for an actual Faddeev solve. / [CN] 主流程使用的公开入口：要么为当前波包基生成并
// 存储稀疏置换矩阵，要么在实际 Faddeev 求解时把它重新读回内存。
void fill_P123_arrays(double** P123_sparse_val_array,
						  int**    P123_sparse_row_array,
						  int**    P123_sparse_col_array,
					  size_t&  P123_sparse_dim,
					  bool     production_run,
					  fwp_statespace fwp_states,
					  pw_3N_statespace pw_states,
					  run_params run_parameters,
					  std::string P123_folder);

