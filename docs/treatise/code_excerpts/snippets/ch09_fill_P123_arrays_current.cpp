// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_permutation_matrix.cpp
// 行号区段：1379..1491
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void fill_P123_arrays(double** P123_sparse_val_array,
					  int**    P123_sparse_row_array,
					  int**    P123_sparse_col_array,
					  size_t&  P123_sparse_dim,
					  bool     production_run,
					  fwp_statespace fwp_states,
					  pw_3N_statespace pw_states,
					  run_params run_parameters,
					  std::string P123_folder){
	
	int     Np_WP 			  = fwp_states.Np_WP;
	int     Nq_WP 			  = fwp_states.Nq_WP;
	double* p_array_WP_bounds = fwp_states.p_WP_array;
	double* q_array_WP_bounds = fwp_states.q_WP_array;

	int J_2N_max = pw_states.J_2N_max;
	int two_J_3N = pw_states.two_J_3N_array[0];
	int P_3N 	 = pw_states.P_3N_array[0];

	int Nphi = run_parameters.Nphi;
	int Nx   = run_parameters.Nx;
	
	/* Default filename for current chn_3N - used for storage and reading P123 */
	std::string P123_filename = run_parameters.P123_folder + "/" + "P123_sparse_JP_"
										+ to_string(two_J_3N) + "_" + to_string(P_3N)
										+ "_Np_" + to_string(Np_WP) + "_Nq_" + to_string(Nq_WP)
										+ "_J2max_" + to_string(J_2N_max) + ".h5";

	if (run_parameters.calculate_and_store_P123){

		// [EN] This is the full "construct once, reuse many times" path from the WPCD workflow: identify the allowed
		// packet couplings, build P123, then persist it for later multi-energy production runs. / [CN] 这里走的是 WPCD
		// 工作流里的“先构造一次、后面多次复用”路径：先识别允许的波包耦合，再构造 P123，最后存盘供之后的多能量生产计算重复使用。
		/* Calculate momentum conservation enforced by P123 */
		printf("Precalculating momentum conservation enforced by P123 ... \n");
		fflush(stdout);
		double* x_array  = new double [Nx];
		double* wx_array = new double [Nx];
		gauss(x_array, wx_array, Nx);

		double*  phi_array  = new double [(size_t)Nq_WP*Np_WP*Nphi];
		double*  wphi_array = new double [(size_t)Nq_WP*Np_WP*Nphi];

		bool* pq_WP_overlap_array = new bool [Nq_WP*Nq_WP*Np_WP*Np_WP];

		calculate_WP_overlap(pq_WP_overlap_array,
							 Np_WP, p_array_WP_bounds,
							 Nq_WP, q_array_WP_bounds,
							 Nx,   x_array,   wx_array,
							 Nphi, phi_array, wphi_array);

		printf("Calculating P123 ... \n");
		auto timestamp_P123_calc_start = chrono::system_clock::now();
		calculate_permutation_matrices_for_all_3N_channels(P123_sparse_val_array,
														   P123_sparse_row_array,
														   P123_sparse_col_array,
														   P123_sparse_dim,
														   run_parameters.production_run,
														   Np_WP, p_array_WP_bounds,
														   Nq_WP, q_array_WP_bounds,
														   Nx, x_array, wx_array,
														   Nphi, phi_array, wphi_array,
														   pq_WP_overlap_array,
														   J_2N_max,
														   pw_states,
														   run_parameters,
														   run_parameters.P123_folder);
		auto timestamp_P123_calc_end = chrono::system_clock::now();
		chrono::duration<double> time_P123_calc = timestamp_P123_calc_end - timestamp_P123_calc_start;
		printf(" - Done. Time used: %.6f\n", time_P123_calc.count());

		printf("Storing P123 to h5 ... \n");
		auto timestamp_P123_store_start = chrono::system_clock::now();
		store_sparse_permutation_matrix_for_3N_channel_h5(*P123_sparse_val_array,
														  *P123_sparse_row_array,
														  *P123_sparse_col_array,
														  P123_sparse_dim,
														  Np_WP, p_array_WP_bounds,
														  Nq_WP, q_array_WP_bounds,
														  pw_states,
														  P123_filename,
														  true);
		auto timestamp_P123_store_end = chrono::system_clock::now();
		chrono::duration<double> time_P123_store = timestamp_P123_store_end - timestamp_P123_store_start;
		printf(" - Done. Time used: %.6f\n", time_P123_store.count());

		delete [] x_array;
		delete [] wx_array;
		delete [] phi_array;
		delete [] wphi_array;
		delete [] pq_WP_overlap_array;
	}
	else if (run_parameters.solve_faddeev){
		// [EN] Once the packet basis is fixed, later runs can skip the expensive permutation build and load the
		// precomputed sparse matrix directly. / [CN] 一旦波包基固定，后续运行就可以跳过昂贵的置换矩阵构造，直接读取
		// 预先计算好的稀疏矩阵。
		printf("Reading P123 from h5 ... \n");

		auto timestamp_P123_read_start = chrono::system_clock::now();
		read_sparse_permutation_matrix_for_3N_channel_h5(P123_sparse_val_array,
														 P123_sparse_row_array,
														 P123_sparse_col_array,
														 P123_sparse_dim,
														 Np_WP, p_array_WP_bounds,
														 Nq_WP, q_array_WP_bounds,
														 pw_states,
														 P123_filename,
														 true);
		auto timestamp_P123_read_end = chrono::system_clock::now();
		chrono::duration<double> time_P123_read = timestamp_P123_read_end - timestamp_P123_read_start;
		printf(" - Done. Time used: %.6f\n", time_P123_read.count());
	}
}
