// ===============================================================
// 抽取自仓库 [current]: src/core/faddeev_solver/solve_faddeev.cpp
// 行号区段：1304..1395
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	for (size_t NM=0; NM<NM_max+1; NM++){
		if (num_converged_elements==num_EL_A_vals){
			printf("     - Convergence reached for all on-shell elements! \n"); fflush(stdout);
			break;
		}

		if (NM!=0){
			printf("     - Working on Pade approximant P[N,M] for N=%ld, M=%ld \n",NM,NM); fflush(stdout);
		}
		
		size_t counter_array [100];
		for (size_t i=0; i<100; i++){
			counter_array[i] = 0;
		}

		/* Time-keeper array for parallel environment */
		double*  times_array = new double [3*num_threads];
		
		for (int n=2*NM-1; n<2*NM+1; n++){
			printf("       - Working on Neumann-terms for n=%d. \n", n); fflush(stdout);
			/* We've already done n=0 above */
			if (n<=0){
				continue;
			}
			
			/* Initialise time-profile variables */
			double time_resolvent        = 0;
			double time_CPVC_cols        = 0;
			double time_An_CPVC_multiply = 0;
			double time_neumann          = 0;

			double timestamp_neumann_start = omp_get_wtime();

			// [EN] The Neumann recursion is implemented exactly as in the lecture notes: first multiply the previous
			// term by the diagonal channel resolvent G, then apply the CPVC kernel to generate the next rescattering
			// contribution. / [CN] Neumann 递推严格按照讲稿里的顺序实现：先把上一阶乘上对角的通道分辨算符 G，再施加 CPVC
			// 核，生成下一阶再散射贡献。
			double timestamp_resolvent_start = omp_get_wtime();
			printf("       - Multiplying in resolvent with An. \n"); fflush(stdout);
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
					size_t idx_row_NDOS = idx_d_row*num_q_com + idx_q_com;

					/* Multiply An by G */
					for (size_t idx=0; idx<dense_dim; idx++){
						double re_An = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx];
						double im_An = im_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx];
						double re_G  = G_array[idx_q_com*dense_dim + idx].real();
						double im_G  = G_array[idx_q_com*dense_dim + idx].imag();
						re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx] = re_An*re_G - im_An*im_G;
						im_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx] = re_An*im_G + im_An*re_G;
					}
				}
			}
			double timestamp_resolvent_end    = omp_get_wtime();
			time_resolvent = timestamp_resolvent_end - timestamp_resolvent_start;

			// [EN] Once some on-shell amplitudes have converged, there is no value in pushing their full dense rows
			// through later rescattering steps. Compacting only the non-converged rows is a pure algebraic shortcut:
			// it preserves the exact iteration on the active rows while reducing GEMM cost. / [CN] 当部分 on-shell
			// 振幅已经收敛后，就没有必要再把它们对应的整条稠密行送入后续再散射步骤；这里只压缩未收敛的行是纯代数层面的加速，在保持活动行迭代完全一致的同时降低了 GEMM 成本。
			size_t num_non_conv_rows = 0;
			for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
				for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
					size_t idx_row_NDOS = elastic_row_storage_index(idx_d_row, idx_q_com, num_q_com);
					if (row_has_only_converged_targets(idx_d_row,
													   idx_q_com,
													   num_deuteron_states,
													   num_q_com,
													   pade_approximants_conv_array,
													   pade_approximants_BU_conv_array,
													   chn_os_indexing,
													   run_parameters.include_breakup_channels)==false){
						for (size_t i=0; i<dense_dim; i++){
							re_A_An_row_array_comp[num_non_conv_rows*dense_dim + i] = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + i];
							im_A_An_row_array_comp[num_non_conv_rows*dense_dim + i] = im_A_An_row_array_prev[idx_row_NDOS*dense_dim + i];
						}
						A_An_indexing_array[num_non_conv_rows] = idx_row_NDOS;
						num_non_conv_rows += 1;
					}
				}
			}
			
			printf("       - Calculating on-shell rows of A*K^n for n=%d. \n", n); fflush(stdout);
			if (keep_CPVC_in_mem==false){
				// [EN] The kernel columns are regenerated in chunks so the dense GEMM path can stream through the
				// active part of CPVC without materializing the full dense matrix. / [CN] 这里按块重建核列，这样稠密 GEMM
				// 路径就能流式处理 CPVC 的活动部分，而不必把整个稠密矩阵完整落在内存中。
				for (size_t idx_col_chunk=0; idx_col_chunk<num_col_chunks; idx_col_chunk++){

					double timestamp_CPVC_chunk_start = omp_get_wtime();

