// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_permutation_matrix.cpp
// 行号区段：742..956
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	// [EN] The outer OpenMP loop walks row packets |X'> and computes the non-zero couplings to column packets |X>.
	// Sparse thread-local buffers are flushed to disk chunk-by-chunk so large production runs stay memory safe.
	// / [CN] 最外层 OpenMP 循环遍历的是行波包 |X'>，并计算其与列波包 |X> 的非零耦合；线程局部稀疏缓冲区会分块刷盘，
	// 以便大规模生产计算仍然保持内存可控。
	/* Start of P123 parallel calculation */
	printf("   - Initiating P123-matrix calculation ... \n");
	printf("     - Running OpenMP on %d threads \n", P123_omp_num_threads);
	fflush(stdout);
	auto timestamp_P123_tot_start = chrono::system_clock::now();
	
	std::vector<int> num_rows_calculated(P123_omp_num_threads, 0);
	
	omp_set_num_threads(P123_omp_num_threads);
	#pragma omp parallel
	{
	int thread_idx = omp_get_thread_num();

	//#pragma omp for
	//for (int row=0; row<P123_dense_dim; row+=stplngth){
	//	for (int col=0; col<P123_dense_dim; col+=stplngth){
	//		double prob_nnz = (double) rand() / RAND_MAX;
	//		int P123_dim_omp = P123_dim_array_omp[thread_idx];;
	//		double P123_element = (double) rand() / RAND_MAX;
	//		(P123_val_array_omp[thread_idx])[P123_dim_omp] = P123_element;
	//		(P123_row_array_omp[thread_idx])[P123_dim_omp] = row;
	//		(P123_col_array_omp[thread_idx])[P123_dim_omp] = col;
	//		P123_dim_omp += 1;
	//		P123_dim_array_omp[thread_idx] = P123_dim_omp;
	//		if ( P123_dim_omp>=current_array_dim ){  // This should occur a small amount of the time
	//			increase_sparse_array_size(&P123_val_array_omp[thread_idx], current_array_dim, sparse_step_length);
	//			increase_sparse_array_size(&P123_row_array_omp[thread_idx], current_array_dim, sparse_step_length);
	//			increase_sparse_array_size(&P123_col_array_omp[thread_idx], current_array_dim, sparse_step_length);
	//			current_array_dim += sparse_step_length;
	//		}
	//		if (P123_dim_omp>50806056/16){
	//			break;
	//		}
	//	}
	//}
	//}

	size_t	P123_row_idx = 0;
	size_t	P123_col_idx = 0;
	size_t  pq_WP_idx	 = 0;
	bool    WP_overlap	 = false;
	double  P123_element = 0;

	int L_2N	   = 0;
	int L_1N	   = 0;
	int L_2N_prime = 0;
	int L_1N_prime = 0;

	std::vector<double> Gtilde_subarray(Gtilde_subarray_size, 0.0);

	#pragma omp for
	/* <X_i'j'^alpha'| - loops (rows of P123) */
	for (size_t qp_idx_WP = 0; qp_idx_WP < Nq_WP; qp_idx_WP++){
		for (size_t pp_idx_WP = 0; pp_idx_WP < Np_WP; pp_idx_WP++){

			/* Progress printout by thread 0 */
			if (thread_idx==0){
				int num_rows_count = 0;
				for (size_t i=0; i<P123_omp_num_threads; i++){
					num_rows_count += num_rows_calculated[i];
				}

				auto timestamp_P123_tot_inter = chrono::system_clock::now();
				chrono::duration<double> time_P123_inter = timestamp_P123_tot_inter - timestamp_P123_tot_start;
				double P123_current_time = time_P123_inter.count();
				double P123_av_row_time  = P123_current_time/num_rows_count;
				double P123_est_completion_time = (P123_dense_dim-num_rows_count)*P123_av_row_time/3600.;

				printf("\r     - Calculated %d of %d rows. Av. time per row: %.3f s. Est. completion time: %.1f h", num_rows_count, P123_dense_dim, P123_av_row_time, P123_est_completion_time);
				fflush(stdout);
			}

			for (size_t alphap_idx = 0; alphap_idx < Nalpha; alphap_idx++){

				/* |X_ij^alpha> - loops (columns of P123) */
				for (size_t alpha_idx = 0; alpha_idx < Nalpha; alpha_idx++){

					L_2N = L_2N_array[alphap_idx];
					L_1N = L_1N_array[alphap_idx];

					L_2N_prime = L_2N_array[alpha_idx];
					L_1N_prime = L_1N_array[alpha_idx];
					const size_t phi_base_idx = phi_packet_index(qp_idx_WP, pp_idx_WP, 0, Np_WP, Nphi);

					// ONLY USE S-WAVE (handy for Malfliet-Tjon debugging)
					//if (L_2N!=0 || L_2N_prime!=0){
					//	continue;
					//}

					if (production_run){
						// [EN] For a fixed row packet, all candidate column packets share the same transformed angular
						// mesh. We therefore build Gtilde once per (alpha',alpha,q',p') block and reuse it for every
						// overlapping (q,p) cell. / [CN] 对固定行波包而言，所有候选列波包共享同一套变换后的角网格；因此这里对
						// 每个 (alpha',alpha,q',p') 块只构造一次 Gtilde，然后在所有重叠的 (q,p) 单元上复用。
						calculate_Gtilde_subarray_polar(Gtilde_subarray.data(),
														&Atilde_store[alphap_idx*Nalpha*(Lmax+1) + alpha_idx*(Lmax+1)],
														Nx, x_array,
														Nphi,
														&sin_phi_array[phi_base_idx],
														&cos_phi_array[phi_base_idx],
														L_2N, L_2N_prime, max_L12,
														L_1N, L_1N_prime, max_l3,
														two_J_3N,
														ClebschGordan_data,
														gsl_Plm_1_array, gsl_Plm_1_stplen,
														&gsl_Plm_2_array[(qp_idx_WP*Np_WP*Nphi*Nx + pp_idx_WP*Nphi*Nx)*gsl_Plm_2_stplen], gsl_Plm_2_stplen,
														&gsl_Plm_3_array[(qp_idx_WP*Np_WP*Nphi*Nx + pp_idx_WP*Nphi*Nx)*gsl_Plm_3_stplen], gsl_Plm_3_stplen,
														prefac_L_array,
														prefac_l_array,
														two_jmax_Clebsch);
					}

					for (size_t q_idx_WP = 0; q_idx_WP < Nq_WP; q_idx_WP++){
						for (size_t p_idx_WP = 0; p_idx_WP < Np_WP; p_idx_WP++){

							size_t pq_WP_idx = wp_overlap_index(qp_idx_WP,
																pp_idx_WP,
																q_idx_WP,
																p_idx_WP,
																Nq_WP,
																Np_WP);
							WP_overlap = pq_WP_overlap_array[pq_WP_idx];

							/* Only calculate P123 if there is WP bin-overlap in Heaviside functions */
							if (WP_overlap){
								if (production_run){
									// [EN] This is the actual packet average <X'|P|X>: once the geometric overlap and
									// angular kernel are known, the remaining work is a finite quadrature over the packet
									// cell. / [CN] 这里才是真正的波包平均矩阵元 <X'|P|X>：在几何重叠和角核都已知之后，
									// 剩下的就是对该波包单元做有限维求积。
									P123_element = calculate_P123_element_in_WP_basis_mod(Gtilde_subarray.data(),
																						   (int) p_idx_WP,  (int) q_idx_WP,
																						   (int) pp_idx_WP, (int) qp_idx_WP,
																						   Np_WP, p_array_WP_bounds,
																						   Nq_WP, q_array_WP_bounds,
																						   Nx, x_array, wx_array,
																						   Nphi,
																						   &sin_phi_array[phi_base_idx],
																						   &cos_phi_array[phi_base_idx],
																						   &wphi_array[phi_base_idx]);
								}
								else{
									P123_row_idx = packet_state_index(alphap_idx, qp_idx_WP, pp_idx_WP, Nq_WP, Np_WP);
									P123_col_idx = packet_state_index(alpha_idx, q_idx_WP, p_idx_WP, Nq_WP, Np_WP);

									P123_element = cos(P123_row_idx)*cos(P123_col_idx);
								}

								if (P123_element!=0){
									P123_row_idx = packet_state_index(alphap_idx, qp_idx_WP, pp_idx_WP, Nq_WP, Np_WP);
									P123_col_idx = packet_state_index(alpha_idx, q_idx_WP, p_idx_WP, Nq_WP, Np_WP);

									if (use_dense_format){
										size_t P123_mat_idx = P123_row_idx*P123_dense_dim + P123_col_idx;
										(*P123_val_dense_array)[P123_mat_idx] = P123_element;
									}
									else{
										size_t P123_dim_omp = P123_dim_array_omp[thread_idx];

										(P123_val_array_omp[thread_idx])[P123_dim_omp] = P123_element;
										(P123_row_array_omp[thread_idx])[P123_dim_omp] = P123_row_idx;
										(P123_col_array_omp[thread_idx])[P123_dim_omp] = P123_col_idx;

										/* Increment sparse dimension (num of non-zero elements) */
										P123_dim_omp += 1;
										P123_dim_array_omp[thread_idx] = P123_dim_omp;

										/* Check if we have filled buffer array.
										 * If so, store array to disk and continue calculations.
										 * Each thread stores to their own file (filename includes thread index)
										 * such that there are no race hazard. Each write-to-disk creates a new file */
										/* Thread File Count (TFC) */
										int current_TFC = max_TFC_array[thread_idx];
										if (P123_dim_omp>=tread_buffer_size){

											std::string thread_filename = generate_subarray_file_name(two_J_3N, P_3N,
																									  Np_WP, Nq_WP,
																									  J_2N_max,
																									  thread_idx,
																									  current_TFC,
																									  P123_folder);

											/* Store array */
											store_sparse_permutation_matrix_for_3N_channel_h5(P123_val_array_omp[thread_idx],
																							  P123_row_array_omp[thread_idx],
																							  P123_col_array_omp[thread_idx],
																							  P123_dim_omp,
																							  Np_WP, p_array_WP_bounds,
																							  Nq_WP, q_array_WP_bounds,
																							  pw_states,
																							  thread_filename,
																							  false);

											/* Re-set sparse-dimension. Old sparse-elements will be rewritten and/or not stored */
											P123_dim_array_omp[thread_idx] = 0;

											/* Increment the number of files stored for current thread */
											max_TFC_array[thread_idx] += 1;
										}
									}
								}
							}
						}
					}
				}

				num_rows_calculated[thread_idx] += 1;
			}
		}
	}
	}
