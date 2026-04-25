// ===============================================================
// 抽取自仓库 [origin]: CPP/make_permutation_matrix.cpp
// 行号区段：720..886
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	size_t	P123_row_idx = 0;
	size_t	P123_col_idx = 0;
	size_t  pq_WP_idx	 = 0;
	bool    WP_overlap	 = false;
	double  P123_element = 0;

	int L_2N	   = 0;
	int L_1N	   = 0;
	int L_2N_prime = 0;
	int L_1N_prime = 0;
	
	double Gtilde_subarray [Gtilde_subarray_size];

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

					// ONLY USE S-WAVE (handy for Malfliet-Tjon debugging)
					//if (L_2N!=0 || L_2N_prime!=0){
					//	continue;
					//}
					
					if (production_run){
						calculate_Gtilde_subarray_polar(Gtilde_subarray,
												  	    &Atilde_store[alphap_idx*Nalpha*(Lmax+1) + alpha_idx*(Lmax+1)],
								   					    Nx, x_array,
								   					    Nphi,
								   					    &sin_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
								   					    &cos_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
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
	
							/* Unique index for current combination of WPs */
							size_t step_length_1 = (size_t) Np_WP*Nq_WP*Np_WP;
							size_t step_length_2 = (size_t) 	  Nq_WP*Np_WP;
							size_t step_length_3 = (size_t)		        Np_WP;
							size_t pq_WP_idx = (size_t) qp_idx_WP*step_length_1
									  	  	 		  + pp_idx_WP*step_length_2
									  	  	 		  +  q_idx_WP*step_length_3
									  	  	 		  +  p_idx_WP;
							WP_overlap = pq_WP_overlap_array[pq_WP_idx];

							/* Only calculate P123 if there is WP bin-overlap in Heaviside functions */
							if (WP_overlap){
								if (production_run){
									P123_element = calculate_P123_element_in_WP_basis_mod (Gtilde_subarray,
											   												(int) p_idx_WP,  (int) q_idx_WP,
											   												(int) pp_idx_WP, (int) qp_idx_WP,
											   												Np_WP,p_array_WP_bounds,
											   												Nq_WP,q_array_WP_bounds,
											   												Nx, x_array, wx_array,
											   												Nphi,
											   												&sin_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
											   												&cos_phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
											   												&wphi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi]);
								}
								else{
									P123_row_idx = alphap_idx*Nq_WP*Np_WP + qp_idx_WP*Np_WP +  pp_idx_WP;
									P123_col_idx = alpha_idx*Nq_WP*Np_WP + q_idx_WP*Np_WP + p_idx_WP;

									P123_element = cos(P123_row_idx)*cos(P123_col_idx);
								}
								
								if (P123_element!=0){
									P123_row_idx = alphap_idx*Nq_WP*Np_WP + qp_idx_WP*Np_WP +  pp_idx_WP;
									P123_col_idx = alpha_idx*Nq_WP*Np_WP + q_idx_WP*Np_WP + p_idx_WP;

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
										if ( P123_dim_omp>=tread_buffer_size ){

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
