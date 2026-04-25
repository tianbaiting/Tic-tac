// ===============================================================
// 抽取自仓库 [current]: src/core/faddeev_solver/solve_faddeev.cpp
// 行号区段：1609..1679
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
		printf("       - Calculating Pade approximants PA[%ld,%ld]. \n", NM, NM); fflush(stdout);
		// [EN] Padé resummation is what turns a slowly convergent or even divergent Neumann history into stable
		// amplitudes. Each on-shell element is treated independently because different channels can converge at very
		// different rates. / [CN] Padé 重求和是把收敛缓慢甚至发散的 Neumann 历史转化为稳定振幅的关键。这里每个 on-shell
		// 元素独立处理，因为不同通道的收敛速度可能相差很大。
		/* Calculate Pade approximants (PA) for elastic amplitudes */
		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
				for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){

					size_t idx_NDOS = elastic_value_storage_index(idx_d_row,
															   idx_d_col,
															   idx_q_com,
															   num_deuteron_states,
															   num_q_com);

					/* Check and skip if we've already reached convergence for this on-shell element */
					if (pade_approximants_conv_array[idx_NDOS]==true){
						continue;
					}

					/* Calculate and append PA */
					cdouble PA = pade_approximant(&a_coeff_array[idx_NDOS*num_neumann_terms], NM, NM, 1);

					pade_approximants_array[idx_NDOS*(NM_max+1) + NM] = PA;
					
					/* See if we've reached convergence with this iteration */
					size_t idx_best_PA = 0;
					bool convergence_reached = false;

					/* Find minimum PA from previous calculations */
					double min_PA_diff = 1;
					for (int NM_prev=0; NM_prev<NM; NM_prev++){
						cdouble PA_prev = pade_approximants_array[idx_NDOS*(NM_max+1) + NM_prev];

						/* Calculate difference between PAs from previous PA-calculations */
						double PA_diff_prev = std::abs(PA_prev - pade_approximants_array[idx_NDOS*(NM_max+1) + NM_prev-1]);

						/* Ignore PA_diff_prev if numerically equal to the previous PA_diff, overwrite if smaller than min_PA_diff */
						if (PA_diff_prev<min_PA_diff && PA_diff_prev>1e-15){
							idx_best_PA = NM_prev;
							min_PA_diff = PA_diff_prev;
						}
					}
					/* See if current PA is better/worse than previous minimum */
					double PA_diff_curr = std::abs(PA - pade_approximants_array[idx_NDOS*(NM_max+1) + NM - 1]);
					
					/* Ignore PA_diff_curr if numerically equal to the previous PA_diff_curr, overwrite if smaller than min_PA_diff */
					if (PA_diff_curr<min_PA_diff && PA_diff_curr>1e-15){
						idx_best_PA = NM;
						min_PA_diff = PA_diff_curr;
					}

					/* Criterias for convergence¨
					 * If any are fulfilled, we set convergence to true for current on-shell element */
					bool convergence_criteria_0 = (NM==NM_max);													// Cannot go past max NM
					bool convergence_criteria_1 = (NM-idx_best_PA>4);												// If nothing better is found in the last 3 PAs, we assume we found the best
					bool convergence_criteria_2 = (min_PA_diff<1e-6*std::abs(pade_approximants_array[idx_NDOS*(NM_max+1) + idx_best_PA]));	// If the difference is less than the 4th significant digit, we assume "good enough"
					bool convergence_criteria_3 = (min_PA_diff<1e-7);												// If we are below single precision resolution, assume convergence
					
					if (convergence_criteria_0 ||
						convergence_criteria_1 ||
						convergence_criteria_2 ||
						convergence_criteria_3){
						pade_approximants_conv_array[idx_NDOS] = true;
						pade_approximants_idx_array[idx_NDOS]  = idx_best_PA;
						num_converged_elements += 1;
					}
				}
			}
		}
