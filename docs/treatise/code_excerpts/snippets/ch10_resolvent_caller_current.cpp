// ===============================================================
// 抽取自仓库 [current]: src/main.cpp
// 行号区段：454..473
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
				/* Start of code segment for resolvent matrix (diagonal array) construction */

				// [EN] Once the SWP spectrum is known, the channel resolvent is diagonal in this basis; step 8 of the
				// workflow therefore reduces to filling one complex propagator value per (alpha, q, p) cell and energy.
				// / [CN] 一旦 SWP 谱已知，通道分辨算符在该基中就是对角的；因此工作流的第 8 步只需为每个 (alpha,q,p) 单元
				// 和能量填入一个复传播子值。
				/* Resolvent array */
				size_t dense_dim = pw_substates.Nalpha * swp_states.Np_WP * swp_states.Nq_WP;
				std::vector<cdouble> G_array(dense_dim * chn_os_indexing.num_T_lab);
				
				printf(" - Constructing 3N resolvents ... \n");
				for (int j=0; j<chn_os_indexing.num_T_lab; j++){
					double E = solve_config.E_com_array[j] + swp_states.E_bound;
					calculate_resolvent_array_in_SWP_basis(&G_array[j*dense_dim],
															E,
															swp_states,
															pw_substates,
															run_parameters);
				}
				printf("   - Done \n");
