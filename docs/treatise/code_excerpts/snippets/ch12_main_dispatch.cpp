// ===============================================================
// 抽取自仓库 [current]: src/main.cpp
// 行号区段：524..555
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
				/* Start of code segment for storing on-shell U-matrix solutions */

				std::string U_mat_filename_t = run_parameters.output_folder + "/" + "U_PW_elements"
																			+ file_identification
																			+ "_PSI_" + std::to_string(idx_param_set)
																			+ ".txt";

				store_U_matrix_elements_txt(U_array.data(),
											solve_config,
											chn_os_indexing,
											run_parameters,
											swp_states,
											pw_substates,
											U_mat_filename_t);
				
				if (run_parameters.include_breakup_channels){
					std::string U_mat_BU_filename_t = run_parameters.output_folder + "/" + "U_PW_breakup_elements"
																				+ file_identification
																				+ "_PSI_" + std::to_string(idx_param_set)
																				+ ".txt";

					store_U_BU_matrix_elements_txt(U_BU_ptr,
												solve_config,
												chn_os_indexing,
												run_parameters,
												fwp_states,
												swp_states,
												pw_substates,
												U_mat_BU_filename_t);
				}

				/* End of code segment for storing on-shell U-matrix solutions */
