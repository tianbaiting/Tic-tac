// ===============================================================
// 抽取自仓库 [current]: src/core/potential/make_potential_matrix.cpp
// 行号区段：400..429
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
						/* Write element to potential matrix V_array */
						if (coupled_matrix){
							if (coupled_via_L_2N){
								V_WP_coup_array[idx_V_WP_upper_left]  = extract_potential_element_from_array(J_r-1, J_r-1, J_r, S_r, coupled_matrix, V_WP_elements);
								V_WP_coup_array[idx_V_WP_upper_right] = extract_potential_element_from_array(J_r-1, J_r+1, J_r, S_r, coupled_matrix, V_WP_elements);
								V_WP_coup_array[idx_V_WP_lower_left]  = extract_potential_element_from_array(J_r+1, J_r-1, J_r, S_r, coupled_matrix, V_WP_elements);
								V_WP_coup_array[idx_V_WP_lower_right] = extract_potential_element_from_array(J_r+1, J_r+1, J_r, S_r, coupled_matrix, V_WP_elements);
							}
							else if (coupled_via_T_3N){
								V_WP_coup_array[idx_V_WP_upper_left]  = V_WP_elements[2];
								V_WP_coup_array[idx_V_WP_upper_right] = V_WP_elements[3];
								V_WP_coup_array[idx_V_WP_lower_left]  = V_WP_elements[4];
								V_WP_coup_array[idx_V_WP_lower_right] = V_WP_elements[5];
							}
							//if (J_r==1 and L_r==0){
							//    std::cout << idx_bin_r <<" "<< idx_bin_c <<" "<< V_WP_coup_array[idx_V_WP_upper_left] << std::endl;
							//    std::cout << idx_bin_r <<" "<< idx_bin_c <<" "<< V_WP_coup_array[idx_V_WP_upper_right] << std::endl;
							//    std::cout << idx_bin_r <<" "<< idx_bin_c <<" "<< V_WP_coup_array[idx_V_WP_lower_left] << std::endl;
							//    std::cout << idx_bin_r <<" "<< idx_bin_c <<" "<< V_WP_coup_array[idx_V_WP_lower_right] << std::endl;
							//}
						}
						else{
							V_WP_unco_array[idx_V_WP_uncoupled]   = extract_potential_element_from_array(L_r, L_c, J_r, S_r, coupled_matrix, V_WP_elements);
						}
					}
				}
				//}
			}
		}
	}
