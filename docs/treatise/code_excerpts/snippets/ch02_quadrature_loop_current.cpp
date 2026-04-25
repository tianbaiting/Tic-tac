// ===============================================================
// 抽取自仓库 [current]: src/core/potential/make_potential_matrix.cpp
// 行号区段：240..330
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
				int idx_V_WP_lower_left  = 0;
				int idx_V_WP_lower_right = 0;
				//#pragma for
				for (int idx_bin_r=0; idx_bin_r<Np_WP; idx_bin_r++){
					/* WP normalization */
					double N_r = norm_p_array[idx_bin_r];

					/* Column p-momentum index loop */
					for (int idx_bin_c=0; idx_bin_c<Np_WP; idx_bin_c++){
						/* WP normalization */
						double N_c = norm_p_array[idx_bin_c];

						/* Potential matrix indexing
						 * Indexing format: (channel index)*(num rows)*(num columns) + (row index)*(row length) + (column index) */
						if (coupled_matrix){
							idx_V_WP_upper_left  = chn_idx_V_coup*4*Np_WP*Np_WP +  idx_bin_r         *2*Np_WP + idx_bin_c;
							idx_V_WP_upper_right = chn_idx_V_coup*4*Np_WP*Np_WP +  idx_bin_r         *2*Np_WP + idx_bin_c + Np_WP;
							idx_V_WP_lower_left  = chn_idx_V_coup*4*Np_WP*Np_WP + (idx_bin_r + Np_WP)*2*Np_WP + idx_bin_c;
							idx_V_WP_lower_right = chn_idx_V_coup*4*Np_WP*Np_WP + (idx_bin_r + Np_WP)*2*Np_WP + idx_bin_c + Np_WP;
						
						}
						else{
							idx_V_WP_uncoupled   = chn_idx_V_unco*Np_WP*Np_WP + idx_bin_r*Np_WP + idx_bin_c;
						}
						
						/* Reset quadrature summation array V_WP_elements */
						for (int idx_element=0; idx_element<6; idx_element++){
							V_WP_elements[idx_element] = 0;
						}

						/* Calculate matrix element, using either mid-point approximation or quadrature */
						if (run_parameters.midpoint_approx==true){
							
							/* Average momenta */
							double p_in  = p_array[idx_bin_c];
							double p_out = p_array[idx_bin_r];

							/* Weigth-function terms */
							double f_in  = fp_array[idx_bin_c];
							double f_out = fp_array[idx_bin_r];

							/* Momentum bin-width */
							double d_c   = p_WP_array[idx_bin_c + 1] - p_WP_array[idx_bin_c];
							double d_r   = p_WP_array[idx_bin_r + 1] - p_WP_array[idx_bin_r];

							/* Row and column indices of p_out and p_in, respectively */
							int i = idx_bin_r;
							int j = idx_bin_c;

							/* We create an isoscalar potential */
							if (isoscalar_type==0){ // Interaction must be np
								pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_np, V_IS_elements);
							}
							else if (isoscalar_type==1){ // Interaction can be either nn or np with IS symmetry-conservation
								pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_nn, V_nn_elements);
								pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_np, V_np_elements);
								for (int idx_element=0; idx_element<6; idx_element++){
									V_IS_elements[idx_element] = (2./3)*V_nn_elements[idx_element] + (1./3)*V_np_elements[idx_element];
								}
							}
							else if (isoscalar_type==2 || isoscalar_type==3){ // Interaction can be either np or nn with IS symmetry-breaking
								pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_nn, V_nn_elements);
								pot_ptr->V(i, j, p_in, p_out, coupled_model, S_r, J_r, T_r, Tz_np, V_np_elements);
								
								/* Extract 1S0-element, using coupling as according to tensor-force only (i.e. no coupling for 1S0) */
								double V_nn_1S0_element = extract_potential_element_from_array(L_r, L_c, J_r, S_r, false, V_nn_elements);
								double V_np_1S0_element = extract_potential_element_from_array(L_r, L_c, J_r, S_r, false, V_np_elements);

								/* Define our own, convenient ordering of V_IS_elements:
								 * idx 2: top left block		(isoscalar type 1)
								 * idx 3: top right block		(isoscalar type 2)
								 * idx 4: bottom left block		(isoscalar type 2)
								 * idx 5: bottom right block	(isoscalar type 3) */
								V_IS_elements[2] = 		 (2./3)*  V_nn_1S0_element + (1./3)*V_np_1S0_element;
								V_IS_elements[3] = (sqrt(2)/3.)* (V_nn_1S0_element - 	    V_np_1S0_element);
								V_IS_elements[4] = (sqrt(2)/3.)* (V_nn_1S0_element - 	    V_np_1S0_element);
								V_IS_elements[5] = 		 (1./3)*  V_nn_1S0_element + (2./3)*V_np_1S0_element;
								
							}
							else{
								raise_error("Unknown isoscalar potential encountered.");
							}
							
							/* We integrate into wave-packet potential and normalize */
							for (int idx_element=0; idx_element<6; idx_element++){
								V_WP_elements[idx_element] += p_in*p_out *f_in*f_out *d_r*d_c * V_IS_elements[idx_element] / (N_r*N_c);
							}
						}
						else{
							/* Loop over quadrature inside row and column cells */
							for (int idx_p_r=0; idx_p_r<Np_per_WP; idx_p_r++){
