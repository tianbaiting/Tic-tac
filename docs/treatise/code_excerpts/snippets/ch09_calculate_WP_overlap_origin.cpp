// ===============================================================
// 抽取自仓库 [origin]: CPP/make_permutation_matrix.cpp
// 行号区段：159..322
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
/* Calculates all overlapping bins <p'q'|pq> in wave-packet representation */
void calculate_WP_overlap(bool* pq_WP_overlap_array,
						  int   Np_WP, double *p_array_WP_bounds,
						  int   Nq_WP, double *q_array_WP_bounds,
						  int   Nx,    double* x_array,   double* wx_array,
						  int   Nphi,  double* phi_array, double* wphi_array){
	if (true){
	#pragma omp parallel
	{
		#pragma omp for
	for (size_t qp_idx_WP=0; qp_idx_WP<Nq_WP; qp_idx_WP++){
		double qp_l = q_array_WP_bounds[qp_idx_WP];
		double qp_u = q_array_WP_bounds[qp_idx_WP+1];
		for (size_t pp_idx_WP=0; pp_idx_WP<Np_WP; pp_idx_WP++){
			double pp_l = p_array_WP_bounds[pp_idx_WP];
			double pp_u = p_array_WP_bounds[pp_idx_WP+1];
	
			double phi_lower = atan(pp_l/qp_u);
			double phi_upper = atan(pp_u/qp_l);
			
			/* Create phi-mesh */
			calc_gauss_points ( &phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
							   &wphi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
							   phi_lower, phi_upper,
							   Nphi);
							   
			/* Verify that on-shell elements can exist for given phi-boundaries */
			for (size_t q_idx_WP=0; q_idx_WP<Nq_WP; q_idx_WP++){
				double q_l = q_array_WP_bounds[q_idx_WP];
				double q_u = q_array_WP_bounds[q_idx_WP+1];
				for (size_t p_idx_WP=0; p_idx_WP<Np_WP; p_idx_WP++){
					double p_l = p_array_WP_bounds[p_idx_WP];
					double p_u = p_array_WP_bounds[p_idx_WP+1];
	
					bool WP_overlap = false;
					/* Ensure possible phi boundaries */
					int hit_counter = 0;
					if (phi_lower<phi_upper){
						/* Search for on-shell elements */
						for (size_t phi_idx=0; phi_idx<Nphi; phi_idx++){
							double phi = phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi + phi_idx];
							double sin_phi = std::sin(phi);
							double cos_phi = std::cos(phi);
	
							double kmin = std::max(pp_l/sin_phi, qp_l/cos_phi);
							double kmax = std::min(pp_u/sin_phi, qp_u/cos_phi);
	
							/* Ensure possible k-boundaries */
							if (kmin<kmax){
								for (int x_idx=0; x_idx<Nx; x_idx++){
									double x = x_array[x_idx];
	
									double zeta_1 = pi1_tilde(sin_phi, cos_phi, x);
									double zeta_2 = pi2_tilde(sin_phi, cos_phi, x);
	
									double Heaviside_lower = std::max(p_l/zeta_1, q_l/zeta_2);
									double Heaviside_upper = std::min(p_u/zeta_1, q_u/zeta_2);
	
									/* Ensure overlapping Heaviside boundaries */
									if (Heaviside_lower<Heaviside_upper){
										double kpmin = std::max(kmin, Heaviside_lower);
										double kpmax = std::min(kmax, Heaviside_upper);
	
										/* Skip momentum-violating integral boundaries */
										if (kpmin<kpmax){
											WP_overlap = true;
											
											//if (std::abs(kpmax-kpmin)>1e-15){
											//	hit_counter += 1;
											//}
											break;
										}
									}
								}
							}
							if (WP_overlap){
								break;
							}
						}
					}
	
					/* Unique index for current combination of WPs */
					size_t step_length_1 = (size_t) Np_WP*Nq_WP*Np_WP;
					size_t step_length_2 = (size_t) 	  Nq_WP*Np_WP;
					size_t step_length_3 = (size_t)		        Np_WP;
					size_t pq_WP_idx = (size_t) qp_idx_WP*step_length_1
							  	  	 		  + pp_idx_WP*step_length_2
							  	  	 		  +  q_idx_WP*step_length_3
							  	  	 		  +  p_idx_WP;
					//std::cout << pq_WP_idx << " " << qp_idx_WP << " " << pp_idx_WP << " " << q_idx_WP << " " << p_idx_WP << std::endl;
					pq_WP_overlap_array[pq_WP_idx] = WP_overlap;
				}
			}
		}
	}
	}
	}
	else{
		#pragma omp parallel
		{
		#pragma omp for
		for (size_t qp_idx_WP=0; qp_idx_WP<Nq_WP; qp_idx_WP++){
			double qp_l = q_array_WP_bounds[qp_idx_WP];
			double qp_u = q_array_WP_bounds[qp_idx_WP+1];
			for (size_t pp_idx_WP=0; pp_idx_WP<Np_WP; pp_idx_WP++){
				double pp_l = p_array_WP_bounds[pp_idx_WP];
				double pp_u = p_array_WP_bounds[pp_idx_WP+1];

				double phi_lower = atan(pp_l/qp_u);
				double phi_upper = atan(pp_u/qp_l);

				/* Create phi-mesh */
				calc_gauss_points ( &phi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
								   &wphi_array[(qp_idx_WP*Np_WP + pp_idx_WP)*Nphi],
								   phi_lower, phi_upper,
								   Nphi);

				for (size_t q_idx_WP=0; q_idx_WP<Nq_WP; q_idx_WP++){
					double q_l = q_array_WP_bounds[q_idx_WP];
					double q_u = q_array_WP_bounds[q_idx_WP+1];
					for (size_t p_idx_WP=0; p_idx_WP<Np_WP; p_idx_WP++){
						double p_l = p_array_WP_bounds[p_idx_WP];
						double p_u = p_array_WP_bounds[p_idx_WP+1];
	
						double pi1_min = pi1_tilde(pp_l, qp_l, -1.0);
						double pi1_max = pi1_tilde(pp_u, qp_u, +1.0);
						bool p_in_pi1 = ( (pi1_min<=p_l && p_l<=pi1_max) || (pi1_min<=p_u && p_u<=pi1_max) );

						double pi2_min = pi2_tilde(pp_l, qp_l, +1.0);
						double pi2_max = pi2_tilde(pp_u, qp_u, -1.0);
						bool q_in_pi2 = ( (pi2_min<=q_l && q_l<=pi2_max) || (pi2_min<=q_u && q_u<=pi2_max) );

						bool WP_overlap = false;
						if ( p_in_pi1 && q_in_pi2 ){
							WP_overlap = true;
						}
					
						/* Unique index for current combination of WPs */
						size_t step_length_1 = (size_t) Np_WP*Nq_WP*Np_WP;
						size_t step_length_2 = (size_t) 	  Nq_WP*Np_WP;
						size_t step_length_3 = (size_t)		        Np_WP;
						size_t pq_WP_idx = (size_t) qp_idx_WP*step_length_1
								  	  	 		  + pp_idx_WP*step_length_2
								  	  	 		  +  q_idx_WP*step_length_3
								  	  	 		  +  p_idx_WP;
						//std::cout << pq_WP_idx << " " << qp_idx_WP << " " << pp_idx_WP << " " << q_idx_WP << " " << p_idx_WP << std::endl;
						pq_WP_overlap_array[pq_WP_idx] = WP_overlap;
					}
				}
			}
		}
		}
	}
	//double sparsity = 0.99881;
	//for (int i=0; i<Nq_WP*Nq_WP*Np_WP*Np_WP; i++){
	//	double prob_nnz = (double) rand() / RAND_MAX;
    //    if (prob_nnz>sparsity){
	//		pq_WP_overlap_array[i] = true;
	//	}
	//	else{
	//		pq_WP_overlap_array[i] = false;
	//	}
	//}
}
