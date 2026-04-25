// ===============================================================
// 抽取自仓库 [current]: src/utils/auxiliary.cpp
// 行号区段：660..746
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
double calculate_P123_element_in_WP_basis_mod (double* Gtilde_subarray,
											   int  p_idx_WP, int  q_idx_WP,
											   int pp_idx_WP, int qp_idx_WP,
											   int Np_WP,     double *p_array_WP_bounds,
											   int Nq_WP,     double *q_array_WP_bounds,
											   int Nx, double *x_array, double *wx_array,
											   int Nphi,
											   double* sin_phi_subarray,
											   double* cos_phi_subarray,
											   double* wphi_subarray){
	
	/* Bin-boundaris (short-hand notation: l/u="lower"/"upper", pp="p-prime") */
	/* Ket bin boundaries */
	double  p_l = p_array_WP_bounds[ p_idx_WP];
	double  p_u = p_array_WP_bounds[ p_idx_WP + 1];
	double  q_l = q_array_WP_bounds[ q_idx_WP];
	double  q_u = q_array_WP_bounds[ q_idx_WP + 1];
	/* Bra bin boundaries */
	double pp_l = p_array_WP_bounds[pp_idx_WP];
	double pp_u = p_array_WP_bounds[pp_idx_WP + 1];
	double qp_l = q_array_WP_bounds[qp_idx_WP];
	double qp_u = q_array_WP_bounds[qp_idx_WP + 1];

	/* Variables used in integration (see project documentation) */
	double x       = 0;
	double wx      = 0;
	double sin_phi = 0;
	double cos_phi = 0;
	double wphi    = 0;
	double zeta_1  = 0;
	double zeta_2  = 0;
	double kmin    = 0;
	double kmax    = 0;
	double Gtilde  = 0;

	int idx_i = qp_idx_WP*Np_WP + pp_idx_WP;
	int idx_j =  q_idx_WP*Np_WP +  p_idx_WP;

	/* Loop over quadrature points */
	double integral_sum = 0;
	for (int x_idx=0; x_idx<Nx; x_idx++){
		x  =  x_array[x_idx];
		wx = wx_array[x_idx];
		for (int phi_idx=0; phi_idx<Nphi; phi_idx++){
			sin_phi = sin_phi_subarray[phi_idx];
			cos_phi = cos_phi_subarray[phi_idx];
			//cos_phi = sin_phi_subarray[phi_idx];
			//sin_phi = cos_phi_subarray[phi_idx];
			wphi    = wphi_subarray[phi_idx];

			zeta_1 = pi1_tilde(sin_phi, cos_phi, x);
			zeta_2 = pi2_tilde(sin_phi, cos_phi, x);

			double kmin_array[] = {pp_l/sin_phi, qp_l/cos_phi, p_l/zeta_1, q_l/zeta_2};
			double kmax_array[] = {pp_u/sin_phi, qp_u/cos_phi, p_u/zeta_1, q_u/zeta_2};

			kmin = *std::max_element(kmin_array, kmin_array+4);
			kmax = *std::min_element(kmax_array, kmax_array+4);

			//if (idx_i==63 && idx_j==58){
			//	std::cout << "qp_lo: " << qp_l << " qp_hi: " << qp_u << std::endl;
			//	std::cout << "pp_lo: " << pp_l << " pp_hi: " << pp_u << std::endl;
			//	std::cout << "q_lo:  " << q_l << "  q_hi:  " << q_u << std::endl;
			//	std::cout << "p_lo:  " << p_l << "  p_hi:  " << p_u << std::endl;
			//	std::cout << "zeta1: " << zeta_1 << "  zeta2: " << zeta_2 << std::endl;
			//	std::cout <<"x: "<<x_idx<<": "<<x<<" | phi: "<<phi_idx<<": "<<asin(sin_phi)<<" | Q: " << kmin << " " << kmax << std::endl;
			//}

			/* Skip momentum-violating integral boundaries */
			if (kmin>kmax){
				continue;
			}

			Gtilde = Gtilde_subarray[phi_idx*Nx + x_idx];

			integral_sum += wx * wphi * cos_phi*sin_phi * Gtilde * (kmax*kmax - kmin*kmin) / (2*zeta_1*zeta_2);
		}
	}

	/* Normalization */
	double norm_WP = q_normalization( q_l,  q_u)
					*p_normalization( p_l,  p_u)
					*q_normalization(qp_l, qp_u)
					*p_normalization(pp_l, pp_u);

	return integral_sum / norm_WP;
}
