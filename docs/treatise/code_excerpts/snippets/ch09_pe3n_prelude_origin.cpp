// ===============================================================
// 抽取自仓库 [origin]: CPP/make_permutation_matrix.cpp
// 行号区段：324..460
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void calculate_permutation_elements_for_3N_channel(double** P123_val_dense_array,
												   int*		max_TFC_array,
												   bool     use_dense_format,
												   bool     production_run,
												   int      Np_WP, double *p_array_WP_bounds,
												   int      Nq_WP, double *q_array_WP_bounds,
												   int      Nx, double* x_array, double* wx_array,
												   int      Nphi, double* phi_array, double* wphi_array,
												   bool*    pq_WP_overlap_array,
												   int      J_2N_max,
												   pw_3N_statespace pw_states,
												   run_params run_parameters,
												   std::string P123_folder){
	
	/* Make local pointers & variables */
	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* L_1N_array		= pw_states.L_1N_array;
	int* two_J_1N_array = pw_states.two_J_1N_array;
	int* two_T_3N_array	= pw_states.two_T_3N_array;
	int  two_J_3N    	= pw_states.two_J_3N_array[0];
	int  P_3N   	    = pw_states.P_3N_array[0];
	
	bool print_content = true;

	/* Notation change */
	int two_J = two_J_3N;
	//int two_T = two_T_3N;
	
	int Nx_Gtilde = Nx;
	int Jj_dim = Nalpha;

	int* L12_Jj    = L_2N_array;
	int* S12_Jj    = S_2N_array;
	int* J12_Jj    = J_2N_array;
	int* T12_Jj    = T_2N_array;
	int* l3_Jj     = L_1N_array;
	int* two_j3_Jj = two_J_1N_array;

	/* End of notation change */

	/* START OF OLD CODE SEGMENT WITH OLD VARIABLE-NOTATION */
	/* This code calculates the geometric function Gtilde_{alpha,alpha'}(p',q',x) as an array */

	// determine optimized Lmax: Lmax = max(get_L)+max(get_l)
	int max_L12 = 0;
	int max_l3 = 0;
	int max_J12 = 0;

	for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
		if (J12_Jj[alpha] > max_J12) max_J12 = J12_Jj[alpha];
		if (L12_Jj[alpha] > max_L12) max_L12 = L12_Jj[alpha];
		if (l3_Jj[alpha] > max_l3) max_l3 = l3_Jj[alpha];
	}

	// for F_local_matrix prestorage and F_interpolate
	int lmax = GSL_MAX_INT(max_l3, max_L12) + 3; // for C4 it is possible to couple l=lmax with THREE Y_{1}^{mu}
	
	//int l_interpolate_max = l_interpolate_max = 2 * (lmax - 3) + 3;
	//if (print_content){
	//	std::cout << "   - lmax = " << lmax << ", l_interpolate_max = " << l_interpolate_max << "\n";
	//}

	int Lmax = max_L12 + max_l3;
	int two_jmax_SixJ = 2 * lmax; // do we need to prestore 6j??
	
	// for angular integration in Gtilde
	double x_Gtilde  [Nx_Gtilde];
	double wx_Gtilde [Nx_Gtilde];

	calc_gauss_points (x_Gtilde, wx_Gtilde, -1.0, 1.0, Nx_Gtilde);

	if (print_content){std::cout << "   - Nalpha    = " <<  Jj_dim << std::endl;}
	if (print_content){std::cout << "   - Np_WP     = " <<  Np_WP << std::endl;}
	if (print_content){std::cout << "   - Nq_WP     = " <<  Nq_WP << std::endl;}
	if (print_content){std::cout << "   - Nphi      = " <<  Nphi << std::endl;}
	if (print_content){std::cout << "   - Nx        = " <<  Nx_Gtilde << std::endl;}
	if (print_content){std::cout << "   - lmax      = " <<  lmax << std::endl;}

	if (print_content){
		printf("   - Precalculating Wigner 6j-symbols \n");
	}
	int SixJ_size = (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1);
	if (SixJ_size < 0){
		raise_error("SixJ_array in make_permutation_matrix had negative size, likely an integer overflow. Check your dimensions.");
	}
	printf("     - Total prestore requirement is %zu doubles. Allocating arrays ... \n", SixJ_size);
	double *SixJ_array = new double [SixJ_size];
	printf("     - Success. Calculating ... \n");
	#pragma omp parallel for collapse(3)
	for (int two_l1 = 0; two_l1 <= two_jmax_SixJ; two_l1++){
		for (int two_l2 = 0; two_l2 <= two_jmax_SixJ; two_l2++){
			for (int two_l3 = 0; two_l3 <= two_jmax_SixJ; two_l3++){
				for (int two_l4 = 0; two_l4 <= two_jmax_SixJ; two_l4++){
					for (int two_l5 = 0; two_l5 <= two_jmax_SixJ; two_l5++){
						for (int two_l6 = 0; two_l6 <= two_jmax_SixJ; two_l6++){
							SixJ_array[
								  two_l1 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
								+ two_l2 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
								+ two_l3 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
								+ two_l4 * (two_jmax_SixJ + 1) * (two_jmax_SixJ + 1)
								+ two_l5 * (two_jmax_SixJ + 1)
								+ two_l6
							]
							// implement checks for quantum numbers because of bug in gsl library
								= gsl_sf_coupling_6j(two_l1, two_l2, two_l3, two_l4, two_l5, two_l6);
						}
					}
				}
			}
		}
	}
	printf("     - Done \n");

	if (print_content){
		printf("   - Precalculating Atilde \n");
	}
	MKL_INT64 Atilde_N = Jj_dim * Jj_dim * (Lmax + 1);
	printf("     - Total prestore requirement is %zu doubles. Allocating arrays ... \n", Atilde_N);
	double *Atilde_store = new double[Atilde_N];
	printf("     - Success. Calculating ... \n");
	for (MKL_INT64 i=0; i<Atilde_N; i++){
		Atilde_store[i] = 0.0;
	}
	for (int alpha = 0; alpha <= Jj_dim - 1; alpha++){
		int two_T_3N_alpha = two_T_3N_array[alpha];
		for (int alphaprime = 0; alphaprime <= Jj_dim - 1; alphaprime++){
			int two_T_3N_alphaprime = two_T_3N_array[alphaprime];
			if (two_T_3N_alpha==two_T_3N_alphaprime){
				int two_T = two_T_3N_alpha;
				for (int Ltotal = 0; Ltotal <= Lmax; Ltotal++){
					Atilde_store[alpha * Jj_dim * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal] = Atilde (alpha, alphaprime, Ltotal, Jj_dim, L12_Jj, l3_Jj, J12_Jj, two_j3_Jj, S12_Jj, T12_Jj, two_J, two_T, SixJ_array, two_jmax_SixJ);
				}
			}
