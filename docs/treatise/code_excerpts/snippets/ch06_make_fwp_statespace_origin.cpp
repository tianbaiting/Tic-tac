// ===============================================================
// 抽取自仓库 [origin]: CPP/make_wp_states.cpp
// 行号区段：116..211
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void make_fwp_statespace(fwp_statespace& fwp_states, run_params run_parameters){
	printf("Constructing wave-packet (WP) state space ... \n");

	/* Copy space dimensions from input struct */
	fwp_states.Np_WP 	 = run_parameters.Np_WP;
	fwp_states.Nq_WP 	 = run_parameters.Nq_WP;
	fwp_states.Np_per_WP = run_parameters.Np_per_WP;
	fwp_states.Nq_per_WP = run_parameters.Nq_per_WP;

	/* Temporary dimensional integer for p_array, wp_array, fp_array (and same for q_...) */
	int Np = 0;
	int Nq = 0;
	if (run_parameters.midpoint_approx==true){
		Np = fwp_states.Np_WP;
		Nq = fwp_states.Nq_WP;
	}
	else{
		Np = fwp_states.Np_per_WP * fwp_states.Np_WP;
		Nq = fwp_states.Nq_per_WP * fwp_states.Nq_WP;
	}

	/* Make bin boundaries */
	printf(" - Constructing wave-packet (WP) p-momentum bin boundaries ... \n");
	fwp_states.p_WP_array = new double [fwp_states.Np_WP+1];
	make_p_bin_grid(fwp_states, run_parameters);
	printf("   - Done \n");
	printf(" - Constructing wave-packet (WP) q-momentum bin boundaries ... \n");
	fwp_states.q_WP_array = new double [fwp_states.Nq_WP+1];
	make_q_bin_grid(fwp_states, run_parameters);
	printf("   - Done \n");

	/* Make mesh for calculating operators in fWP basis.
	 * if midpoint_approx==true: Make Gauss-Legendre quadrature meshes inside each bin
	   else:					 Extract bin-midpoints from bin boundaries */
	if (run_parameters.midpoint_approx==true){
		printf(" - Constructing p bin-midpoint mesh ... \n");
		fwp_states.p_array  = new double [Np];
		for (int idx_p=0; idx_p<Np; idx_p++){
			fwp_states.p_array[idx_p] = 0.5*(fwp_states.p_WP_array[idx_p] + fwp_states.p_WP_array[idx_p+1]);
		}
		printf("   - Done \n");
		printf(" - Constructing q bin-midpoint mesh ... \n");
		fwp_states.q_array  = new double [Nq];
		for (int idx_q=0; idx_q<Nq; idx_q++){
			fwp_states.q_array[idx_q] = 0.5*(fwp_states.q_WP_array[idx_q] + fwp_states.q_WP_array[idx_q+1]);
		}
		printf("   - Done \n");
	}
	else{
		printf(" - Constructing p quadrature mesh per WP, for all WPs ... \n");
		fwp_states.p_array  = new double [Np];
		fwp_states.wp_array = new double [Np];
		make_p_bin_quadrature_grids(fwp_states);
		printf("   - Done \n");
		printf(" - Constructing q quadrature mesh per WP, for all WPs ... \n");
		fwp_states.q_array  = new double [Nq];
		fwp_states.wq_array = new double [Nq];
		make_q_bin_quadrature_grids(fwp_states);
		printf("   - Done \n");
	}

	/* Calculate weight-functions */
	printf(" - Calculating weight-functions for p-momentum WPs ... \n");
	fwp_states.fp_array  = new double [Np];
	for (int idx_p=0; idx_p<Np; idx_p++){
		fwp_states.fp_array[idx_p] = p_weight_function(fwp_states.p_array[idx_p]);
	}
	printf("   - Done \n");
	printf(" - Calculating weight-functions for q-momentum WPs ... \n");
	fwp_states.fq_array  = new double [Nq];
	for (int idx_q=0; idx_q<Nq; idx_q++){
		fwp_states.fq_array[idx_q] = p_weight_function(fwp_states.q_array[idx_q]);
	}
	printf("   - Done \n");

	/* Calculate normalization constants */
	printf(" - Normalizing p-momentum WPs ... \n");
	fwp_states.norm_p_array  = new double [fwp_states.Np_WP];
	for (int idx_p=0; idx_p<fwp_states.Np_WP; idx_p++){
		fwp_states.norm_p_array[idx_p] = p_normalization(fwp_states.p_WP_array[idx_p], fwp_states.p_WP_array[idx_p+1]);
	}
	printf("   - Done \n");
	printf(" - Normalizing q-momentum WPs ... \n");
	fwp_states.norm_q_array  = new double [fwp_states.Nq_WP];
	for (int idx_q=0; idx_q<fwp_states.Nq_WP; idx_q++){
		fwp_states.norm_q_array[idx_q] = q_normalization(fwp_states.q_WP_array[idx_q], fwp_states.q_WP_array[idx_q+1]);
	}
	printf("   - Done \n");

	/* Store boundaries for post-processing of output  */
	printf(" - Storing q boundaries to CSV-file ... \n");
	std::string q_boundaries_filename = run_parameters.output_folder + "/" + "q_boundaries_Nq_" + std::to_string(fwp_states.Nq_WP) + ".csv";
	store_q_WP_boundaries_csv(fwp_states, q_boundaries_filename);
	printf("   - Done \n");

	printf(" - Done \n");
