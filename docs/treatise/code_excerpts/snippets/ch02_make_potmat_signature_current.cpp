// ===============================================================
// 抽取自仓库 [current]: src/core/potential/make_potential_matrix.cpp
// 行号区段：50..110
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
/* Construct 2N potential matrices <k|v|k_p> for all 3N partial wave channels */
void calculate_potential_matrices_array_in_WP_basis(double*  V_WP_unco_array, int num_2N_unco_states,
													double*  V_WP_coup_array, int num_2N_coup_states,
													fwp_statespace 	 fwp_states,
													pw_3N_statespace pw_states,
													potential_model* pot_ptr,
													run_params run_parameters){
	
	/* Make local pointers & variables for pw-statespace */
	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* two_T_3N_array	= pw_states.two_T_3N_array;
	/* Make local pointers & variables for FWP-statespace */
	int 	Np_WP		 = fwp_states.Np_WP;
	double* p_WP_array	 = fwp_states.p_WP_array;
	int 	Np_per_WP	 = fwp_states.Np_per_WP;
	double* p_array		 = fwp_states.p_array;
	double* wp_array	 = fwp_states.wp_array;
	double* norm_p_array = fwp_states.norm_p_array;
	double* fp_array 	 = fwp_states.fp_array;

	/* This test will be reused several times */
	bool tensor_force_true = (run_parameters.tensor_force==true);

	int J_2N_max = run_parameters.J_2N_max;

	int Tz_nn = -1;
	int Tz_np =  0;

	/* Potential-model input arrays */
	//double V_WP_elements [6];	// Isoscalar wave-packet potential elements (WP)
	//double V_IS_elements [6];	// Isoscalar (IS)
	//double V_nn_elements [6];	// neutron-neutron (nn)
	//double V_np_elements [6];	// neutron-proton (np)

	/* Potential matrix indexing */
	//int idx_V_WP_uncoupled   = 0;
	//int idx_V_WP_upper_left  = 0;
	//int idx_V_WP_upper_right = 0;
	//int idx_V_WP_lower_left  = 0;
	//int idx_V_WP_lower_right = 0;

	/* Temporary track-keeping arrays for avoiding repeat matrix-calculations (initialize to false) */
	int unco_array_size = num_2N_unco_states;
	int coup_array_size = num_2N_coup_states;
	
	bool* matrix_calculated_unco_array = new bool [unco_array_size];
	bool* matrix_calculated_coup_array = new bool [coup_array_size];
	for (int i=0; i<unco_array_size; i++){
		matrix_calculated_unco_array[i] = false;
	}
	for (int i=0; i<coup_array_size; i++){
		matrix_calculated_coup_array[i] = false;
	}

	printf("   - There are %d uncoupled 2N-channels and %d coupled 2N-channels \n", unco_array_size, coup_array_size);

	/* Row state */
