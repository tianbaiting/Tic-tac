// ===============================================================
// 抽取自仓库 [origin]: CPP/make_swp_states.cpp
// 行号区段：222..449
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
					 double* e_SWP_coup_array,
					 double* C_WP_unco_array,
					 double* C_WP_coup_array,
					 double* V_WP_unco_array,
                     double* V_WP_coup_array,
					 int num_2N_unco_states,
					 int num_2N_coup_states,
					 double& E_bound,
					 fwp_statespace  fwp_states,
					 swp_statespace& swp_states,
					 pw_3N_statespace pw_states,
					 run_params run_parameters){

	/* Make local pointers & variables for pw-statespace */
	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* two_T_3N_array	= pw_states.two_T_3N_array;
	/* Make local pointers & variables for WP-statespace */
	int 	Np_WP		= fwp_states.Np_WP;
	double* p_WP_array	= fwp_states.p_WP_array;
	/* Make local pointers & variables for SWP-statespace */
	//double* e_SWP_unco_array   = swp_states.e_SWP_unco_array;
	//double* e_SWP_coup_array   = swp_states.e_SWP_coup_array;
	//double* C_WP_unco_array	   = swp_states.C_WP_unco_array;
	//double* C_WP_coup_array	   = swp_states.C_WP_coup_array;
	//int 	num_2N_unco_states = swp_states.num_2N_unco_states;
	//int 	num_2N_coup_states = swp_states.num_2N_coup_states;
	//double& E_bound			   = swp_states.E_bound;

	/* Copy pointers and variables from fwp-statespace */
	swp_states.Np_WP 	  = fwp_states.Np_WP;
	swp_states.Nq_WP 	  = fwp_states.Nq_WP;
	swp_states.q_WP_array = fwp_states.q_WP_array;
	
	/* This test will be reused several times */
	bool tensor_force_true = (run_parameters.tensor_force==true);

	/* Number of uncoupled and coupled 2N-channels */
	int num_unco_chns = num_2N_unco_states;
	int num_coup_chns = num_2N_coup_states;

	/* Check-lists to keep track of which 2N Hamiltonian diagonalizations have been done.
	 * This removes excessive work due to distinct 3N channels containing equal
	 * 2N channels, as well as overwriting existing calculations (thus giving wrongful results) */
	bool* check_list_unco = new bool [num_unco_chns];
	bool* check_list_coup = new bool [num_coup_chns];

	/* Set check_list-arrays to false */
	for (int i=0; i<num_unco_chns; i++){
		check_list_unco[i] = false;
	}
	for (int i=0; i<num_coup_chns; i++){
		check_list_coup[i] = false;
	}

	/* Boundaries of scattering wave-packets (SWP) in energy */
	int e_SWP_unco_array_size =   (Np_WP+1) * num_unco_chns;
	int e_SWP_coup_array_size = 2*(Np_WP+1) * num_coup_chns;

	/* Set e_SWP-arrays to zero */
	for (int i=0; i<e_SWP_unco_array_size; i++){
		e_SWP_unco_array[i] = 0;
	}
	for (int i=0; i<e_SWP_coup_array_size; i++){
		e_SWP_coup_array[i] = 0;
	}

	/* The Hamiltonian matrices are symmetric and stored as upper-triangular
	 * for computational efficiency when diagonalizing */
	int H_unco_array_size =   Np_WP*(  Np_WP+1)/2 * num_unco_chns;
    int H_coup_array_size = 2*Np_WP*(2*Np_WP+1)/2 * num_coup_chns;
	double* H_WP_unco_array = new double [H_unco_array_size];
	double* H_WP_coup_array = new double [H_coup_array_size];

	/* Set H_WP-arrays to zero */
	for (int i=0; i<H_unco_array_size; i++){
		H_WP_unco_array[i] = 0;
	}
	for (int i=0; i<H_coup_array_size; i++){
		H_WP_coup_array[i] = 0;
	}

	/* Free Hamiltonian arrays */
	int H0_unco_array_size =   Np_WP;
    int H0_coup_array_size = 2*Np_WP;
	double* H0_WP_unco_array = new double [H0_unco_array_size];
	double* H0_WP_coup_array = new double [H0_coup_array_size];

	/* Set H0_WP-arrays to zero */
	for (int i=0; i<H0_unco_array_size; i++){
		H0_WP_unco_array[i] = 0;
	}
	for (int i=0; i<H0_coup_array_size; i++){
		H0_WP_coup_array[i] = 0;
	}

	/* Fill free Hamiltonian arrays. All 3 calls will produce the same,
	 * but it's quite fast and makes the code clearer to interpret */
	construct_free_hamiltonian( H0_WP_unco_array,        Np_WP, p_WP_array);
	construct_free_hamiltonian(&H0_WP_coup_array[0],     Np_WP, p_WP_array);
	construct_free_hamiltonian(&H0_WP_coup_array[Np_WP], Np_WP, p_WP_array);
	
	/* Row state */
	printf("   - Diagonalizing 2N Hamiltonians and constructing SWP boundaries ... \n");
    for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
        int L_r = L_2N_array[idx_alpha_r];
        int S_r = S_2N_array[idx_alpha_r];
        int J_r = J_2N_array[idx_alpha_r];
        int T_r = T_2N_array[idx_alpha_r];

		int two_T_3N_r = two_T_3N_array[idx_alpha_r];
        
        /* Column state */
        for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
            int L_c = L_2N_array[idx_alpha_c];
            int S_c = S_2N_array[idx_alpha_c];
            int J_c = J_2N_array[idx_alpha_c];
            int T_c = T_2N_array[idx_alpha_c];

			int two_T_3N_c = two_T_3N_array[idx_alpha_c];

            /* Check if possible channel through interaction */
            if (T_r==T_c and J_r==J_c and S_r==S_c and abs(L_r-L_c)<=2){
				
                /* Detemine if this is a coupled channel.
				 * !!! With isospin symmetry-breaking we count 1S0 as a coupled matrix via T_3N-coupling !!! */
				bool coupled_matrix = false;
				bool state_1S0 = (S_r==0 && J_r==0 && L_r==0);
				bool coupled_via_L_2N = (tensor_force_true && (L_r!=L_c || (L_r==L_c && L_r!=J_r && J_r!=0)));
				bool coupled_via_T_3N = (state_1S0==true && run_parameters.isospin_breaking_1S0==true);
				if (coupled_via_L_2N && coupled_via_T_3N){
					raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
				}
				if (coupled_via_L_2N || coupled_via_T_3N){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}

				/* Hamiltonian matrix pointer and dimension
                 * Indexing format of Hamitonian arrays: (channel index)*(num rows)*(num columns) + (row index)*(row length) + (column index) */
				int     mat_dim   		= 0;
				int		chn_idx	  		= 0;
				double* mat_ptr_H 		= NULL;
				double* mat_ptr_H0 		= NULL;
				double* mat_ptr_V 		= NULL;
				double* mat_ptr_C 		= NULL;
				double* e_SWP_array_ptr = NULL;
                if (coupled_matrix){
					mat_dim = 2*Np_WP;
					chn_idx = unique_2N_idx(L_r, S_r, J_r, T_r, coupled_matrix, run_parameters);

					/* Check if 2N channels diagonalization has already
					 * been performed in previous loop-iterations,
					 * and if not then set to true */
					if (check_list_coup[chn_idx]==true){
						continue;
					}
					else{
						check_list_coup[chn_idx] = true;
					}

					/* H-matrices are stored as upper-triangular -> special indexing and step-length */
					mat_ptr_H 		= &H_WP_coup_array [chn_idx * mat_dim*(mat_dim+1)/2];
					mat_ptr_V 		= &V_WP_coup_array [chn_idx * mat_dim*mat_dim];
					mat_ptr_C 		= &C_WP_coup_array [chn_idx * mat_dim*mat_dim];
					e_SWP_array_ptr = &e_SWP_coup_array[chn_idx * 2*(Np_WP+1)];
					mat_ptr_H0 		= H0_WP_coup_array;
                }
			    else{
					mat_dim = Np_WP;
					chn_idx = unique_2N_idx(L_r, S_r, J_r, T_r, coupled_matrix, run_parameters);

					/* Check if 2N channels diagonalization has already
					 * been performed in previous loop-iterations,
					 * and if not then set to true */
					if (check_list_unco[chn_idx]==true){
						continue;
					}
					else{
						check_list_unco[chn_idx] = true;
					}

					/* H-matrices are stored as upper-triangular -> special indexing and step-length */
					mat_ptr_H 		= &H_WP_unco_array [chn_idx * mat_dim*(mat_dim+1)/2];
					mat_ptr_V 		= &V_WP_unco_array [chn_idx * mat_dim*mat_dim];
					mat_ptr_C 		= &C_WP_unco_array [chn_idx * mat_dim*mat_dim];
					e_SWP_array_ptr = &e_SWP_unco_array[chn_idx * (Np_WP+1)];
					mat_ptr_H0 		= H0_WP_unco_array;
                }

				/* Construct channel Hamiltonian */
				construct_full_hamiltonian(mat_ptr_H,
										   mat_ptr_H0,
										   mat_ptr_V,
                                		   mat_dim);

				/* Hamiltonian eigenvalue array */
				double eigenvalues [mat_dim];

				/* Diagonalize channel Hamiltonian - fill eigenvalues and C_array coefficients */
				diagonalize_real_symm_matrix(mat_ptr_H, eigenvalues, mat_ptr_C, mat_dim);

				/* Abort if unphysical bound states are found in eigenvalues,
				   or if 3S1-bound state is missing */
				bool chn_3S1 = (L_r==0 && S_r==1 && J_r==1 && T_r==0);
				look_for_unphysical_bound_states(eigenvalues, mat_dim, chn_3S1, E_bound);
				
				/* The eigenspectrum of coupled channels is returned in ascending
				 * order from the diagonalization routine, so we reorder the spectrum here */
				if (coupled_matrix){
					reorder_coupled_eigenspectrum(eigenvalues,
												  mat_ptr_C,
												  Np_WP);
				}
				
				/* Construct energy bin boundaries for swp states */
				make_swp_bin_boundaries(eigenvalues,
										e_SWP_array_ptr,
										Np_WP,
										coupled_matrix,
										chn_3S1);
			}
		}
	}
	printf("     - Done \n");

