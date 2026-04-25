// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_swp_states.cpp
// 行号区段：287..480
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void make_swp_states(double* e_SWP_unco_array,
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

	// [EN] This is the WPCD basis-change step: for each two-body channel we diagonalize H=H0+V in the free packet
	// basis, interpret negative eigenvalues in the 3S1 sector as the deuteron bound state, and reuse the resulting
	// eigenvectors as the C matrices that map fWP states into interacting SWP states. / [CN] 这里执行的是 WPCD 的基变换步骤：对每个两体通道在自由波包基上对角化 H=H0+V，把 3S1 扇区中的负本征值解释为氘核束缚态，并把得到的本征向量作为 C 矩阵，把 fWP 态映射到相互作用的 SWP 态。

	/* Make local pointers & variables for pw-statespace */
	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
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
	
	/* Number of uncoupled and coupled 2N-channels */
	int num_unco_chns = num_2N_unco_states;
	int num_coup_chns = num_2N_coup_states;

	/* Check-lists to keep track of which 2N Hamiltonian diagonalizations have been done.
	 * This removes excessive work due to distinct 3N channels containing equal
	 * 2N channels, as well as overwriting existing calculations (thus giving wrongful results) */
	std::vector<bool> check_list_unco(num_unco_chns, false);
	std::vector<bool> check_list_coup(num_coup_chns, false);

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
	std::vector<double> H_WP_unco_array(H_unco_array_size, 0.0);
	std::vector<double> H_WP_coup_array(H_coup_array_size, 0.0);

	/* Free Hamiltonian arrays */
	int H0_unco_array_size =   Np_WP;
    int H0_coup_array_size = 2*Np_WP;
	std::vector<double> H0_WP_unco_array(H0_unco_array_size, 0.0);
	std::vector<double> H0_WP_coup_array(H0_coup_array_size, 0.0);

	/* Fill free Hamiltonian arrays. All 3 calls will produce the same,
	 * but it's quite fast and makes the code clearer to interpret */
	fill_free_hamiltonian_branches(H0_WP_unco_array, 1, Np_WP, p_WP_array);
	fill_free_hamiltonian_branches(H0_WP_coup_array, 2, Np_WP, p_WP_array);
	
	// [EN] Each unique 2N partial wave defines one Hamiltonian block. Several 3N basis states can reference the same
	// block, so we diagonalize each block once and reuse the resulting SWP spectrum and C matrix wherever it appears.
	// / [CN] 每个唯一的两体分波都定义了一个 Hamiltonian 块。多个三体基态可能引用同一个块，因此这里每个块只对角化一次，
	// 然后在所有出现它的地方复用对应的 SWP 谱和 C 矩阵。
	/* Row state */
	printf("   - Diagonalizing 2N Hamiltonians and constructing SWP boundaries ... \n");
    for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
        int L_r = L_2N_array[idx_alpha_r];
        int S_r = S_2N_array[idx_alpha_r];
        int J_r = J_2N_array[idx_alpha_r];
        int T_r = T_2N_array[idx_alpha_r];

        /* Column state */
        for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
            int L_c = L_2N_array[idx_alpha_c];
            int S_c = S_2N_array[idx_alpha_c];
            int J_c = J_2N_array[idx_alpha_c];
            int T_c = T_2N_array[idx_alpha_c];

            /* Check if possible channel through interaction */
            if (T_r==T_c and J_r==J_c and S_r==S_c and abs(L_r-L_c)<=2){
				
                /* Detemine if this is a coupled channel.
				 * !!! With isospin symmetry-breaking we count 1S0 as a coupled matrix via T_3N-coupling !!! */
				const bool coupled_matrix = uses_coupled_storage(L_r,
																	 L_c,
																	 S_r,
																	 J_r,
																	 run_parameters);

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
					if (claim_channel(check_list_coup, chn_idx)==false){
						continue;
					}

					/* H-matrices are stored as upper-triangular -> special indexing and step-length */
					mat_ptr_H 		= &H_WP_coup_array [chn_idx * mat_dim*(mat_dim+1)/2];
					mat_ptr_V 		= &V_WP_coup_array [chn_idx * mat_dim*mat_dim];
					mat_ptr_C 		= &C_WP_coup_array [chn_idx * mat_dim*mat_dim];
					e_SWP_array_ptr = &e_SWP_coup_array[chn_idx * 2*(Np_WP+1)];
					mat_ptr_H0 		= H0_WP_coup_array.data();
                }
			    else{
					mat_dim = Np_WP;
					chn_idx = unique_2N_idx(L_r, S_r, J_r, T_r, coupled_matrix, run_parameters);

					/* Check if 2N channels diagonalization has already
					 * been performed in previous loop-iterations,
					 * and if not then set to true */
					if (claim_channel(check_list_unco, chn_idx)==false){
						continue;
					}

					/* H-matrices are stored as upper-triangular -> special indexing and step-length */
					mat_ptr_H 		= &H_WP_unco_array [chn_idx * mat_dim*(mat_dim+1)/2];
					mat_ptr_V 		= &V_WP_unco_array [chn_idx * mat_dim*mat_dim];
					mat_ptr_C 		= &C_WP_unco_array [chn_idx * mat_dim*mat_dim];
					e_SWP_array_ptr = &e_SWP_unco_array[chn_idx * (Np_WP+1)];
					mat_ptr_H0 		= H0_WP_unco_array.data();
                }

				/* Construct channel Hamiltonian */
				construct_full_hamiltonian(mat_ptr_H,
										   mat_ptr_H0,
										   mat_ptr_V,
                                		   mat_dim);

				/* Hamiltonian eigenvalue array */
				std::vector<double> eigenvalues(mat_dim);

				/* Diagonalize channel Hamiltonian - fill eigenvalues and C_array coefficients */
				diagonalize_real_symm_matrix(mat_ptr_H, eigenvalues.data(), mat_ptr_C, mat_dim);

				/* Abort if unphysical bound states are found in eigenvalues,
				   or if 3S1-bound state is missing */
				const bool chn_3S1 = is_triplet_s_wave_channel(L_r, S_r, J_r, T_r);
				look_for_unphysical_bound_states(eigenvalues.data(), mat_dim, chn_3S1, E_bound);
				
				/* The eigenspectrum of coupled channels is returned in ascending
				 * order from the diagonalization routine, so we reorder the spectrum here */
				if (coupled_matrix){
					reorder_coupled_eigenspectrum(eigenvalues.data(),
												  mat_ptr_C,
												  Np_WP);
				}
				
				// [EN] Neighboring eigenvalues define the packet cell boundaries in the interacting basis. This is the
				// discrete counterpart of replacing the continuum spectrum by averaged scattering packets in the notes.
				// / [CN] 相邻本征值共同定义相互作用基中的波包单元边界；这正是讲稿里“用平均散射波包替代连续谱”的离散对应。
				make_swp_bin_boundaries(eigenvalues.data(),
										e_SWP_array_ptr,
										Np_WP,
										coupled_matrix,
										chn_3S1);
			}
		}
	}
	printf("     - Done \n");
