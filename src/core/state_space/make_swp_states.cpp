#include "make_swp_states.h"
#include <vector>

namespace {

// [EN] The deuteron is the unique bound subsystem carried by the triplet S-wave channel. / [CN] 氘核是三重态
// S 波通道携带的唯一束缚子系统。
bool is_triplet_s_wave_channel(int L_2N, int S_2N, int J_2N, int T_2N){
	return L_2N==0 && S_2N==1 && J_2N==1 && T_2N==0;
}

// [EN] 1S0 needs special bookkeeping because optional isospin breaking reuses the coupled-channel storage pattern.
// / [CN] 1S0 需要特殊 bookkeeping，因为可选的同位旋破缺会复用耦合通道的存储模式。
bool is_singlet_s_wave_channel(int L_2N, int S_2N, int J_2N){
	return L_2N==0 && S_2N==0 && J_2N==0;
}

// [EN] Decide whether a two-body block occupies the uncoupled Np x Np storage or the coupled 2Np x 2Np storage.
// / [CN] 判断某个两体块应放入非耦合的 Np x Np 存储，还是耦合的 2Np x 2Np 存储。
bool uses_coupled_storage(int L_2N_row,
						  int L_2N_col,
						  int S_2N,
						  int J_2N,
						  const run_params& run_parameters){
	const bool coupled_via_L_2N = run_parameters.tensor_force && (L_2N_row!=L_2N_col || (L_2N_row==L_2N_col && L_2N_row!=J_2N && J_2N!=0));
	const bool coupled_via_T_3N = is_singlet_s_wave_channel(L_2N_row, S_2N, J_2N) && run_parameters.isospin_breaking_1S0;
	if (coupled_via_L_2N && coupled_via_T_3N){
		raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
	}
	return coupled_via_L_2N || coupled_via_T_3N;
}

// [EN] Claim one distinct 2N Hamiltonian block the first time it is encountered in the enclosing 3N basis loop.
// / [CN] 在外层三体基循环里首次遇到某个不同的 2N Hamiltonian 块时，将其标记为“已领取”。
bool claim_channel(std::vector<bool>& channel_done_flags, int channel_index){
	if (channel_done_flags[channel_index]){
		return false;
	}
	channel_done_flags[channel_index] = true;
	return true;
}

// [EN] The free packet Hamiltonian is identical for every uncoupled branch and for each leg of a coupled branch, so
// we build the repeated H0 tables once here. / [CN] 自由波包 Hamiltonian 对所有非耦合分支以及每条耦合分支腿都是相同的，
// 因此在这里一次性构建并复用这些重复的 H0 表。
void fill_free_hamiltonian_branches(std::vector<double>& free_hamiltonian_branches,
									int num_branches,
									int Np_WP,
									double* p_WP_array){
	for (int idx_branch=0; idx_branch<num_branches; idx_branch++){
		construct_free_hamiltonian(&free_hamiltonian_branches[idx_branch * Np_WP],
								   Np_WP,
								   p_WP_array);
	}
}

} // namespace

/* Finds the eigenvalues and eigenvectors
 * of a real, symmetric matrix A.
 * For simplicity here, A must be
 * stored as an upper triangle. w will
 * be filled with the eigenvalues in
 * ascending order, and z will be a matrix
 * with the corresponding eigenvectors
 * (column by column). Lastly, we work
 * with row major matrices, as is usual
 * with C & C++ */
 
void diagonalize_real_symm_matrix(float *A, float *w, float *z, int N){
	char jobz = 'V';
	char uplo = 'U';
	
	LAPACKE_sspevd(LAPACK_ROW_MAJOR, jobz, uplo, N, A, w, z, N);
}

void diagonalize_real_symm_matrix(double *A, double *w, double *z, int N){
	char jobz = 'V';
	char uplo = 'U';
	
	LAPACKE_dspevd(LAPACK_ROW_MAJOR, jobz, uplo, N, A, w, z, N);
}

void construct_free_hamiltonian(double* H0_WP_array,
							    int Np_WP, double* p_WP_array){
	
    /* Loop over p-momenta */
    for (int idx_p=0; idx_p<Np_WP; idx_p++){
		/* Kinetic energy at bin-boundaries */
		double p1 = p_WP_array[idx_p];
		double p2 = p_WP_array[idx_p+1];
	
		/* Free Hamiltonian for 2-nucleon pair for momentum WPs */
		H0_WP_array[idx_p] = (p2*p2 + p2*p1 + p1*p1)/(6*mu23);
		/* Free Hamiltonian for 2-nucleon pair for energy WPs */
        //H0_WP_array[idx_p] = (p2*p2 + p1*p1)/(2*MN);
	}
}

/* Constructs a NN-pair full Hamiltonians as upper-triangular from given arrays */
void construct_full_hamiltonian(double* mat_ptr_H,
								double* mat_ptr_H0,
								double* mat_ptr_V,
                                int     mat_dim){

    /* Row p-momentum index loop */
    for (int idx_bin_r=0; idx_bin_r<mat_dim; idx_bin_r++){
					
		/* Because we only store the upper triangular part we need special indexing for H */
		int idx_r_r = mat_dim*(mat_dim-1)/2 - (mat_dim-idx_bin_r)*(mat_dim-idx_bin_r-1)/2 + idx_bin_r;

		/* Set diagonal to free kinetic energy given by H0_array */
		mat_ptr_H[idx_r_r] += mat_ptr_H0[idx_bin_r];

        /* Column p-momentum index loop (note it starts on idx_bin_r since H is upper triangular) */
        for (int idx_bin_c=idx_bin_r; idx_bin_c<mat_dim; idx_bin_c++){
						
			/* Because we only store the upper triangular part we need special indexing for H */
			int idx_r_c = mat_dim*(mat_dim-1)/2 - (mat_dim-idx_bin_r)*(mat_dim-idx_bin_r-1)/2 + idx_bin_c;

			/* Add potential element to Hamiltonian */
			mat_ptr_H[idx_r_c] += mat_ptr_V[idx_bin_r*mat_dim + idx_bin_c];
		}
	}

	//for (int idx_bin_c=0; idx_bin_c<mat_dim; idx_bin_c++){
	//	printf(" %.8e \n", mat_ptr_H0[idx_bin_c]);
	//}
	//printf("\n");
	//for (int idx_bin_r=0; idx_bin_r<mat_dim; idx_bin_r++){
    //    for (int idx_bin_c=0; idx_bin_c<mat_dim; idx_bin_c++){
	//		printf(" %.8e", mat_ptr_V[idx_bin_r*mat_dim + idx_bin_c]);
	//	}
	//	printf("\n");
	//}
	//printf("\n");
}

// [EN] Coupled-channel diagonalization returns the two eigenbranches interleaved. We reorder them into contiguous
// blocks so the later SWP storage layout matches the solver's expectation that each branch can be addressed by a
// simple offset, which keeps the resolvent and on-shell indexing logic branch-local. / [CN] 耦合通道对角化后，两条本征分支通常是交错返回的；这里把它们重排成连续块，使后续 SWP 存储可以通过简单偏移访问各分支，并让分辨算符和 on-shell 索引逻辑保持分支局部化。
void reorder_coupled_eigenspectrum(double* eigenvalues,
								   double* eigenvectors,
								   int     Np_WP){
	
	/* Temporary arrays to hold values as we reorder indices */
	std::vector<double> buffer_array_vals(2*Np_WP);
	std::vector<double> buffer_array_vecs(4*Np_WP*Np_WP);

	for (int i=0; i<Np_WP; i++){
		/* Temporarily store reorder of eigenvalues */
		buffer_array_vals[i] 	   = eigenvalues[2*i];
		buffer_array_vals[i+Np_WP] = eigenvalues[2*i+1];

		for (int j=0; j<2*Np_WP; j++){
			/* Temporarily store reorder of eigenvectors */
			buffer_array_vecs[j*2*Np_WP + i] 	     = eigenvectors[j*2*Np_WP + 2*i];
			buffer_array_vecs[j*2*Np_WP + i + Np_WP] = eigenvectors[j*2*Np_WP + 2*i+1];
		}
	}

	/* Write reordered elements back into original arrays */
	for (int i=0; i<2*Np_WP; i++){
		eigenvalues[i] = buffer_array_vals[i];
		for (int j=0; j<2*Np_WP; j++){
			eigenvectors[i*2*Np_WP + j] = buffer_array_vecs[i*2*Np_WP + j];
		}
	}
}

void look_for_unphysical_bound_states(double* eigenvalues,
						   			  int     mat_dim,
						   			  bool    chn_3S1,
									  double& E_bound){
											 
	int num_bound_states_found = 0;

	// [EN] In the packet discretization a physical bound state still appears as a negative eigenvalue of H=H0+V.
	// The rest of the code assumes exactly one such state in 3S1 and none elsewhere. / [CN] 在波包离散化里，物理束缚态
	// 仍然表现为 H=H0+V 的负本征值；后续代码假定只有 3S1 中恰好存在一个，其他通道则一个都没有。
	/* Count number of bound states in eigenspectrum */
	for (int idx=0; idx<mat_dim; idx++){
		if (eigenvalues[idx]<0){
			num_bound_states_found += 1;
			printf("   - Found bound state with energy %.10f MeV \n", eigenvalues[idx]);

			/* Save bound-state energy to special variable. This is okay since the code halts
			 * if more than one bound state is located */
			E_bound = eigenvalues[idx];
		}
	}
	
	/* See if we find the expected number of bounds states (1 for 3S1, 0 otherwise),
	 * if not then we abort program (unphysical scenario) */
	
	if (num_bound_states_found==0 and chn_3S1==false);
	else if (num_bound_states_found==1 and chn_3S1==true);
	else if (num_bound_states_found!=0){
		raise_error("Found unphysical bound states in NN-pair Hamiltonian eigenspectrum!");
	}
	else{
		raise_error("Didn't find any bound states in NN-pair Hamiltonian eigenspectrum!");
	}
}

void make_swp_bin_boundaries(double* eigenvalues,
							 double* e_SWP_array,
							 int	 Np_WP,
							 bool    coupled,
							 bool    chn_3S1){

	// [EN] Neighboring eigenvalues define the edges of one interacting packet cell. For the deuteron channel the
	// first cell is split so the bound state sits on its own branch below threshold. / [CN] 相邻本征值共同定义一个
	// 相互作用波包单元的边界；对氘核通道而言，第一个单元会被特殊处理，使束缚态单独占据阈值以下的一条分支。
	/* Set first boundary */
	e_SWP_array[0] = 0;
	if (coupled){
		e_SWP_array[Np_WP+1] = 0;
	}
	
	/* Set all the middle boundaries */
	for (int idx_p=0; idx_p<Np_WP-1; idx_p++){
		e_SWP_array[idx_p+1] = 0.5*(eigenvalues[idx_p+1] + eigenvalues[idx_p]);

		if (coupled){
			e_SWP_array[idx_p+1+Np_WP+1] = 0.5*(eigenvalues[idx_p+1+Np_WP] + eigenvalues[idx_p+Np_WP]);
		}
	}

	/* Set end boundary */
	e_SWP_array[Np_WP] = eigenvalues[Np_WP-1] + 0.5*(e_SWP_array[Np_WP-1] - e_SWP_array[Np_WP-2]);
	if (coupled){
		e_SWP_array[2*Np_WP+1] = eigenvalues[2*Np_WP-1] + 0.5*(e_SWP_array[2*Np_WP] - e_SWP_array[2*Np_WP-1]);
	}

	/* Modify first boundary for 3S1 */
	if (chn_3S1){
		e_SWP_array[0] = eigenvalues[0];
		e_SWP_array[1] = 0;
	}

	//for (int idx_p=0; idx_p<Np_WP+1; idx_p++){
	//	std::cout << e_SWP_array[idx_p] << std::endl;
	//}
	//std::cout << std::endl;
}

void store_swp_kinematics(swp_statespace swp_states,
						  run_params run_parameters){
	
	/* Local pointer */
	double* q_WP_array = swp_states.q_WP_array;

	// [EN] These tables expose the discrete spectator kinematics seen by the solver so downstream scripts can map
	// packet indices back to the physical E_cm and T_lab scales. / [CN] 这些表把求解器实际使用的离散 spectator 运动学
	// 导出出来，使后处理脚本能够把波包索引重新映射回物理的 E_cm 与 T_lab 标度。
	std::vector<double> Eq_WP_boundaries(swp_states.Nq_WP+1);
	std::vector<double> Tlab_WP_boundaries(swp_states.Nq_WP+1);
	for (size_t q_WP_idx=0; q_WP_idx<swp_states.Nq_WP+1; q_WP_idx++){
		Eq_WP_boundaries[q_WP_idx]   = com_q_momentum_to_com_energy(q_WP_array[q_WP_idx]);
		Tlab_WP_boundaries[q_WP_idx] = com_momentum_to_lab_energy(q_WP_array[q_WP_idx], swp_states.E_bound);
	}
	std::vector<double> q_WP_midpoints(swp_states.Nq_WP);
	std::vector<double> Eq_WP_midpoints(swp_states.Nq_WP);
	std::vector<double> Tlab_WP_midpoints(swp_states.Nq_WP);
	for (size_t q_WP_idx=0; q_WP_idx<swp_states.Nq_WP; q_WP_idx++){
		double Eq_lower = 0.5*(q_WP_array[q_WP_idx]   * q_WP_array[q_WP_idx])  /mu1(swp_states.E_bound);
		double Eq_upper = 0.5*(q_WP_array[q_WP_idx+1] * q_WP_array[q_WP_idx+1])/mu1(swp_states.E_bound);
		double E_com = 0.5*(Eq_upper + Eq_lower);
		q_WP_midpoints[q_WP_idx]    = com_energy_to_com_q_momentum(E_com);
		Eq_WP_midpoints[q_WP_idx]   = E_com;
		Tlab_WP_midpoints[q_WP_idx] = com_momentum_to_lab_energy(q_WP_midpoints[q_WP_idx], swp_states.E_bound);
	}
	
	std::string q_kinematics_filename = run_parameters.output_folder + "/" + "q_kinematics_Nq_" + std::to_string(swp_states.Nq_WP) + ".txt";
	store_q_WP_kinematics_txt(swp_states.Nq_WP,
							  swp_states.q_WP_array,
							  Eq_WP_boundaries.data(),
							  Tlab_WP_boundaries.data(),
							  q_WP_midpoints.data(),
							  Eq_WP_midpoints.data(),
							  Tlab_WP_midpoints.data(),
							  q_kinematics_filename);
}


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

	/* TEMP */
	swp_states.e_SWP_unco_array		= e_SWP_unco_array;
	swp_states.e_SWP_coup_array		= e_SWP_coup_array;
	swp_states.C_SWP_unco_array		= C_WP_unco_array;
	swp_states.C_SWP_coup_array		= C_WP_coup_array;
	swp_states.num_2N_unco_states	= num_2N_unco_states;
	swp_states.num_2N_coup_states	= num_2N_coup_states;
	swp_states.E_bound				= E_bound;

	printf("   - Storing kinematic values of WP statespace to txt-file ... \n");
	store_swp_kinematics(swp_states, run_parameters);
	printf("     - Done \n");
}
