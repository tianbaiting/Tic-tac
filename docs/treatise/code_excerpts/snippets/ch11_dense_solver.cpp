// ===============================================================
// 抽取自仓库 [current]: src/core/faddeev_solver/solve_faddeev.cpp
// 行号区段：676..776
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
/* Solves the Faddeev equations
 * U = P*V + P*V*G*U
 * on the form L*U = R, where L and R are the left-
 * and right-handed sides of the equations, given by 
 * L = 1 - P*V*G
 * R = P*V
 * Since G is expressed in an SWP basis, we also must include the basis-transormation matrices C */
void faddeev_dense_solver(cdouble*  U_array,
					      cdouble*  G_array,
					      int*		q_com_idx_array,	size_t num_q_com,
					      int*      deuteron_idx_array, size_t num_deuteron_states,
					      size_t    Nalpha,
					      size_t 	Nq_WP,
					      size_t 	Np_WP,
					      double**  CT_RM_array,
					      double**  VC_CM_array,
					      double*   P123_sparse_val_array,
					      int*      P123_sparse_row_array,
					      size_t*   P123_sparse_col_array,
					      size_t    P123_sparse_dim,
					      const tnf_kernel_context& tnf_ctx){
	
	/* Stores A and K arrays */
	bool store_A_array = true;
	bool store_K_array = true;
	bool store_U_array = true;

	/* Dense dimension of 3N-channel */
	size_t dense_dim = Nalpha * Nq_WP * Np_WP;
	
	std::complex<double>* L_array = new cdouble [dense_dim*dense_dim];
	std::complex<double>* R_array = new cdouble [dense_dim*dense_dim];

	double* CPVC_col_array 		   = new double [dense_dim];
	int*    CPVC_row_to_nnz_array  = new int    [dense_dim];
	int*    CPVC_nnz_to_row_array  = new int    [dense_dim];
	
	for (size_t j=0; j<num_q_com; j++){

		/* Reset L- and R-arrays */
		for (size_t idx=0; idx<dense_dim*dense_dim; idx++){
			L_array[idx] = 0;
			R_array[idx] = 0;
		}

		/* Construct L- and R-arrays */
		//#pragma omp parallel for
		for (size_t idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
			for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
				for (size_t idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){
					size_t col_idx = idx_alpha_c*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_c;

					/* Reset CPVC-column array */
					for (size_t row_idx=0; row_idx<dense_dim; row_idx++){
						CPVC_col_array[row_idx] = 0;
						CPVC_row_to_nnz_array[row_idx] = -1;
						CPVC_nnz_to_row_array[row_idx] = -1;
					}

					/* Calculate CPVC-column */
					size_t CPVC_num_nnz = 0;
					calculate_CPVC_col(CPVC_col_array,
									   CPVC_row_to_nnz_array,
									   CPVC_nnz_to_row_array,
									   CPVC_num_nnz,
									   idx_alpha_c, idx_p_c, idx_q_c,
									   Nalpha, Nq_WP, Np_WP,
									   CT_RM_array,
									   VC_CM_array,
									   P123_sparse_val_array,
									   P123_sparse_row_array,
									   P123_sparse_col_array,
									   P123_sparse_dim,
									   tnf_ctx);

    	    		for (size_t row_idx=0; row_idx<dense_dim; row_idx++){
						L_array[row_idx*dense_dim + col_idx] = -CPVC_col_array[row_idx]*G_array[j*dense_dim + col_idx];
						R_array[row_idx*dense_dim + col_idx] =  CPVC_col_array[row_idx];
    	    		}

    	    		L_array[col_idx*dense_dim + col_idx] += 1;
				}
			}
    	}

		if (store_A_array){
			std::string A_arr_filename = "A_array_E_idx_" + std::to_string(j) + ".txt";
			store_array(R_array, dense_dim*dense_dim, A_arr_filename);
		}
		if (store_K_array){
			std::string K_arr_filename = "K_array_E_idx_" + std::to_string(j) + ".txt";
			store_array(L_array, dense_dim*dense_dim, K_arr_filename);
		}

		/* Solve */
		solve_MM(L_array, R_array, dense_dim);

		if (store_U_array){
			std::string U_arr_filename = "U_array_E_idx_" + std::to_string(j) + ".txt";
			store_array(R_array, dense_dim*dense_dim, U_arr_filename);
		}
