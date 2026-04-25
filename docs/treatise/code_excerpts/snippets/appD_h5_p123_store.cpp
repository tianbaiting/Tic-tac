// ===============================================================
// 抽取自仓库 [current]: src/io/disk_io_routines.cpp
// 行号区段：1078..1165
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void store_sparse_permutation_matrix_for_3N_channel_h5(double* P123_sparse_val_array,
													   int*    P123_sparse_row_array,
													   int*    P123_sparse_col_array,
													   size_t  P123_sparse_dim,
													   int     Np_WP, double* p_WP_array,
													   int     Nq_WP, double* q_WP_array,
													   pw_3N_statespace pw_states,
													   std::string filename_in,
													   bool print_content){
	
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

	if (print_content){
		printf("   - Setting up h5-file \n");
	}

	/* Convert filename_in to char-array */
	char filename[300];
	std::strcpy(filename, filename_in.c_str());

	if (print_content){
		cout << " - Write to: " << filename << "\n";
	}
	
	/* Create and open file */
	hid_t file_id = H5Fcreate(filename,
							  H5F_ACC_TRUNC,
							  H5P_DEFAULT,
							  H5P_DEFAULT);

	/* Write number of mesh points (Nalpha, Np_WP, Nq_WP, and P123_sparse_dim) */
	write_integer_to_h5(Nalpha,          "Nalpha",          file_id);
	write_integer_to_h5(Np_WP,           "Np_WP",           file_id);
	write_integer_to_h5(Nq_WP,           "Nq_WP",           file_id);

	unsigned long long int P123_sparse_dim_temp = P123_sparse_dim;
	write_ULL_integer_to_h5(P123_sparse_dim_temp, "P123_sparse_dim", file_id);

	/* Write p-momentum WP boundaries */
	if (print_content){
		printf("   - Writing p-momentum bins \n");
	}
	write_WP_boundaries_to_h5(p_WP_array, Np_WP, "p boundaries", file_id);
	/* Write q-momentum WP boundaries */
	if (print_content){
		printf("   - Writing q-momentum bins \n");
	}
	write_WP_boundaries_to_h5(q_WP_array, Nq_WP, "q boundaries", file_id);

	/* PW quantum numbers */
	if (print_content){
		printf("   - Writing partial-wave state space \n");
	}
	write_PW_statespace_to_h5(Nalpha,
							  L_2N_array,
							  S_2N_array,
							  J_2N_array,
							  T_2N_array,
							  L_1N_array, 
							  two_J_1N_array,
							  two_T_3N_array,
							  two_J_3N,
							  P_3N,
							  file_id);
	
	/* Sparse matrix elements */
	if (print_content){
		printf("   - Writing P123 sparse matrix elements and indices \n");
	}
	write_sparse_permutation_matrix_h5(P123_sparse_val_array,
									   P123_sparse_row_array,
									   P123_sparse_col_array,
									   P123_sparse_dim,
									   file_id);

	herr_t status = H5Fclose(file_id);
}

