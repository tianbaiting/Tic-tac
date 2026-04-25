// ===============================================================
// 抽取自仓库 [current]: src/io/disk_io_routines.cpp
// 行号区段：1209..1284
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	/* Verify mesh points are equal to current program run, exit if not */
	if (Nalpha_file!=Nalpha || Np_WP_file!=Np_WP || Nq_WP_file!=Nq_WP){
		raise_error("File-read P123 state-space dimensions (Nalpha, Nq_WP, Np_WP) mismatch.");
	}

	/* Read p-momentum WP boundaries */
	double p_WP_array_file [Np_WP+1];
	read_WP_boundaries_from_h5(p_WP_array_file, Np_WP, "p boundaries", filename);
	/* Verify boundaries match current program run, exit if not */
	for (int i=0; i<Np_WP+1; i++){
		double p_boundary_prog = p_WP_array[i];
		double p_boundary_read = p_WP_array_file[i];
		if (p_boundary_read!=p_boundary_prog){
			raise_error("File-read P123 p-momentum boundaries mismatch.");
		}
	}

	/* Read q-momentum WP boundaries */
	double q_WP_array_file [Nq_WP+1];
	read_WP_boundaries_from_h5(q_WP_array_file, Nq_WP, "q boundaries", filename);
	/* Verify boundaries match current program run, exit if not */
	for (int i=0; i<Nq_WP+1; i++){
		double q_boundary_prog = q_WP_array[i];
		double q_boundary_read = q_WP_array_file[i];
		if (q_boundary_read!=q_boundary_prog){
			raise_error("File-read P123 q-momentum boundaries mismatch.");
		}
	}

	/* Read PW state space */
	int L_2N_array_file     [Nalpha];
	int S_2N_array_file     [Nalpha];
	int J_2N_array_file     [Nalpha];
	int T_2N_array_file     [Nalpha];
	int L_1N_array_file     [Nalpha];
	int two_J_1N_array_file [Nalpha];
	int two_T_3N_array_file [Nalpha];
	int two_J_3N_array_file [Nalpha];
	int P_3N_array_file     [Nalpha];
	read_PW_statespace_to_h5(Nalpha,
							 L_2N_array_file,
							 S_2N_array_file,
							 J_2N_array_file,
							 T_2N_array_file,
							 L_1N_array_file, 
							 two_J_1N_array_file,
							 two_T_3N_array_file,
							 two_J_3N_array_file,
							 P_3N_array_file,
							 filename);
	/* Verify PW statespace match current program run, exit if not */
	for (int i=0; i<Nalpha; i++){
		if (L_2N_array_file[i]     != L_2N_array[i] ||
			S_2N_array_file[i]     != S_2N_array[i] ||
			J_2N_array_file[i]     != J_2N_array[i] ||
			T_2N_array_file[i]     != T_2N_array[i] ||
			L_1N_array_file[i]     != L_1N_array[i] ||
			two_J_1N_array_file[i] != two_J_1N_array[i] ||
			two_T_3N_array_file[i] != two_T_3N_array[i] ||
			two_J_3N_array_file[i] != two_J_3N ||
			P_3N_array_file[i]     != P_3N){
			raise_error("File-read P123 PW state-space mismatch.");
		}
	}

	/* Read P123 sparse matrix elements and indices */
	*P123_sparse_row_array = new int    [P123_sparse_dim];
	*P123_sparse_col_array = new int    [P123_sparse_dim];
	*P123_sparse_val_array = new double [P123_sparse_dim];

	read_sparse_permutation_matrix_h5(*P123_sparse_val_array,
									  *P123_sparse_row_array,
									  *P123_sparse_col_array,
									  P123_sparse_dim,
									  filename);
}
