// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_permutation_matrix.cpp
// 行号区段：1015..1220
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void read_and_merge_thread_files_to_single_array(double** P123_val_sparse_array,
												 int**    P123_row_array,
												 int**    P123_col_array,
												 size_t&  P123_dim,
												 int*	  max_TFC_array,
												 int	  P123_omp_num_threads,
												 int      Np_WP, double *p_array_WP_bounds,
												 int      Nq_WP, double *q_array_WP_bounds,
												 int      Nx, double* x_array, double* wx_array,
												 int      Nphi,
												 int      J_2N_max,
												 pw_3N_statespace pw_states,
												 std::string P123_folder){

	int Nalpha 	 = pw_states.Nalpha;
	int two_J_3N = pw_states.two_J_3N_array[0];
	int P_3N  	 = pw_states.P_3N_array[0];

	size_t P123_dense_dim = Nalpha * Nq_WP * Np_WP;

	std::vector<std::string> TF_filenames_vector;
	
	// [EN] During construction each OpenMP worker flushes its non-zero triplets independently, so the global sparse
	// matrix exists temporarily as many unordered COO fragments. This merge step restores a single reusable P123 object.
	// / [CN] 在构造阶段，每个 OpenMP 线程都会独立把自己的非零三元组刷盘，因此全局稀疏矩阵会暂时分裂成许多无序的 COO
	// 片段；这里的合并步骤负责把它们恢复成一个可复用的 P123 对象。
	printf("   - Checking parallel thread files exist and reading dimensions \n"); fflush(stdout);
	/* Read dimensions from sparse subarray files */
	for (int thread_idx=0; thread_idx<P123_omp_num_threads; thread_idx++){
		int TFC_max = max_TFC_array[thread_idx];
		for (int current_TFC=0; current_TFC<TFC_max; current_TFC++){
			/* Generate filename for current thread_idx and TFC */
			std::string thread_filename = generate_subarray_file_name(two_J_3N, P_3N,
																	  Np_WP, Nq_WP,
																	  J_2N_max,
																	  thread_idx,
																	  current_TFC,
																	  P123_folder);
			
			/* Verify that file exists */
			if (std::experimental::filesystem::exists(thread_filename)==false){
				printf("     - WARNING file not found: %s \n", thread_filename.c_str()); 
				printf("       Are you certain input P123 thread files are complete? \n");
				printf("       Program will ignore file and continue ... \n");
				fflush(stdout);
				continue;
			}

			/* Append filename since it exists */
			TF_filenames_vector.push_back(thread_filename);
				
			/* Convert std::string filename to char */
			char filename_char[300];
			std::strcpy(filename_char, thread_filename.c_str());
				
			/* Retrieve number of P123-elements in current file */
			unsigned long long int current_P123_sparse_dim = 0;
			read_ULL_integer_from_h5(current_P123_sparse_dim, "P123_sparse_dim", filename_char);

			P123_dim += current_P123_sparse_dim;
		}
	}
	printf("   - Total number of non-zero elements: %zu \n", P123_dim); fflush(stdout);
		
	/* Allocate required memory to fite whole P123-matrix */
	double required_mem = P123_dim * (sizeof(double) + 2*sizeof(int))/std::pow(2.0,30);
	printf("   - Allocating necessary arrays (requires %.2f GB) \n", required_mem); fflush(stdout);
	try{
		*P123_val_sparse_array = new double [P123_dim];
		*P123_row_array        = new int    [P123_dim];
		*P123_col_array        = new int    [P123_dim];
	}
	catch (...) {
		raise_error("Failed. Memory exceeded in sparse memory allocation");
	}

	/* Read elements from sparse subarray files */
	printf("   - Reading elements (in random ordering) from parallel thread files \n"); fflush(stdout);
	size_t nnz_counter = 0;
	for (size_t idx_TF=0; idx_TF<TF_filenames_vector.size(); idx_TF++){
		
		/* Retrieve filename from vector */
		std::string thread_filename = TF_filenames_vector[idx_TF];
				
		/* Retrieve P123-elements in current file and write to temporary arrays*/
		double* P123_sparse_val_subarray = NULL;
		int* 	P123_sparse_row_subarray = NULL;
		int* 	P123_sparse_col_subarray = NULL;
		size_t  current_P123_sparse_dim  = 0;
		read_sparse_permutation_matrix_for_3N_channel_h5(&P123_sparse_val_subarray,
													     &P123_sparse_row_subarray,
													     &P123_sparse_col_subarray,
													     current_P123_sparse_dim,
													     Np_WP, p_array_WP_bounds,
													     Nq_WP, q_array_WP_bounds,
													     pw_states,
											   		     thread_filename,
														 false);

		/* Append RANDOM ORDER elements to input arrays */
		for (size_t idx=0; idx<current_P123_sparse_dim; idx++){
			(*P123_val_sparse_array)[idx+nnz_counter] = P123_sparse_val_subarray[idx];
			(*P123_row_array)       [idx+nnz_counter] = P123_sparse_row_subarray[idx];
			(*P123_col_array)       [idx+nnz_counter] = P123_sparse_col_subarray[idx];
		}

		nnz_counter += current_P123_sparse_dim;

		delete [] P123_sparse_val_subarray;
		delete [] P123_sparse_row_subarray;
		delete [] P123_sparse_col_subarray;
	}

	// [EN] Downstream sparse kernels assume a canonical row-major COO ordering. Sorting changes only storage order,
	// not the numerical values of P123, but it makes later scans deterministic and cache-friendlier. / [CN] 后续
	// 稀疏核默认使用规范的 row-major COO 顺序；这里的排序只改变存储顺序，不改变 P123 的数值，但能让后续扫描更稳定、
	// 更符合缓存访问。
	printf("   - Sorting random-order P123-elements to row-major COO-format \n"); fflush(stdout);
	/* Index array to sort */
	printf("     - Creating index-vector \n"); fflush(stdout);
	std::vector<size_t> P123_idx_vector (P123_dim, 0);
	#pragma omp parallel
	for (size_t idx=0; idx<P123_dim; idx++){
		size_t coo_array_idx = (size_t)(*P123_row_array)[idx]*P123_dense_dim + (size_t)(*P123_col_array)[idx];
		if ((size_t)(*P123_row_array)[idx] > P123_dense_dim){
			raise_error("Encountered row-index larger than dense dimension!");
		}
		if ((size_t)(*P123_col_array)[idx] > P123_dense_dim){
			raise_error("Encountered column-index larger than dense dimension!");
		}
		if(coo_array_idx>P123_dense_dim*P123_dense_dim){
			raise_error("coo_array_idx=row_idx*dense_dim + col_idx was larger than dense_dim^2!");
		}
		P123_idx_vector[idx] = coo_array_idx;
	}

	/* Sort indices using template (Warning - lambda expression - requires C++11 or newer compiler) */
	printf("     - Creating sorting-vector \n"); fflush(stdout);
	std::vector<size_t> sorted_indices;
	try{
		sorted_indices = sort_indexes(P123_idx_vector);
	}
	catch(...){
		raise_error("Sorting failed!");
	}
	printf("       - Done \n"); fflush(stdout);
	printf("     - Verifying no index-overflow has occured ... \n"); fflush(stdout);
	for (size_t i=0; i<P123_dim; i++){
		if (sorted_indices[i]>=P123_dim){
			raise_error("Non-zero index of P123-array was larger than P123_dim");
		}
	}

	/* Free up memory - the next few steps are very memory intensive */
	P123_idx_vector.clear();

	///* Sorting test */
	//std::vector<int> test_vec (10, 0);;
	//for (int idx=0; idx<10; idx++){
	//	test_vec[idx] = -(5-idx);
	//	std::cout << "test_vec: " << test_vec[idx] << std::endl;
	//}
	//auto sorted_test = sort_indexes(test_vec);
	//for (int idx=0; idx<10; idx++){
	//	std::cout << "test_vec: " << test_vec[sorted_test[idx]] << std::endl;
	//}

	/* Sort row indices */
	printf("     - Sorting random-order row indices \n"); fflush(stdout);
	int* sparse_idx_array_temp = new int [P123_dim];
	#pragma omp parallel
	for (size_t i=0; i<P123_dim; i++){
		sparse_idx_array_temp[i] = (*P123_row_array)[sorted_indices[i]];
	}
	std::copy(sparse_idx_array_temp, sparse_idx_array_temp + P123_dim, *P123_row_array);

	//for (int i=0; i<P123_dim-1; i++){
	//	if (sparse_idx_array_temp[i]>sparse_idx_array_temp[i+1]){
	//		std::cout << i << " " << sparse_idx_array_temp[i] << " " << sparse_idx_array_temp[i+1] << std::endl;
	//		std::cout << i << " " << (*P123_row_array)[sorted_indices[i]] << " " << (*P123_row_array)[sorted_indices[i+1]] << std::endl;
	//		std::cout << i << " " << (*P123_row_array)[i] << " " << (*P123_row_array)[i+1] << std::endl;
	//		raise_error("row index error");
	//	}
	//}

	/* Sort column indices */
	printf("     - Sorting random-order column indices \n"); fflush(stdout);
	#pragma omp parallel
	for (size_t i=0; i<P123_dim; i++){
		sparse_idx_array_temp[i] = (*P123_col_array)[sorted_indices[i]];
	}
	std::copy(sparse_idx_array_temp, sparse_idx_array_temp + P123_dim, *P123_col_array);

	delete [] sparse_idx_array_temp;
	double* sparse_val_array_temp = new double [P123_dim];

	/* Sort values */
	printf("     - Sorting random-order P123-values \n"); fflush(stdout);
	#pragma omp parallel
	for (size_t i=0; i<P123_dim; i++){
		sparse_val_array_temp[i] = (*P123_val_sparse_array)[sorted_indices[i]];
	}
	std::copy(sparse_val_array_temp, sparse_val_array_temp + P123_dim, *P123_val_sparse_array);
	delete [] sparse_val_array_temp;
	printf("     - Done \n"); fflush(stdout);
}
