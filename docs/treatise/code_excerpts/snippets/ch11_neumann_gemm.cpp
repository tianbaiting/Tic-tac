// ===============================================================
// 抽取自仓库 [current]: src/core/faddeev_solver/solve_faddeev.cpp
// 行号区段：1444..1505
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
					time_CPVC_cols += timestamp_CPVC_chunk_end - timestamp_CPVC_chunk_start;

					double beta  = 0;
					double alpha = 1;
					int M    = num_non_conv_rows;// num_on_shell_A_rows
					int N    = cols_in_chunk;// max_num_cols_in_mem;
					int K    = dense_dim;
					int lda  = dense_dim;
					int ldb  = dense_dim;//max_num_cols_in_mem;
					int ldc  = dense_dim;
					double* re_A = &re_A_An_row_array_comp[0];
					double* im_A = &im_A_An_row_array_comp[0];
					double* B 	 = &CPVC_cols_array[0];
					double* re_C = &re_A_An_row_array_prod[idx_col_start];
					double* im_C = &im_A_An_row_array_prod[idx_col_start];

					double timestamp_gemm_start = omp_get_wtime();
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, re_A, lda, B, ldb, beta, re_C, ldc);	// real multiplication
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, im_A, lda, B, ldb, beta, im_C, ldc);	// imag multiplication
					double timestamp_gemm_end   = omp_get_wtime();
					time_An_CPVC_multiply += timestamp_gemm_end - timestamp_gemm_start;
				}
			}
			else{
				double timestamp_gemm_start = omp_get_wtime();
				const char   ordering = 'R';
				const char   trans 	  = 'T';
				const double alpha 	  = 1.0;
				/* Transpose An before sparse multiplication */
				inplace_transpose(re_A_An_row_array_comp, num_non_conv_rows, dense_dim);
				inplace_transpose(im_A_An_row_array_comp, num_non_conv_rows, dense_dim);
				/* Multiply CPVC with An using sparse multiplication */
				//dot_MM_sparse(CPVC_v_array, CPVC_c_array_LL, CPVC_csc_array_LL, re_A_An_row_array_comp, re_A_An_row_array_prod, dense_dim, dense_dim, num_non_conv_rows, true);
				//dot_MM_sparse(CPVC_v_array, CPVC_c_array_LL, CPVC_csc_array_LL, re_A_An_row_array_comp, im_A_An_row_array_prod, dense_dim, dense_dim, num_non_conv_rows, true);
				/* Transpose An+1 after sparse multiplication */
				inplace_transpose(re_A_An_row_array_prod, dense_dim, num_non_conv_rows);
				inplace_transpose(im_A_An_row_array_prod, dense_dim, num_non_conv_rows);
				double timestamp_gemm_end   = omp_get_wtime();
				time_An_CPVC_multiply += timestamp_gemm_end - timestamp_gemm_start;
			}

			/* Write compact format back to full format */
			for (size_t r=0; r<num_non_conv_rows; r++){
				size_t idx_row_NDOS = A_An_indexing_array[r];
				for (size_t i=0; i<dense_dim; i++){
					re_A_An_row_array[idx_row_NDOS*dense_dim + i] = re_A_An_row_array_prod[r*dense_dim + i];
					im_A_An_row_array[idx_row_NDOS*dense_dim + i] = im_A_An_row_array_prod[r*dense_dim + i];
				}
			}
			double timestamp_neumann_end = omp_get_wtime();
			time_neumann = timestamp_neumann_end - timestamp_neumann_start;

			printf("         - Time multiplying An with G:    %.6f s \n", time_resolvent);
			printf("         - Time generating CPVC-cols:     %.6f s \n", time_CPVC_cols);
			printf("         - Time multiplying An with CPVC: %.6f s \n", time_An_CPVC_multiply);
			printf("         - Total time:                    %.6f s \n", time_neumann);
			printf("         - Done \n"); fflush(stdout);

			/* Rewrite previous A_An with current A_An */
			for (size_t i=0; i<num_on_shell_A_rows*dense_dim; i++){
				re_A_An_row_array_prev[i] = re_A_An_row_array[i];
				im_A_An_row_array_prev[i] = im_A_An_row_array[i];
