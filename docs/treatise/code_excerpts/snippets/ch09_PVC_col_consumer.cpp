// ===============================================================
// 抽取自仓库 [current]: src/core/faddeev_solver/solve_faddeev.cpp
// 行号区段：146..193
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
// [EN] This is the raw driving column PVC that appears in the AGS kernel before the left basis rotation by C^T.
// Applying the sparse permutation first keeps the expensive part of the kernel sparse. / [CN] 这里计算的是 AGS 核中
// 左乘 C^T 之前的原始驱动列 PVC；先施加稀疏置换算符，可以把核中最昂贵的部分保持为稀疏结构。
void calculate_PVC_col(double*  col_array,
					   size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
					   size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
					   double** VC_CM_array,
					   double*  P123_val_array,
					   int*  	P123_row_array,
					   size_t*  P123_col_array,
					   size_t   P123_dim){

	double* VC_subarray     = NULL;

	size_t dense_dim = Nalpha*Np_WP*Nq_WP;

	for (size_t idx_alpha_j=0; idx_alpha_j<Nalpha; idx_alpha_j++){
		VC_subarray = VC_CM_array[idx_alpha_c*Nalpha + idx_alpha_j];

		/* Only do inner-product if VC is not zero due to conservation laws */
		if (VC_subarray!=NULL){
			for (size_t idx_p_j=0; idx_p_j<Np_WP; idx_p_j++){

				/* Access VC element */
				double VC_element = VC_subarray[idx_p_c*Np_WP + idx_p_j];

				/* Inner-product index */
				size_t idx_j =  idx_alpha_j*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_j;

				/* CSC-format indexing */
				size_t idx_i_lower = P123_col_array[idx_j    ];
				size_t idx_i_upper = P123_col_array[idx_j + 1];

				/* Loop through rows of column we're calculating, and append */
				for (size_t idx_i=idx_i_lower; idx_i<idx_i_upper; idx_i++){

					/* Access arrays, this whole function is written to minimize these two calls */
					/* NOTE THAT P = P123 + P132 = 2*P123 FOR ANTISYMMETRIC PAIR-STATES */
					double P_element  = 2*P123_val_array[idx_i];

					size_t idx_row = P123_row_array[idx_i];
					/* Write inner product to col_array of PVC-product */
					col_array[idx_row] += P_element * VC_element;
				}
			}
		}
	}
}
