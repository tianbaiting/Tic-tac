// ===============================================================
// 抽取自仓库 [origin]: CPP/solve_faddeev.cpp
// 行号区段：271..296
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
cdouble pade_approximant(cdouble* a_coeff_array, size_t N, size_t M, cdouble z){

	/* a_coeff_array must have length N+M+1 */
	cdouble P_array [(M+1)*(M+1)];
	cdouble Q_array [(M+1)*(M+1)];

	for (size_t row_idx=0; row_idx<M; row_idx++){
		for (size_t col_idx=0; col_idx<M+1; col_idx++){
			P_array[row_idx*(M+1) + col_idx] = a_coeff_array[N-M+1 + row_idx + col_idx];
			Q_array[row_idx*(M+1) + col_idx] = a_coeff_array[N-M+1 + row_idx + col_idx];
		}
	}

	for (size_t col_idx=0; col_idx<M+1; col_idx++){
		Q_array[M*(M+1) + col_idx] = std::pow(z, M-col_idx);
		P_array[M*(M+1) + col_idx] = 0;
		for (size_t j=M-col_idx; j<N+1; j++){
			P_array[M*(M+1) + col_idx] += a_coeff_array[j - (M-col_idx)] * std::pow(z, j);
		}
	}

	cdouble P_det = determinant(P_array, M+1);
	cdouble Q_det = determinant(Q_array, M+1);

	return P_det/Q_det;
}
