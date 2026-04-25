// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_swp_states.cpp
// 行号区段：70..82
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
