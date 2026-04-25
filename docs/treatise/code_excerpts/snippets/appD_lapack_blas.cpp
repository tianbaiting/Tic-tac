// ===============================================================
// 抽取自仓库 [current]: CPP/General_functions/matrix_routines.cpp
// 行号区段：4..75
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void dot_MV(double *A, double *B, double *C, int N, int M){
	double  beta  = 0.0;
	double  alpha = 1.0;
	int incrx = 1;
	int incry = 1;

    // Assuming A is M x N matrix, B is vector of size N, C is vector of size M
    // C = alpha*A*B + beta*C
	cblas_dgemv(CblasRowMajor, CblasNoTrans, M, N, alpha, A, N, B, incrx, beta, C, incry);
}


std::complex<double> cdot_VV(std::complex<float> *X, std::complex<float> *Y, int N, int INCR_X, int INCR_Y){
	std::complex<float> dot_product = 0;
	cblas_cdotu_sub(N, reinterpret_cast<const float*>(X), INCR_X, reinterpret_cast<const float*>(Y), INCR_Y, &dot_product);
	return dot_product;
}
std::complex<double> cdot_VV(std::complex<double> *X, std::complex<double> *Y, int N, int INCR_X, int INCR_Y){
	std::complex<double> dot_product = 0;
	cblas_zdotu_sub(N, reinterpret_cast<const double*>(X), INCR_X, reinterpret_cast<const double*>(Y), INCR_Y, &dot_product);
	return dot_product;
}

void dot_MM(float *A, float *B, float *C, int N, int K, int M){
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, 1.0, A, K, B, M, 0.0, C, M);
}
void dot_MM(double *A, double *B, double *C, int N, int K, int M){
	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, 1.0, A, K, B, M, 0.0, C, M);
}

void cdot_MM(std::complex<float> *A, std::complex<float> *B, std::complex<float> *C, int N, int K, int M){
	std::complex<float> beta = {0,0};
	std::complex<float> alpha = {1,0};
	cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, &alpha, A, K, B, M, &beta, C, M);
}
void cdot_MM(std::complex<double> *A, std::complex<double> *B, std::complex<double> *C, int N, int K, int M){
	std::complex<double> beta = {0,0};
	std::complex<double> alpha = {1,0};
	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, M, K, &alpha, A, K, B, M, &beta, C, M);
}

void solve_MM(float* A, float* B, int dim){
	char trans = 'N';
	lapack_int ipiv [dim];
	
	LAPACKE_sgetrf(LAPACK_ROW_MAJOR, dim, dim, A, dim, ipiv);
	LAPACKE_sgetrs(LAPACK_ROW_MAJOR, trans, dim, dim, A, dim, ipiv, B, dim);
}
void solve_MM(double* A, double* B, int dim){
	char trans = 'N';
	lapack_int ipiv [dim];
	
	LAPACKE_dgetrf(LAPACK_ROW_MAJOR, dim, dim, A, dim, ipiv);
	LAPACKE_dgetrs(LAPACK_ROW_MAJOR, trans, dim, dim, A, dim, ipiv, B, dim);
}
void solve_MM(std::complex<float> *A, std::complex<float> *B, int N){
	
	char trans = 'N';
	lapack_int ipiv [N];
	
	LAPACKE_cgetrf(LAPACK_ROW_MAJOR, N, N, reinterpret_cast<lapack_complex_float*>(A), N, ipiv);
	LAPACKE_cgetrs(LAPACK_ROW_MAJOR, trans, N, N, reinterpret_cast<const lapack_complex_float*>(A), N, ipiv, reinterpret_cast<lapack_complex_float*>(B), N);
}
void solve_MM(std::complex<double> *A, std::complex<double> *B, int N){
	
	char trans = 'N';
	lapack_int ipiv [N];
	
	LAPACKE_zgetrf(LAPACK_ROW_MAJOR, N, N, reinterpret_cast<lapack_complex_double*>(A), N, ipiv);
	LAPACKE_zgetrs(LAPACK_ROW_MAJOR, trans, N, N, reinterpret_cast<const lapack_complex_double*>(A), N, ipiv, reinterpret_cast<lapack_complex_double*>(B), N);
}

