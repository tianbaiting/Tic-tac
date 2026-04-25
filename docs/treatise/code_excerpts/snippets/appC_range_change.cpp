// ===============================================================
// 抽取自仓库 [current]: CPP/General_functions/gauss_legendre.cpp
// 行号区段：76..95
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void rangeChange_0_inf(double* x, double* w, double scale, int N){
	for (int i=0; i<N; i++){
		w[i] = scale*w[i]*M_PI/(4*(cos(M_PI*(x[i]+1)/4))*(cos(M_PI*(x[i]+1)/4)));
		x[i] = scale*tan(M_PI*(x[i]+1)/4);
	}
}

// changes range from (-1,1) to (a,b)
void updateRange_a_b(float* x, float* w, float a, float b, int N){
	for (int i=0; i<N; i++){
		x[i] = 0.5*(b-a)*x[i] + 0.5*(b+a);
		w[i] = 0.5*(b-a)*w[i];
	}
}
void updateRange_a_b(double* x, double* w, double a, double b, int N){
	for (int i=0; i<N; i++){
		x[i] = 0.5*(b-a)*x[i] + 0.5*(b+a);
		w[i] = 0.5*(b-a)*w[i];
	}
}
