// ===============================================================
// 抽取自仓库 [current]: CPP/General_functions/gauss_legendre.cpp
// 行号区段：36..67
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void gauss(double* x, double* w, int N){
	
	double p3	= 0;								// Legendre polynomial, P(xi), for i=3
	double eps = 3.e-16;
	int m	= (N+1)/2;
	
	for (int i=1; i<m+1; i++){
		double t  = cos(M_PI*(i-0.25)/(N+0.5));		// xi
		double t1 = 1;								// old xi
		double pp = 0;
		
		while (std::abs(t-t1) >= eps){
			double p1 = 1;							// Legendre polynomial, P(xi), for i=1
			double p2 = 0;							// Legendre polynomial, P(xi), for i=2
			
			for (int j=1; j<N+1; j++){
				p3 = p2;
				p2 = p1;
				p1 = ((2*j-1)*t*p2 - (j-1)*p3)/j;	// recurrence relation for Legendre polynomials
			}
			pp = N*(t*p1-p2)/(t*t-1);				// identity for P'(xi)
			t1 = t;
			t  = t1 - p1/pp; 						// Newton's method for finding roots
		}
		
		x[i-1] = -t;
		x[N-i] = t;
		
		w[i-1] = 2./((1-t*t)*pp*pp);
		w[N-i] = w[i-1];
	}
}
