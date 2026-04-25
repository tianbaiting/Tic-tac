// ===============================================================
// 抽取自仓库 [current]: src/utils/auxiliary.cpp
// 行号区段：863..948
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
double Atilde (int alpha, int alphaprime, int Ltotal, int Jj_dim, int *L12_Jj, int *l3_Jj, int *J12_Jj, int *two_j3_Jj, int *S12_Jj, int *T12_Jj, int two_J, int two_T, double *SixJ_array, int two_jmax_SixJ)
{

	double ret = 0.0;

	for (int two_Stotal = 1; two_Stotal <= 3; two_Stotal += 2)
	{

		ret +=          sqrt(    (2 * J12_Jj[alpha] + 1) *     (two_j3_Jj[alpha] + 1) *     (2 * S12_Jj[alpha] + 1) *     (2 * T12_Jj[alpha] + 1)
								 * (2 * J12_Jj[alphaprime] + 1) * (two_j3_Jj[alphaprime] + 1) * (2 * S12_Jj[alphaprime] + 1) * (2 * T12_Jj[alphaprime] + 1)
							)
						* gsl_sf_pow_int(-1, S12_Jj[alphaprime] + T12_Jj[alphaprime])
						* (two_Stotal + 1)
						* gsl_sf_coupling_9j(2 * L12_Jj[alpha],     2 * S12_Jj[alpha],     2 * J12_Jj[alpha],     2 * l3_Jj[alpha],     1, two_j3_Jj[alpha],     2 * Ltotal, two_Stotal, two_J)
						* gsl_sf_coupling_9j(2 * L12_Jj[alphaprime], 2 * S12_Jj[alphaprime], 2 * J12_Jj[alphaprime], 2 * l3_Jj[alphaprime], 1, two_j3_Jj[alphaprime], 2 * Ltotal, two_Stotal, two_J)
						* SixJSymbol(SixJ_array, two_jmax_SixJ, 1, 1, 2 * S12_Jj[alpha], 1, two_Stotal, 2 * S12_Jj[alphaprime])
						* SixJSymbol(SixJ_array, two_jmax_SixJ, 1, 1, 2 * T12_Jj[alpha], 1, two_T, 2 * T12_Jj[alphaprime]);

	}

	return ret;
}

double Gtilde_new (double p, double q, double x, int alpha, int alphaprime, int N_alpha, int Lmax, int *L12_Jj, int *l3_Jj, double *A_store, int two_Jtotal)
{

	double ret = 0.0;

	double fac1, fac2;

	double pi1 = pi1_tilde(p, q, x);
	double pi2 = pi2_tilde(p, q, x);

	double costheta1 = -(0.5 * p + 0.75 * q * x) / pi1;
	double costheta2 = (p - 0.5 * q * x) / pi2;

	/* Prevent numerical error in Plm */
	if ( costheta1>1 ){
		costheta1 = 1;
	}
	else if( costheta1<-1 ){
		costheta1 = -1;
	}
	
	/* Prevent numerical error in Plm */
	if ( costheta2>1 ){
		costheta2 = 1;
	}
	else if( costheta2<-1 ){
		costheta2 = -1;
	}

	int L12 = L12_Jj[alpha];
	int l3 = l3_Jj[alpha];
	int L12prime = L12_Jj[alphaprime];
	int l3prime = l3_Jj[alphaprime];

	for (int Ltotal = max(abs(L12 - l3), abs(L12prime - l3prime)); Ltotal <= min((two_Jtotal + 5) / 2, min(L12 + l3, L12prime + l3prime)); Ltotal++)
	{
		fac1 = 8.0 * M_PI * M_PI * A_store[alpha * N_alpha * (Lmax + 1) + alphaprime * (Lmax + 1) + Ltotal];
		
		for (int Mtotal = -min(l3, Ltotal); Mtotal <= min(l3, Ltotal); Mtotal++)
		{

			fac2 = ClebschGordan(2 * L12, 2 * l3, 2 * Ltotal, 0, 2 * Mtotal, 2 * Mtotal)
				   * sqrt((2.0 * L12 + 1) / (4 * M_PI))
				   * gsl_sf_pow_int(-1, Mtotal)
				   * Plm(l3, Mtotal, x); // -1^M phase since azimutal angles of p' and q' = pi

			for (int M12primesum = -L12prime; M12primesum <= L12prime; M12primesum++)
			{
				if (abs(Mtotal - M12primesum) <= l3prime)
				{
					ret += fac1
						   * fac2
						   * ClebschGordan(2 * L12prime, 2 * l3prime, 2 * Ltotal, 2 * M12primesum, 2 * Mtotal - 2 * M12primesum, 2 * Mtotal)
						   * Plm(L12prime, M12primesum, costheta1)
						   * Plm(l3prime, Mtotal - M12primesum, costheta2);
				}
			}
		}
	}

	return ret;

}
