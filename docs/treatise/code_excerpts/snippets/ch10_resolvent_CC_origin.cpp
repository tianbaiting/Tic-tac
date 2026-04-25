// ===============================================================
// 抽取自仓库 [origin]: CPP/make_resolvent.cpp
// 行号区段：38..80
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
/* See header-file, commentary (A), for explanation of notation and equations */
cdouble resolvent_continuum_continuum(double E, double Eb,
									  double q_bin_upper,
									  double q_bin_lower,
									  double e_bin_upper,
									  double e_bin_lower){

	double mu1 = Mn*(Mn+Mp+Eb)/(Mn + Mn + Mp + Eb);

	/* Bin boundaries in energy */
	double Eq_lower = q_bin_lower*q_bin_lower/(2*mu1);
	double Eq_upper = q_bin_upper*q_bin_upper/(2*mu1);

	double Ep_lower = e_bin_lower;
	double Ep_upper = e_bin_upper;

	double Ep = 0.5*(Ep_lower + Ep_upper);
	double Eq = 0.5*(Eq_lower + Eq_upper);

	/* Energy widths of SWPs (D is short for Delta) for p and q momenta */
	double Dq = Eq_upper - Eq_lower;
	double Dp = Ep_upper - Ep_lower;

	/* Temporary variables */
	double DM = 0.5*(Dp - Dq);
	double DP = 0.5*(Dp + Dq);
	double D  = Ep + Eq - E;
	
	/* Real part of the CC resolvent */
	double Re_Q = (  (D+DM) * std::log( std::abs(D+DM) )
				   + (D-DM) * std::log( std::abs(D-DM) )
				   - (D+DP) * std::log( std::abs(D+DP) )
				   - (D-DP) * std::log( std::abs(D-DP) ) ) / (Dp*Dq); 
	
	/* Imaginary part of the CC resolvent */
	double Im_Q = (  (D+DM) * heaviside_step_function( D+DM )
				   + (D-DM) * heaviside_step_function( D-DM )
				   - (D+DP) * heaviside_step_function( D+DP )
				   - (D-DP) * heaviside_step_function( D-DP ) ) * M_PI/ (Dp*Dq); 
	
	/* Return the complex CC resolvent term Q */
	return {Re_Q, Im_Q};
}
