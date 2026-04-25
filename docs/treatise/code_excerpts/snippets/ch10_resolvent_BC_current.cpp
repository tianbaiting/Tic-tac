// ===============================================================
// 抽取自仓库 [current]: src/core/resolvent/make_resolvent.cpp
// 行号区段：60..85
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
// [EN] In the SWP basis the channel resolvent becomes diagonal, and the only remaining work is to integrate the
// resolvent kernel over the packet energy cells. For a deuteron-like bound channel this yields a bound-continuum
// term with the expected logarithmic real part and step-function imaginary part. / [CN] 在 SWP 基中，通道分辨算符变为对角形式；剩余工作只是把分辨核在波包能量元上积分。对于氘核这类束缚分支，会得到具有对数实部和阶跃函数虚部的束缚-连续项。
cdouble resolvent_bound_continuum(double E, double Eb,
								  double q_bin_upper,
								  double q_bin_lower){

	double mu1 = Mn*(Mn+Mp+Eb)/(Mn + Mn + Mp + Eb);

	/* Bin boundaries in energy */
	double Eq_lower = q_bin_lower*q_bin_lower/(2*mu1);
	double Eq_upper = q_bin_upper*q_bin_upper/(2*mu1);

	/* Energy width of SWPs (D is short for Delta) for q momenta */
	double Dq = Eq_upper - Eq_lower;

	/* Real part of the BC resolvent */
	double Re_R = std::log( std::abs( (Eq_lower + Eb - E)/(Eq_upper + Eb - E) ) ) / (Dq); 
	
	/* Imaginary part of the BC resolvent */
	double Im_R = (  heaviside_step_function( Eq_upper + Eb - E )
				   - heaviside_step_function( Eq_lower + Eb - E ) ) * (-M_PI)/ (Dq); 
	
	/* Return the complex BC resolvent term Q */
	return {Re_R, Im_R};
}
