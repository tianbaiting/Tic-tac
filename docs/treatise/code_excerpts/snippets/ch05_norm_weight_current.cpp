// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_wp_states.cpp
// 行号区段：64..85
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
double p_normalization(double p0, double p1){
	//double Ep0 = p0*p0/MN;
	//double Ep1 = p1*p1/MN;
	//return sqrt(Ep1-Ep0);	// energy WPs
	return sqrt(p1-p0);	// momentum WPs
}
double q_normalization(double q0, double q1){
	//double Eq0 = com_q_momentum_to_com_energy(q0);
	//double Eq1 = com_q_momentum_to_com_energy(q1);
	//return sqrt(Eq1-Eq0);	// energy WPs
	//return sqrt( (q1*q1-q0*q0)/MN );	// energy WPs (old)
	return sqrt(q1-q0);	// momentum WPs
}

double p_weight_function(double p){
	//return sqrt(2*p/MN);	// energy WPs
	return 1;				// momentum WPs
}
double q_weight_function(double q){
	//return sqrt(2*q/MN);	// energy WPs
	return 1;				// momentum WPs
}
