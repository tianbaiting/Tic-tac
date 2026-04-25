// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_swp_states.cpp
// 行号区段：84..98
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void construct_free_hamiltonian(double* H0_WP_array,
							    int Np_WP, double* p_WP_array){
	
    /* Loop over p-momenta */
    for (int idx_p=0; idx_p<Np_WP; idx_p++){
		/* Kinetic energy at bin-boundaries */
		double p1 = p_WP_array[idx_p];
		double p2 = p_WP_array[idx_p+1];
	
		/* Free Hamiltonian for 2-nucleon pair for momentum WPs */
		H0_WP_array[idx_p] = (p2*p2 + p2*p1 + p1*p1)/(6*mu23);
		/* Free Hamiltonian for 2-nucleon pair for energy WPs */
        //H0_WP_array[idx_p] = (p2*p2 + p1*p1)/(2*MN);
	}
}
