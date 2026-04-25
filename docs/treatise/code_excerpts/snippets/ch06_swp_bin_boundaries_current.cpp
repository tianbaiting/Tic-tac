// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_swp_states.cpp
// 行号区段：206..246
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void make_swp_bin_boundaries(double* eigenvalues,
							 double* e_SWP_array,
							 int	 Np_WP,
							 bool    coupled,
							 bool    chn_3S1){

	// [EN] Neighboring eigenvalues define the edges of one interacting packet cell. For the deuteron channel the
	// first cell is split so the bound state sits on its own branch below threshold. / [CN] 相邻本征值共同定义一个
	// 相互作用波包单元的边界；对氘核通道而言，第一个单元会被特殊处理，使束缚态单独占据阈值以下的一条分支。
	/* Set first boundary */
	e_SWP_array[0] = 0;
	if (coupled){
		e_SWP_array[Np_WP+1] = 0;
	}
	
	/* Set all the middle boundaries */
	for (int idx_p=0; idx_p<Np_WP-1; idx_p++){
		e_SWP_array[idx_p+1] = 0.5*(eigenvalues[idx_p+1] + eigenvalues[idx_p]);

		if (coupled){
			e_SWP_array[idx_p+1+Np_WP+1] = 0.5*(eigenvalues[idx_p+1+Np_WP] + eigenvalues[idx_p+Np_WP]);
		}
	}

	/* Set end boundary */
	e_SWP_array[Np_WP] = eigenvalues[Np_WP-1] + 0.5*(e_SWP_array[Np_WP-1] - e_SWP_array[Np_WP-2]);
	if (coupled){
		e_SWP_array[2*Np_WP+1] = eigenvalues[2*Np_WP-1] + 0.5*(e_SWP_array[2*Np_WP] - e_SWP_array[2*Np_WP-1]);
	}

	/* Modify first boundary for 3S1 */
	if (chn_3S1){
		e_SWP_array[0] = eigenvalues[0];
		e_SWP_array[1] = 0;
	}

	//for (int idx_p=0; idx_p<Np_WP+1; idx_p++){
	//	std::cout << e_SWP_array[idx_p] << std::endl;
	//}
	//std::cout << std::endl;
}
