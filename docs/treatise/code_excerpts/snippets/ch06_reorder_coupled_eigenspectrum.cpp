// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_swp_states.cpp
// 行号区段：142..169
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void reorder_coupled_eigenspectrum(double* eigenvalues,
								   double* eigenvectors,
								   int     Np_WP){
	
	/* Temporary arrays to hold values as we reorder indices */
	std::vector<double> buffer_array_vals(2*Np_WP);
	std::vector<double> buffer_array_vecs(4*Np_WP*Np_WP);

	for (int i=0; i<Np_WP; i++){
		/* Temporarily store reorder of eigenvalues */
		buffer_array_vals[i] 	   = eigenvalues[2*i];
		buffer_array_vals[i+Np_WP] = eigenvalues[2*i+1];

		for (int j=0; j<2*Np_WP; j++){
			/* Temporarily store reorder of eigenvectors */
			buffer_array_vecs[j*2*Np_WP + i] 	     = eigenvectors[j*2*Np_WP + 2*i];
			buffer_array_vecs[j*2*Np_WP + i + Np_WP] = eigenvectors[j*2*Np_WP + 2*i+1];
		}
	}

	/* Write reordered elements back into original arrays */
	for (int i=0; i<2*Np_WP; i++){
		eigenvalues[i] = buffer_array_vals[i];
		for (int j=0; j<2*Np_WP; j++){
			eigenvectors[i*2*Np_WP + j] = buffer_array_vecs[i*2*Np_WP + j];
		}
	}
}
