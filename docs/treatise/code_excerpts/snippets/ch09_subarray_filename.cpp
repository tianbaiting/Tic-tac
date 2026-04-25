// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_permutation_matrix.cpp
// 行号区段：185..200
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
std::string generate_subarray_file_name(int two_J_3N, int P_3N,
										int Np_WP, int Nq_WP,
										int J_2N_max,
										int thread_idx,
										int current_TFC,
										std::string P123_folder){

	std::string filename = P123_folder + "/P123_subsparse_JP_"
						 + to_string(two_J_3N) + "_" + to_string(P_3N)
						 + "_Np_" + to_string(Np_WP) + "_Nq_" + to_string(Nq_WP)
						 + "_J2max_" + to_string(J_2N_max) + "_TFC_" + to_string(thread_idx)
						 + "_" + to_string(current_TFC) + ".h5";

	return filename;
}

