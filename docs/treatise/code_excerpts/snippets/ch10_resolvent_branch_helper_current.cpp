// ===============================================================
// 抽取自仓库 [current]: src/core/resolvent/make_resolvent.cpp
// 行号区段：4..49
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
namespace {

bool is_singlet_s_wave_channel(int L_2N, int S_2N, int J_2N){
	return L_2N==0 && S_2N==0 && J_2N==0;
}

double* select_swp_energy_branch(int L_2N,
								 int S_2N,
								 int J_2N,
								 int T_2N,
								 int two_T_3N,
								 int Np_WP,
								 double* e_SWP_unco_array,
								 double* e_SWP_coup_array,
								 const run_params& run_parameters){
	const bool coupled_via_L_2N = run_parameters.tensor_force && L_2N!=J_2N && J_2N!=0;
	const bool coupled_via_T_3N = is_singlet_s_wave_channel(L_2N, S_2N, J_2N) && run_parameters.isospin_breaking_1S0;
	if (coupled_via_L_2N && coupled_via_T_3N){
		raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
	}

	if (!coupled_via_L_2N && !coupled_via_T_3N){
		const int chn_2N_idx = unique_2N_idx(L_2N, S_2N, J_2N, T_2N, false, run_parameters);
		return &e_SWP_unco_array[chn_2N_idx * (Np_WP+1)];
	}

	const int chn_2N_idx = unique_2N_idx(L_2N, S_2N, J_2N, T_2N, true, run_parameters);
	double* channel_branch_ptr = &e_SWP_coup_array[chn_2N_idx * 2*(Np_WP+1)];

	// [EN] The coupled SWP storage contains two contiguous branches. Tensor-force coupling uses the lower/upper
	// orbital branch (L<J or L>J), while isospin-breaking in 1S0 uses the T_3N branch label. / [CN] 耦合 SWP 存储里
	// 连续放着两条分支：张量力耦合时按轨道角动量分支选取（L<J 或 L>J），而 1S0 的同位旋破缺则按 T_3N 分支选取。
	if (coupled_via_L_2N){
		if (L_2N>J_2N){
			channel_branch_ptr += Np_WP + 1;
		}
		return channel_branch_ptr;
	}

	if (two_T_3N!=1){
		channel_branch_ptr += Np_WP + 1;
	}
	return channel_branch_ptr;
}

} // namespace
