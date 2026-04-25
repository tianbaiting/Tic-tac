// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_swp_states.cpp
// 行号区段：4..56
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
namespace {

// [EN] The deuteron is the unique bound subsystem carried by the triplet S-wave channel. / [CN] 氘核是三重态
// S 波通道携带的唯一束缚子系统。
bool is_triplet_s_wave_channel(int L_2N, int S_2N, int J_2N, int T_2N){
	return L_2N==0 && S_2N==1 && J_2N==1 && T_2N==0;
}

// [EN] 1S0 needs special bookkeeping because optional isospin breaking reuses the coupled-channel storage pattern.
// / [CN] 1S0 需要特殊 bookkeeping，因为可选的同位旋破缺会复用耦合通道的存储模式。
bool is_singlet_s_wave_channel(int L_2N, int S_2N, int J_2N){
	return L_2N==0 && S_2N==0 && J_2N==0;
}

// [EN] Decide whether a two-body block occupies the uncoupled Np x Np storage or the coupled 2Np x 2Np storage.
// / [CN] 判断某个两体块应放入非耦合的 Np x Np 存储，还是耦合的 2Np x 2Np 存储。
bool uses_coupled_storage(int L_2N_row,
						  int L_2N_col,
						  int S_2N,
						  int J_2N,
						  const run_params& run_parameters){
	const bool coupled_via_L_2N = run_parameters.tensor_force && (L_2N_row!=L_2N_col || (L_2N_row==L_2N_col && L_2N_row!=J_2N && J_2N!=0));
	const bool coupled_via_T_3N = is_singlet_s_wave_channel(L_2N_row, S_2N, J_2N) && run_parameters.isospin_breaking_1S0;
	if (coupled_via_L_2N && coupled_via_T_3N){
		raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
	}
	return coupled_via_L_2N || coupled_via_T_3N;
}

// [EN] Claim one distinct 2N Hamiltonian block the first time it is encountered in the enclosing 3N basis loop.
// / [CN] 在外层三体基循环里首次遇到某个不同的 2N Hamiltonian 块时，将其标记为“已领取”。
bool claim_channel(std::vector<bool>& channel_done_flags, int channel_index){
	if (channel_done_flags[channel_index]){
		return false;
	}
	channel_done_flags[channel_index] = true;
	return true;
}

// [EN] The free packet Hamiltonian is identical for every uncoupled branch and for each leg of a coupled branch, so
// we build the repeated H0 tables once here. / [CN] 自由波包 Hamiltonian 对所有非耦合分支以及每条耦合分支腿都是相同的，
// 因此在这里一次性构建并复用这些重复的 H0 表。
void fill_free_hamiltonian_branches(std::vector<double>& free_hamiltonian_branches,
									int num_branches,
									int Np_WP,
									double* p_WP_array){
	for (int idx_branch=0; idx_branch<num_branches; idx_branch++){
		construct_free_hamiltonian(&free_hamiltonian_branches[idx_branch * Np_WP],
								   Np_WP,
								   p_WP_array);
	}
}

