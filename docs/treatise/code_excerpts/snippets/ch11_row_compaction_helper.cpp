// ===============================================================
// 抽取自仓库 [current]: src/core/faddeev_solver/solve_faddeev.cpp
// 行号区段：86..122
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
// [EN] Row compaction is purely an acceleration device: skip rows whose observed on-shell elements have already
// converged, but leave the active rows untouched. / [CN] 行压缩纯粹是加速手段：跳过那些其可观测 on-shell 元素已经收敛的行，
// 但对仍然活动的行不做任何代数近似。
bool row_has_only_converged_targets(size_t idx_d_row,
									size_t idx_q_com,
									size_t num_deuteron_states,
									size_t num_q_com,
									const bool* pade_approximants_conv_array,
									const bool* pade_approximants_BU_conv_array,
									const channel_os_indexing& chn_os_indexing,
									bool include_breakup_channels){
	for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
		size_t idx_NDOS = elastic_value_storage_index(idx_d_row,
														 idx_d_col,
														 idx_q_com,
														 num_deuteron_states,
														 num_q_com);
		if (pade_approximants_conv_array[idx_NDOS]==false){
			return false;
		}
	}

	if (include_breakup_channels){
		int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[idx_q_com];
		int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[idx_q_com+1];
		for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){
			size_t idx_NDOS = breakup_value_storage_index(idx_d_row,
														 idx_BU_chn,
														 chn_os_indexing.num_BU_chns);
			if (pade_approximants_BU_conv_array[idx_NDOS]==false){
				return false;
			}
		}
	}
	return true;
}

