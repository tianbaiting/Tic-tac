#include "solve_faddeev.h"
#include "interactions/three_nucleon_force_model.h"
#include "interactions/w1_pw_cache.h"
#include <gsl/gsl_cblas.h>
#include <cstring>
#include <vector>

namespace {

// [EN] The deuteron packet always occupies p-index 0 after the SWP eigenvalue ordering used in this solver, so all
// elastic on-shell amplitudes live on one fixed p slice. / [CN] 按当前求解器使用的 SWP 本征值排序，氘核波包总是位于
// p-index 0，因此所有弹性 on-shell 振幅都落在同一个固定的 p 切片上。
constexpr size_t elastic_bound_packet_index = 0;

struct elastic_on_shell_index {
	size_t alpha_row;
	size_t alpha_col;
	size_t q_idx;
	size_t row_storage_idx;
	size_t col_storage_idx;
	size_t value_storage_idx;
};

// [EN] Flatten one (alpha, q, p) packet state into the dense packet-lattice index used everywhere in the kernel.
// / [CN] 把单个 (alpha, q, p) 波包态展平成核内部统一使用的稠密格点索引。
inline size_t dense_packet_index(size_t alpha_idx,
								 size_t q_idx,
								 size_t p_idx,
								 size_t Nq_WP,
								 size_t Np_WP){
	return alpha_idx*Nq_WP*Np_WP + q_idx*Np_WP + p_idx;
}

// [EN] Elastic rows are stored as (deuteron-row, q-shell) blocks because that is the minimum subset needed by the
// Miller/Sean workflow before Padé resummation. / [CN] 弹性行按 (deuteron-row, q-shell) 分块存储，因为在 Sean Miller
// 的工作流里，Padé 重求和前只需要这部分最小子集。
inline size_t elastic_row_storage_index(size_t idx_d_row,
										size_t idx_q_com,
										size_t num_q_com){
	return idx_d_row*num_q_com + idx_q_com;
}

inline size_t elastic_value_storage_index(size_t idx_d_row,
										  size_t idx_d_col,
										  size_t idx_q_com,
										  size_t num_deuteron_states,
										  size_t num_q_com){
	return idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + idx_q_com;
}

inline size_t breakup_value_storage_index(size_t idx_d_row,
										  size_t idx_BU_chn,
										  size_t num_BU_chns){
	return idx_d_row*num_BU_chns + idx_BU_chn;
}

// [EN] Bundle the repeated elastic on-shell bookkeeping so the hot loops can read more like the multiple-scattering
// formulas in the paper and less like raw index arithmetic. / [CN] 把重复出现的弹性 on-shell bookkeeping 打包起来，
// 这样热点循环读起来会更像论文里的多重散射公式，而不是裸索引运算。
inline elastic_on_shell_index make_elastic_on_shell_index(size_t idx_d_row,
														  size_t idx_d_col,
														  size_t idx_q_com,
														  const int* deuteron_idx_array,
														  const int* q_com_idx_array,
														  size_t num_deuteron_states,
														  size_t num_q_com,
														  size_t Nq_WP,
														  size_t Np_WP){
	elastic_on_shell_index index = {};
	index.alpha_row = deuteron_idx_array[idx_d_row];
	index.alpha_col = deuteron_idx_array[idx_d_col];
	index.q_idx = q_com_idx_array[idx_q_com];
	index.row_storage_idx = elastic_row_storage_index(idx_d_row, idx_q_com, num_q_com);
	index.col_storage_idx = dense_packet_index(index.alpha_col,
												 index.q_idx,
												 elastic_bound_packet_index,
												 Nq_WP,
												 Np_WP);
	index.value_storage_idx = elastic_value_storage_index(idx_d_row,
															idx_d_col,
															idx_q_com,
															num_deuteron_states,
															num_q_com);
	return index;
}

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

} // namespace

// Helper function for in-place matrix transpose
void inplace_transpose(double* A, int rows, int cols) {
    if (rows == cols) { // Square matrix
        for (int i = 0; i < rows; ++i) {
            for (int j = i + 1; j < cols; ++j) {
                double temp = A[i * cols + j];
                A[i * cols + j] = A[j * cols + i];
                A[j * cols + i] = temp;
            }
        }
    } else { // Non-square matrix, requires a temporary buffer
        std::vector<double> B(rows * cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                B[j * rows + i] = A[i * cols + j];
            }
        }
        memcpy(A, B.data(), rows * cols * sizeof(double));
    }
}

// [EN] This is the raw driving column PVC that appears in the AGS kernel before the left basis rotation by C^T.
// Applying the sparse permutation first keeps the expensive part of the kernel sparse. / [CN] 这里计算的是 AGS 核中
// 左乘 C^T 之前的原始驱动列 PVC；先施加稀疏置换算符，可以把核中最昂贵的部分保持为稀疏结构。
void calculate_PVC_col(double*  col_array,
					   size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
					   size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
					   double** VC_CM_array,
					   double*  P123_val_array,
					   int*  	P123_row_array,
					   size_t*  P123_col_array,
					   size_t   P123_dim){

	double* VC_subarray     = NULL;

	size_t dense_dim = Nalpha*Np_WP*Nq_WP;

	for (size_t idx_alpha_j=0; idx_alpha_j<Nalpha; idx_alpha_j++){
		VC_subarray = VC_CM_array[idx_alpha_c*Nalpha + idx_alpha_j];

		/* Only do inner-product if VC is not zero due to conservation laws */
		if (VC_subarray!=NULL){
			for (size_t idx_p_j=0; idx_p_j<Np_WP; idx_p_j++){

				/* Access VC element */
				double VC_element = VC_subarray[idx_p_c*Np_WP + idx_p_j];

				/* Inner-product index */
				size_t idx_j =  idx_alpha_j*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_j;

				/* CSC-format indexing */
				size_t idx_i_lower = P123_col_array[idx_j    ];
				size_t idx_i_upper = P123_col_array[idx_j + 1];

				/* Loop through rows of column we're calculating, and append */
				for (size_t idx_i=idx_i_lower; idx_i<idx_i_upper; idx_i++){

					/* Access arrays, this whole function is written to minimize these two calls */
					/* NOTE THAT P = P123 + P132 = 2*P123 FOR ANTISYMMETRIC PAIR-STATES */
					double P_element  = 2*P123_val_array[idx_i];

					size_t idx_row = P123_row_array[idx_i];
					/* Write inner product to col_array of PVC-product */
					col_array[idx_row] += P_element * VC_element;
				}
			}
		}
	}
}

// [EN] CPVC = C^T P V C is the packet-space kernel that drives both the first Neumann term and every later
// rescattering step. We form one dense column at a time because that matches both the sparse P123 access pattern
// and the later chunked GEMM strategy. / [CN] CPVC = C^T P V C 是波包空间中的核，它既驱动第一项 Neumann 项，
// 也驱动之后所有再散射步骤。这里按“每次一列”构造，是为了同时匹配稀疏 P123 的访问模式和后面分块 GEMM 的策略。
void calculate_CPVC_col(double*  col_array,
					    int* 	 row_to_nnz_array,
					    int* 	 nnz_to_row_array,
					    size_t&  num_nnz,
					    size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
					    size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
					    double** CT_RM_array,
					    double** VC_CM_array,
					    double*  P123_val_array,
					    int*     P123_row_array,
					    size_t*  P123_col_array,
					    size_t   P123_dim,
					    const tnf_kernel_context& tnf_ctx){

	/* Generate PVC-column */
	std::vector<double> PVC_col(Nalpha*Nq_WP*Np_WP);
	/* Ensure PVC_col contains only zeroes */
	for (size_t idx=0; idx<Nalpha*Nq_WP*Np_WP; idx++){
		PVC_col[idx] = 0;
	}

	//auto timestamp_start = std::chrono::system_clock::now();
	calculate_PVC_col(PVC_col.data(),
					  idx_alpha_c, idx_p_c, idx_q_c,
					  Nalpha,      Nq_WP,   Np_WP,
					  VC_CM_array,
					  P123_val_array,
					  P123_row_array,
					  P123_col_array,
					  P123_dim);
	//auto timestamp_end = std::chrono::system_clock::now();
	//std::chrono::duration<double>  time1 = timestamp_end - timestamp_start;
	//printf("TIME PVC:  %.6f \n", time1.count()); fflush(stdout);

	// [EN] 3NF contribution (Born level): add W^(1)·(1+P)·C to PVC_col buffer.
	//   Identity part:    W^(1)·C  — direct W1×C sum over intermediate p
	//   Permutation part: W^(1)·P·C — compute P·C column via sparse P123 loop, then apply W^(1)
	// / [CN] 3NF 贡献（Born 级）：把 W^(1)·(1+P)·C 加到 PVC_col 缓冲区。
	//   单位部分：W^(1)·C —— 对中间 p 直接做 W1×C 求和
	//   置换部分：W^(1)·P·C —— 用稀疏 P123 循环计算 P·C 列，然后施加 W^(1)
	if (tnf_ctx.tnf != nullptr && tnf_ctx.tnf->enabled() && tnf_ctx.w1_scale != 0.0){
		const three_nucleon_force_model* tnf = tnf_ctx.tnf;
		const pw_3N_statespace& pw_st = *tnf_ctx.pw_states;
		const double* p_WP = tnf_ctx.p_WP_array;
		const double* q_WP = tnf_ctx.q_WP_array;

		// [EN] WP bin-averaging normalization for the 3NF matrix element W^(1)_WP.
		// V_WP uses p_r × p_c × √dp_r × √dp_c (one momentum factor + √(bin width) per side).
		// P123_WP uses 1/(√dp_r × √dq_r × √dp_c × √dq_c) from WP normalization.
		// For W^(1), which acts in the full (p,q) space, the WP matrix element is:
		//   W^(1)_WP(p_r,q_r; p_c,q_c) = p_r × q_r × p_c × q_c × √dp_r × √dq_r × √dp_c × √dq_c × W^(1)(mids)
		// This follows the reduced-function convention g(p,q) = p×q×ψ(p,q) consistent with V_WP.
		// / [CN] 3NF 矩阵元 W^(1)_WP 的 WP 基平均归一化。
		// V_WP 使用 p_r × p_c × √dp_r × √dp_c（每侧一个动量因子 + √(bin 宽)）。
		// P123_WP 使用 1/(√dp_r × √dq_r × √dp_c × √dq_c) 作为 WP 归一化。
		// 对于作用在完整 (p,q) 空间的 W^(1)，其 WP 矩阵元为：
		//   W^(1)_WP = p_r × q_r × p_c × q_c × √dp_r × √dq_r × √dp_c × √dq_c × W^(1)(mids)

		// [EN] Momentum unit convention: p_WP/q_WP are stored in MeV (consistent with 2NF V, which uses MeV),
		// but the 3NF W1_element API expects Jacobi momenta in fm^{-1}. We convert on the fly and scale
		// the output so the WP matrix element lives in the same MeV-based convention as V_WP.
		// W^(1) returned in fm^5 has natural-unit equivalence (hbarc)^5 MeV^{-5}; multiplying by 1/hbarc^5
		// converts it to MeV^{-5}, which then combines with the four momentum×sqrt(bin-width) factors
		// (MeV^6) to give a W1_WP matrix element in MeV — matching how V_WP (= p p' √(dp dp') V_MeV)
		// combines to give MeV^4 and then integrates against MeV-measure basis functions.
		// / [CN] p_WP/q_WP 为 MeV（与 2NF V 一致），但 W1_element API 约定动量单位为 fm^{-1}，
		// 故在调用处做 MeV→fm^{-1} 换算，并用 1/hbarc^5 把 fm^5 的输出换成 MeV^{-5}，与 V_WP 的
		// MeV 基约定一致。
		const double inv_hbarc  = 1.0 / hbarc;
		const double inv_hbarc5 = inv_hbarc * inv_hbarc * inv_hbarc * inv_hbarc * inv_hbarc;
		const double w1_unit    = inv_hbarc5 * tnf_ctx.w1_scale;  // combined unit-conversion + diagnostic knob

		// q_c bin midpoint and WP normalization factor
		double q_c_mid = 0.5 * (q_WP[idx_q_c] + q_WP[idx_q_c + 1]);
		double dq_c    = q_WP[idx_q_c + 1] - q_WP[idx_q_c];
		double wq_c    = q_c_mid * std::sqrt(dq_c);  // q_c × √dq_c  [MeV^{3/2}]
		double q_c_fm  = q_c_mid * inv_hbarc;

		double* C_block = CT_RM_array[idx_alpha_c * Nalpha + idx_alpha_c];

		if (C_block != nullptr){
			for (size_t idx_alpha_r = 0; idx_alpha_r < Nalpha; idx_alpha_r++){
				if (pw_st.two_J_3N_array[idx_alpha_r] != pw_st.two_J_3N_array[idx_alpha_c]) continue;
				if (pw_st.two_T_3N_array[idx_alpha_r] != pw_st.two_T_3N_array[idx_alpha_c]) continue;
				if (pw_st.P_3N_array[idx_alpha_r]     != pw_st.P_3N_array[idx_alpha_c])     continue;

				for (size_t idx_q_r = 0; idx_q_r < Nq_WP; idx_q_r++){
					double q_r_mid = 0.5 * (q_WP[idx_q_r] + q_WP[idx_q_r + 1]);
					double dq_r    = q_WP[idx_q_r + 1] - q_WP[idx_q_r];
					double wq_r    = q_r_mid * std::sqrt(dq_r);
					double q_r_fm  = q_r_mid * inv_hbarc;

					for (size_t idx_p_r = 0; idx_p_r < Np_WP; idx_p_r++){
						double p_r_mid = 0.5 * (p_WP[idx_p_r] + p_WP[idx_p_r + 1]);
						double dp_r    = p_WP[idx_p_r + 1] - p_WP[idx_p_r];
						double wp_r    = p_r_mid * std::sqrt(dp_r);
						double p_r_fm  = p_r_mid * inv_hbarc;

						// Sum over intermediate p_j (WP basis) with C weight.
						// C_block = CT_RM_array[alpha_c, alpha_c] stores C^T, so C[p_j, p_c] = C^T[p_c, p_j]
						// = C_block[idx_p_c * Np_WP + idx_p_j] — NOT C_block[idx_p_j * Np_WP + idx_p_c].
						double w1c_element = 0.0;
						for (size_t idx_p_j = 0; idx_p_j < Np_WP; idx_p_j++){
							double C_val = C_block[idx_p_c * Np_WP + idx_p_j];
							if (C_val == 0.0) continue;

							double p_j_mid = 0.5 * (p_WP[idx_p_j] + p_WP[idx_p_j + 1]);
							double dp_j    = p_WP[idx_p_j + 1] - p_WP[idx_p_j];
							double wp_j    = p_j_mid * std::sqrt(dp_j);
							double p_j_fm  = p_j_mid * inv_hbarc;

							// Cache value is the WP bin matrix element in MeV (pre-w1_scale);
							// fallback does the legacy 1-point midpoint computation inline.
							double w1_bin;
							if (tnf_ctx.w1_cache != nullptr) {
								w1_bin = tnf_ctx.w1_cache->get(idx_alpha_r, idx_alpha_c,
															   idx_p_r, idx_q_r, idx_p_j, idx_q_c)
								         * tnf_ctx.w1_scale;
							} else {
								double w1_raw = tnf->W1_element(idx_alpha_r, idx_alpha_c,
																 p_r_fm, q_r_fm,
																 p_j_fm, q_c_fm,
																 pw_st);
								w1_bin = (w1_raw * w1_unit) * (wp_r * wq_r * wp_j * wq_c);
							}
							w1c_element += w1_bin * C_val;
						}

						if (w1c_element != 0.0){
							size_t idx_row = idx_alpha_r * Nq_WP * Np_WP + idx_q_r * Np_WP + idx_p_r;
							PVC_col[idx_row] += w1c_element;
						}
					}
				}
			}
		}

		// [EN] W^(1)·P·C permutation contribution.
		// / [CN] W^(1)·P·C 置换贡献。
		std::vector<double> PC_col(Nalpha * Nq_WP * Np_WP, 0.0);
		calculate_PVC_col(PC_col.data(),
						  idx_alpha_c, idx_p_c, idx_q_c,
						  Nalpha, Nq_WP, Np_WP,
						  tnf_ctx.CT_RM_array,
						  P123_val_array,
						  P123_row_array,
						  P123_col_array,
						  P123_dim);

		// Apply W^(1)_WP to PC_col. PC_col elements are already in WP basis (from P·C).
		// W^(1)_WP is a matrix in WP basis, so this is a standard matrix-vector product:
		//   (W1_WP × PC_col)[row] = Σ_k W1_WP[row,k] × PC_col[k]
		// W1_WP[row,k] = p_r q_r p_k q_k √dp_r √dq_r √dp_k √dq_k × W1(mids) × 1/hbarc^5
		// (MeV→fm^{-1} for the W1 call; fm^5 → MeV^{-5} via 1/hbarc^5.)
		for (size_t idx_alpha_k = 0; idx_alpha_k < Nalpha; idx_alpha_k++){
			for (size_t idx_q_k = 0; idx_q_k < Nq_WP; idx_q_k++){
				for (size_t idx_p_k = 0; idx_p_k < Np_WP; idx_p_k++){
					size_t idx_k = idx_alpha_k * Nq_WP * Np_WP + idx_q_k * Np_WP + idx_p_k;
					double pc_val = PC_col[idx_k];
					if (pc_val == 0.0) continue;

					double p_k_mid = 0.5 * (p_WP[idx_p_k] + p_WP[idx_p_k + 1]);
					double q_k_mid = 0.5 * (q_WP[idx_q_k] + q_WP[idx_q_k + 1]);
					double dp_k    = p_WP[idx_p_k + 1] - p_WP[idx_p_k];
					double dq_k    = q_WP[idx_q_k + 1] - q_WP[idx_q_k];
					double wp_k    = p_k_mid * std::sqrt(dp_k);
					double wq_k    = q_k_mid * std::sqrt(dq_k);
					double p_k_fm  = p_k_mid * inv_hbarc;
					double q_k_fm  = q_k_mid * inv_hbarc;

					for (size_t idx_alpha_r = 0; idx_alpha_r < Nalpha; idx_alpha_r++){
						if (pw_st.two_J_3N_array[idx_alpha_r] != pw_st.two_J_3N_array[idx_alpha_k]) continue;
						if (pw_st.two_T_3N_array[idx_alpha_r] != pw_st.two_T_3N_array[idx_alpha_k]) continue;
						if (pw_st.P_3N_array[idx_alpha_r]     != pw_st.P_3N_array[idx_alpha_k])     continue;

						for (size_t idx_q_r = 0; idx_q_r < Nq_WP; idx_q_r++){
							double q_r_mid = 0.5 * (q_WP[idx_q_r] + q_WP[idx_q_r + 1]);
							double dq_r    = q_WP[idx_q_r + 1] - q_WP[idx_q_r];
							double wq_r    = q_r_mid * std::sqrt(dq_r);
							double q_r_fm  = q_r_mid * inv_hbarc;
							for (size_t idx_p_r = 0; idx_p_r < Np_WP; idx_p_r++){
								double p_r_mid = 0.5 * (p_WP[idx_p_r] + p_WP[idx_p_r + 1]);
								double dp_r    = p_WP[idx_p_r + 1] - p_WP[idx_p_r];
								double wp_r    = p_r_mid * std::sqrt(dp_r);
								double p_r_fm  = p_r_mid * inv_hbarc;

								// Cache value is the WP bin matrix element in MeV (pre-w1_scale);
								// fallback does the legacy 1-point midpoint computation inline.
								double w1_wp;
								if (tnf_ctx.w1_cache != nullptr) {
									w1_wp = tnf_ctx.w1_cache->get(idx_alpha_r, idx_alpha_k,
																  idx_p_r, idx_q_r, idx_p_k, idx_q_k)
									        * tnf_ctx.w1_scale;
								} else {
									double w1_raw = tnf->W1_element(idx_alpha_r, idx_alpha_k,
																	 p_r_fm, q_r_fm,
																	 p_k_fm, q_k_fm, pw_st);
									w1_wp = (w1_raw * w1_unit) * (wp_r * wq_r * wp_k * wq_k);
								}
								if (w1_wp != 0.0){
									size_t idx_row = idx_alpha_r * Nq_WP * Np_WP + idx_q_r * Np_WP + idx_p_r;
									PVC_col[idx_row] += w1_wp * pc_val;
								}
							}
						}
					}
				}
			}
		}
	}

	/* THOUGHT:
	 * MOVE ALPHA_I OUTWARDS AND GO BACK TO DIRECT APPEND TO COL_ARRAY.
	 * BUT USE NON-ZERO IF-TESTING. */

	/* TOUGHT (CONTRADICTORY TO DIRECT APPEND):
	 * Use dense mat-vec multiplication for sub-blocks  (Can let q be columns of right-vectors? Appealing use of MM-multiplication)
	 * Use dense vec-vec (or vec-mat-vec?) multiplication for A_An */

	/* Generate (C^T x PVC)-column */
	double* CT_subarray     = NULL;
	double* CT_subarray_row = NULL;
	double* PVC_subcol 		= NULL;

	
	/*  Looping based on first nnz-lookup of CT-matrices */
	/* Loop over rows of col_array */
	//timestamp_start = std::chrono::system_clock::now();
	for (auto idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		/* Beginning of inner-product loops (index "i") */
		for (auto idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
			auto idx_CT_2N_block = idx_alpha_r*Nalpha + idx_alpha_i;
			CT_subarray = CT_RM_array[idx_CT_2N_block];
			/* Only do inner-product if CT is not zero due to conservation laws */
			if (CT_subarray!=NULL){
				PVC_subcol = &PVC_col[idx_alpha_i*Nq_WP*Np_WP];
				for (auto idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
					for (auto idx_p_i=0; idx_p_i<Np_WP; idx_p_i++){
						//size_t idx_PVC     = idx_alpha_i*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_i;
						auto PVC_element = PVC_subcol[idx_q_r*Np_WP + idx_p_i];
						if (PVC_element!=0){
							for (auto idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){
								/* I'm not sure if this is the fastest ordering of the loops */
								//double CT_element  = CT_subarray[idx_p_r*Np_WP + idx_p_i];
								//auto prod = PVC_element  * CT_subarray[idx_p_r*Np_WP + idx_p_i];
								//if (prod!=0){
									auto idx_CPVC = idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_r;
									col_array[idx_CPVC] += PVC_element  * CT_subarray[idx_p_r*Np_WP + idx_p_i];
								//}
							}
						}
					}
				}
			}
		}
	}
	
	//timestamp_end = std::chrono::system_clock::now();
	//std::chrono::duration<double>  time2 = timestamp_end - timestamp_start;
	//printf("TIME CPVC:  %.6f \n", time2.count()); fflush(stdout);
}

// [EN] Only on-shell rows are ever observed in the elastic output, so the initial driving term A is assembled only
// for those rows. This is the first major WPCD reduction compared with a full dense solve. / [CN] 最终弹性输出只会用到
// on-shell 行，因此初始驱动项 A 也只为这些行构造；这是相对于完整稠密求解的第一层 WPCD 降维。
void calculate_all_CPVC_rows(double*  row_arrays,
							 int*	  q_com_idx_array,	  size_t num_q_com,
					   		 int*     deuteron_idx_array, size_t num_deuteron_states,
							 size_t   Nalpha,
							 size_t   Nq_WP,
							 size_t   Np_WP,
							 double** CT_RM_array,
							 double** VC_CM_array,
							 double*  P123_val_array,
							 int*     P123_row_array,
							 size_t*  P123_col_array,
							 size_t   P123_dim,
							 const tnf_kernel_context& tnf_ctx){
	
	size_t dense_dim = Nalpha*Nq_WP*Np_WP;

	/* Loop over cols of row */
	#pragma omp parallel
	{
	/* Generate PVC-column */
	std::vector<double> PVC_col(dense_dim);
	/* Generate (C^T x PVC)-column */
	double* CT_subarray     = NULL;
	double* CT_subarray_row = NULL;
	#pragma omp for
	for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
		for (size_t idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
			for (size_t idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){

				/* Ensure PVC_col contains only zeroes */
				for (size_t idx=0; idx<dense_dim; idx++){
					PVC_col[idx] = 0;
				}
				
				/* Calculate PVC-column for alpha_i, p_i, q_r */
				calculate_PVC_col(PVC_col.data(),
								  idx_alpha_c, idx_p_c, idx_q_c,
								  Nalpha,      Nq_WP,   Np_WP,
								  VC_CM_array,
								  P123_val_array,
								  P123_row_array,
								  P123_col_array,
								  P123_dim);

				// [EN] 3NF contributions: W^(1)·(1+P)·C = W^(1)·C + W^(1)·P·C
				// Momenta converted MeV→fm^{-1} for W1 API; output scaled fm^5→MeV^{-5}.
				// w1_scale==0 short-circuit skips the ~O(Nα² Nq² Np³) W1_element hot path.
				if (tnf_ctx.tnf != nullptr && tnf_ctx.tnf->enabled() && tnf_ctx.w1_scale != 0.0){
					const three_nucleon_force_model* tnf = tnf_ctx.tnf;
					const pw_3N_statespace& pw_st = *tnf_ctx.pw_states;
					const double* p_WP = tnf_ctx.p_WP_array;
					const double* q_WP = tnf_ctx.q_WP_array;

					const double inv_hbarc  = 1.0 / hbarc;
					const double inv_hbarc5 = inv_hbarc * inv_hbarc * inv_hbarc * inv_hbarc * inv_hbarc;
					const double w1_unit    = inv_hbarc5 * tnf_ctx.w1_scale;

					double q_c_mid = 0.5 * (q_WP[idx_q_c] + q_WP[idx_q_c + 1]);
					double dq_c    = q_WP[idx_q_c + 1] - q_WP[idx_q_c];
					double wq_c    = q_c_mid * std::sqrt(dq_c);
					double q_c_fm  = q_c_mid * inv_hbarc;

					// Identity part: W^(1)·C
					double* C_block = tnf_ctx.CT_RM_array[idx_alpha_c * Nalpha + idx_alpha_c];
					if (C_block != nullptr){
						for (size_t idx_alpha_r = 0; idx_alpha_r < Nalpha; idx_alpha_r++){
							if (pw_st.two_J_3N_array[idx_alpha_r] != pw_st.two_J_3N_array[idx_alpha_c]) continue;
							if (pw_st.two_T_3N_array[idx_alpha_r] != pw_st.two_T_3N_array[idx_alpha_c]) continue;
							if (pw_st.P_3N_array[idx_alpha_r]     != pw_st.P_3N_array[idx_alpha_c])     continue;
							for (size_t idx_q_r = 0; idx_q_r < Nq_WP; idx_q_r++){
								double q_r_mid = 0.5 * (q_WP[idx_q_r] + q_WP[idx_q_r + 1]);
								double dq_r    = q_WP[idx_q_r + 1] - q_WP[idx_q_r];
								double wq_r    = q_r_mid * std::sqrt(dq_r);
								double q_r_fm  = q_r_mid * inv_hbarc;
								for (size_t idx_p_r = 0; idx_p_r < Np_WP; idx_p_r++){
									double p_r_mid = 0.5 * (p_WP[idx_p_r] + p_WP[idx_p_r + 1]);
									double dp_r    = p_WP[idx_p_r + 1] - p_WP[idx_p_r];
									double wp_r    = p_r_mid * std::sqrt(dp_r);
									double p_r_fm  = p_r_mid * inv_hbarc;
									double w1c = 0.0;
									for (size_t idx_p_j = 0; idx_p_j < Np_WP; idx_p_j++){
										// CT_RM_array stores C^T in row-major; C[p_j, p_c] = C^T[p_c, p_j] = C_block[p_c*Np+p_j]
										double C_val = C_block[idx_p_c * Np_WP + idx_p_j];
										if (C_val == 0.0) continue;
										double p_j_mid = 0.5 * (p_WP[idx_p_j] + p_WP[idx_p_j + 1]);
										double dp_j    = p_WP[idx_p_j + 1] - p_WP[idx_p_j];
										double wp_j    = p_j_mid * std::sqrt(dp_j);
										double p_j_fm  = p_j_mid * inv_hbarc;
										// Cache value is the WP bin matrix element in MeV (pre-w1_scale);
										// fallback does the legacy 1-point midpoint computation inline.
										double w1_bin;
										if (tnf_ctx.w1_cache != nullptr) {
											w1_bin = tnf_ctx.w1_cache->get(idx_alpha_r, idx_alpha_c,
																		   idx_p_r, idx_q_r, idx_p_j, idx_q_c)
											         * tnf_ctx.w1_scale;
										} else {
											double w1_raw = tnf->W1_element(idx_alpha_r, idx_alpha_c,
																			 p_r_fm, q_r_fm,
																			 p_j_fm, q_c_fm, pw_st);
											w1_bin = (w1_raw * w1_unit) * (wp_r * wq_r * wp_j * wq_c);
										}
										w1c += w1_bin * C_val;
									}
									if (w1c != 0.0){
										PVC_col[idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_r] += w1c;
									}
								}
							}
						}
					}

					// Permutation part: W^(1)·P·C
					std::vector<double> PC_col_local(Nalpha * Nq_WP * Np_WP, 0.0);
					calculate_PVC_col(PC_col_local.data(),
									  idx_alpha_c, idx_p_c, idx_q_c,
									  Nalpha, Nq_WP, Np_WP,
									  tnf_ctx.CT_RM_array,
									  P123_val_array,
									  P123_row_array,
									  P123_col_array,
									  P123_dim);
					for (size_t idx_k = 0; idx_k < Nalpha * Nq_WP * Np_WP; idx_k++){
						double pc_val = PC_col_local[idx_k];
						if (pc_val == 0.0) continue;
						size_t idx_alpha_k = idx_k / (Nq_WP * Np_WP);
						size_t idx_q_k = (idx_k % (Nq_WP * Np_WP)) / Np_WP;
						size_t idx_p_k = idx_k % Np_WP;
						double p_k_mid = 0.5 * (p_WP[idx_p_k] + p_WP[idx_p_k + 1]);
						double q_k_mid = 0.5 * (q_WP[idx_q_k] + q_WP[idx_q_k + 1]);
						double dp_k    = p_WP[idx_p_k + 1] - p_WP[idx_p_k];
						double dq_k    = q_WP[idx_q_k + 1] - q_WP[idx_q_k];
						double wp_k    = p_k_mid * std::sqrt(dp_k);
						double wq_k    = q_k_mid * std::sqrt(dq_k);
						double p_k_fm  = p_k_mid * inv_hbarc;
						double q_k_fm  = q_k_mid * inv_hbarc;
						for (size_t idx_alpha_r = 0; idx_alpha_r < Nalpha; idx_alpha_r++){
							if (pw_st.two_J_3N_array[idx_alpha_r] != pw_st.two_J_3N_array[idx_alpha_k]) continue;
							if (pw_st.two_T_3N_array[idx_alpha_r] != pw_st.two_T_3N_array[idx_alpha_k]) continue;
							if (pw_st.P_3N_array[idx_alpha_r]     != pw_st.P_3N_array[idx_alpha_k])     continue;
							for (size_t idx_q_r = 0; idx_q_r < Nq_WP; idx_q_r++){
								double q_r_mid = 0.5 * (q_WP[idx_q_r] + q_WP[idx_q_r + 1]);
								double dq_r    = q_WP[idx_q_r + 1] - q_WP[idx_q_r];
								double wq_r    = q_r_mid * std::sqrt(dq_r);
								double q_r_fm  = q_r_mid * inv_hbarc;
								for (size_t idx_p_r = 0; idx_p_r < Np_WP; idx_p_r++){
									double p_r_mid = 0.5 * (p_WP[idx_p_r] + p_WP[idx_p_r + 1]);
									double dp_r    = p_WP[idx_p_r + 1] - p_WP[idx_p_r];
									double wp_r    = p_r_mid * std::sqrt(dp_r);
									double p_r_fm  = p_r_mid * inv_hbarc;
									// Cache value is the WP bin matrix element in MeV (pre-w1_scale);
									// fallback does the legacy 1-point midpoint computation inline.
									double w1_wp;
									if (tnf_ctx.w1_cache != nullptr) {
										w1_wp = tnf_ctx.w1_cache->get(idx_alpha_r, idx_alpha_k,
																	  idx_p_r, idx_q_r, idx_p_k, idx_q_k)
										        * tnf_ctx.w1_scale;
									} else {
										double w1_raw = tnf->W1_element(idx_alpha_r, idx_alpha_k,
																		 p_r_fm, q_r_fm,
																		 p_k_fm, q_k_fm, pw_st);
										w1_wp = (w1_raw * w1_unit) * (wp_r * wq_r * wp_k * wq_k);
									}
									if (w1_wp != 0.0){
										PVC_col[idx_alpha_r * Nq_WP * Np_WP + idx_q_r * Np_WP + idx_p_r] += w1_wp * pc_val;
									}
								}
							}
						}
					}
				}

				/* Re-use PVC-column in all relevant calculations */
				for (size_t i=0; i<num_deuteron_states; i++){
					for (size_t j=0; j<num_q_com; j++){
						const elastic_on_shell_index ndos = make_elastic_on_shell_index(i,
																						 0,
																						 j,
																						 deuteron_idx_array,
																						 q_com_idx_array,
																						 num_deuteron_states,
																						 num_q_com,
																						 Nq_WP,
																						 Np_WP);

						double inner_product_CPVC = 0;
						/* Beginning of inner-product loops (index "i") */
						for (size_t idx_alpha_i=0; idx_alpha_i<Nalpha; idx_alpha_i++){
							
							CT_subarray = CT_RM_array[ndos.alpha_row*Nalpha + idx_alpha_i];
		
							/* Only do inner-product if CT is not zero due to conservation laws */
							if (CT_subarray!=NULL){
								CT_subarray_row = &CT_subarray[elastic_bound_packet_index*Np_WP];
		
								for (size_t idx_p_i=0; idx_p_i<Np_WP; idx_p_i++){
								
									size_t idx_PVC     = dense_packet_index(idx_alpha_i,
																			ndos.q_idx,
																			idx_p_i,
																			Nq_WP,
																			Np_WP);
									double PVC_element = PVC_col[idx_PVC];
		
									/* I'm not sure if this is the fastest ordering of the loops */
									double CT_element  = CT_subarray_row[idx_p_i];
		
									inner_product_CPVC += CT_element * PVC_element;
								}
							}
						}
				
						size_t idx_CPVC = dense_packet_index(idx_alpha_c, idx_q_c, idx_p_c, Nq_WP, Np_WP);
						row_arrays[ndos.row_storage_idx*dense_dim + idx_CPVC] = inner_product_CPVC;
					}
				}
			}
		}
	}
	}
}

cdouble pade_approximant(cdouble* a_coeff_array, size_t N, size_t M, cdouble z){

	/* a_coeff_array must have length N+M+1 */
	cdouble P_array [(M+1)*(M+1)];
	cdouble Q_array [(M+1)*(M+1)];

	for (size_t row_idx=0; row_idx<M; row_idx++){
		for (size_t col_idx=0; col_idx<M+1; col_idx++){
			P_array[row_idx*(M+1) + col_idx] = a_coeff_array[N-M+1 + row_idx + col_idx];
			Q_array[row_idx*(M+1) + col_idx] = a_coeff_array[N-M+1 + row_idx + col_idx];
		}
	}

	for (size_t col_idx=0; col_idx<M+1; col_idx++){
		Q_array[M*(M+1) + col_idx] = std::pow(z, M-col_idx);
		P_array[M*(M+1) + col_idx] = 0;
		for (size_t j=M-col_idx; j<N+1; j++){
			P_array[M*(M+1) + col_idx] += a_coeff_array[j - (M-col_idx)] * std::pow(z, j);
		}
	}

	cdouble P_det = determinant(P_array, M+1);
	cdouble Q_det = determinant(Q_array, M+1);

	return P_det/Q_det;
}

/* Solves the Faddeev equations
 * U = P*V + P*V*G*U
 * on the form L*U = R, where L and R are the left-
 * and right-handed sides of the equations, given by 
 * L = 1 - P*V*G
 * R = P*V
 * Since G is expressed in an SWP basis, we also must include the basis-transormation matrices C */
void faddeev_dense_solver(cdouble*  U_array,
					      cdouble*  G_array,
					      int*		q_com_idx_array,	size_t num_q_com,
					      int*      deuteron_idx_array, size_t num_deuteron_states,
					      size_t    Nalpha,
					      size_t 	Nq_WP,
					      size_t 	Np_WP,
					      double**  CT_RM_array,
					      double**  VC_CM_array,
					      double*   P123_sparse_val_array,
					      int*      P123_sparse_row_array,
					      size_t*   P123_sparse_col_array,
					      size_t    P123_sparse_dim,
					      const tnf_kernel_context& tnf_ctx){
	
	/* Stores A and K arrays */
	bool store_A_array = true;
	bool store_K_array = true;
	bool store_U_array = true;

	/* Dense dimension of 3N-channel */
	size_t dense_dim = Nalpha * Nq_WP * Np_WP;
	
	std::complex<double>* L_array = new cdouble [dense_dim*dense_dim];
	std::complex<double>* R_array = new cdouble [dense_dim*dense_dim];

	double* CPVC_col_array 		   = new double [dense_dim];
	int*    CPVC_row_to_nnz_array  = new int    [dense_dim];
	int*    CPVC_nnz_to_row_array  = new int    [dense_dim];
	
	for (size_t j=0; j<num_q_com; j++){

		/* Reset L- and R-arrays */
		for (size_t idx=0; idx<dense_dim*dense_dim; idx++){
			L_array[idx] = 0;
			R_array[idx] = 0;
		}

		/* Construct L- and R-arrays */
		//#pragma omp parallel for
		for (size_t idx_q_c=0; idx_q_c<Nq_WP; idx_q_c++){
			for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
				for (size_t idx_p_c=0; idx_p_c<Np_WP; idx_p_c++){
					size_t col_idx = idx_alpha_c*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_c;

					/* Reset CPVC-column array */
					for (size_t row_idx=0; row_idx<dense_dim; row_idx++){
						CPVC_col_array[row_idx] = 0;
						CPVC_row_to_nnz_array[row_idx] = -1;
						CPVC_nnz_to_row_array[row_idx] = -1;
					}

					/* Calculate CPVC-column */
					size_t CPVC_num_nnz = 0;
					calculate_CPVC_col(CPVC_col_array,
									   CPVC_row_to_nnz_array,
									   CPVC_nnz_to_row_array,
									   CPVC_num_nnz,
									   idx_alpha_c, idx_p_c, idx_q_c,
									   Nalpha, Nq_WP, Np_WP,
									   CT_RM_array,
									   VC_CM_array,
									   P123_sparse_val_array,
									   P123_sparse_row_array,
									   P123_sparse_col_array,
									   P123_sparse_dim,
									   tnf_ctx);

    	    		for (size_t row_idx=0; row_idx<dense_dim; row_idx++){
						L_array[row_idx*dense_dim + col_idx] = -CPVC_col_array[row_idx]*G_array[j*dense_dim + col_idx];
						R_array[row_idx*dense_dim + col_idx] =  CPVC_col_array[row_idx];
    	    		}

    	    		L_array[col_idx*dense_dim + col_idx] += 1;
				}
			}
    	}

		if (store_A_array){
			std::string A_arr_filename = "A_array_E_idx_" + std::to_string(j) + ".txt";
			store_array(R_array, dense_dim*dense_dim, A_arr_filename);
		}
		if (store_K_array){
			std::string K_arr_filename = "K_array_E_idx_" + std::to_string(j) + ".txt";
			store_array(L_array, dense_dim*dense_dim, K_arr_filename);
		}

		/* Solve */
		solve_MM(L_array, R_array, dense_dim);

		if (store_U_array){
			std::string U_arr_filename = "U_array_E_idx_" + std::to_string(j) + ".txt";
			store_array(R_array, dense_dim*dense_dim, U_arr_filename);
		}

		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){

				size_t idx_q_com		  = j;
				const elastic_on_shell_index ndos = make_elastic_on_shell_index(idx_d_row,
																				 idx_d_col,
																				 idx_q_com,
																				 deuteron_idx_array,
																				 q_com_idx_array,
																				 num_deuteron_states,
																				 num_q_com,
																				 Nq_WP,
																				 Np_WP);

				cdouble U_val = R_array[dense_packet_index(ndos.alpha_row,
														 ndos.q_idx,
														 elastic_bound_packet_index,
														 Nq_WP,
														 Np_WP)*dense_dim + ndos.col_storage_idx];

				U_array[ndos.value_storage_idx] = U_val;

				printf("   - U-matrix element for alpha'=%ld, alpha=%ld, q=%ld: %.10e + %.10ei \n", ndos.alpha_row, ndos.alpha_col, ndos.q_idx, U_array[ndos.value_storage_idx].real(), U_array[ndos.value_storage_idx].imag());
			}
		}
	}
	delete [] L_array;
	delete [] R_array;
	delete [] CPVC_col_array;
	delete [] CPVC_row_to_nnz_array;
	delete [] CPVC_nnz_to_row_array;
}

void pade_method_solve(cdouble*  U_array,
					   cdouble*  U_BU_array,
					   cdouble*  G_array,
					   int*		 q_com_idx_array,	 size_t num_q_com,
					   int*      deuteron_idx_array, size_t num_deuteron_states,
					   size_t    Nalpha,
					   size_t 	 Nq_WP,
					   size_t 	 Np_WP,
					   double**  CT_RM_array,
					   double**  VC_CM_array,
					   double*   P123_sparse_val_array,
					   int*      P123_sparse_row_array,
					   size_t*   P123_sparse_col_array,
					   size_t    P123_sparse_dim,
					   channel_os_indexing chn_os_indexing,
					   run_params run_parameters,
					   std::string file_identification,
					   const tnf_kernel_context& tnf_ctx){

	// [EN] This routine implements the matrix version of the WPCD multiple-scattering expansion used in the Miller
	// benchmarks: start from the driving term A=(C^T)PVC, generate Neumann terms A(GA)^n only on the physically
	// needed on-shell rows, and use Padé resummation to recover stable elastic amplitudes when the raw series
	// converges too slowly. / [CN] 该例程实现的是 Miller 基准工作中使用的 WPCD 多重散射矩阵算法：从驱动项 A=(C^T)PVC 出发，只在物理上需要的 on-shell 行上生成 A(GA)^n 的 Neumann 项，并用 Padé 重求和在原级数收敛过慢时恢复稳定的弹性振幅。
						   
	/* Print Pade-approximant convergences */
	bool print_PA_convergences = false;
	/* Print Neumann terms */
	bool print_neumann_terms   = false;
	/* Store Neumann terms */
	bool store_neumann_terms   = true;
	/* Store An-matrices */
	bool store_An_arrays 	   = false;
	/* Store CPVC-kernel to disk */
	bool store_CPVC_array	   = false;
	/* Store CPVC-kernel in memory (faster but much more memory intensive) */
	bool keep_CPVC_in_mem	   = false;
	
	/* Timekeeping variables */
	auto timestamp_start = std::chrono::system_clock::now();
	auto timestamp_end   = std::chrono::system_clock::now();

	/* Use as many threads as possible in MKL-GEMM */
	omp_set_num_threads(omp_get_max_threads());

	/* Number of on-shell nucleon-deuteron channels (deuteron states can mix, hence ^2) */
	size_t num_on_shell_A_rows = num_deuteron_states * num_q_com;
	size_t num_EL_A_vals 	   = num_deuteron_states * num_deuteron_states * num_q_com;	// Elastic elements
	size_t num_BU_A_vals 	   = 0;//num_deuteron_states * chn_os_indexing.num_BU_chns;	// Breakup elements
	
	/* Upper limit on polynomial approximation of Faddeev eq. */
	size_t NM_max = 14;
	size_t num_neumann_terms = 2*NM_max+1;

	// [EN] Chapter 7 rewrites the matrix AGS equation into a finite list of Neumann coefficients plus a Padé
	// resummation step. These arrays are the concrete storage for that coefficient pipeline. / [CN] 讲稿第 7 章把矩阵
	// AGS 方程改写成“有限个 Neumann 系数 + Padé 重求和”；下面这些数组就是这条系数流水线的具体存储。
	/* Coefficients for calculating Pade approximant */
	cdouble* a_coeff_array 	  = new cdouble [ num_neumann_terms * num_EL_A_vals];
	cdouble* a_BU_coeff_array = NULL;//new cdouble [ num_neumann_terms * num_BU_A_vals];

	/* Index holder array for storing Neumann terms */
	cdouble* 	 neumann_store_array 	= new cdouble [num_EL_A_vals];
	std::string* comment_array 	  		= new std::string [num_EL_A_vals];

	/* Dense dimension of 3N-channel */
	size_t dense_dim = Nalpha * Nq_WP * Np_WP;

	/* Allocate row-arrays for A*A^n, where A=(C^T)(P)(VC) 
	 * _array: current iteration of Neumann terms
	 * _array_prev: previous iteration of Neumann terms
	 * _array_comp: compact, past iteration of Neumann terms (compactified to contain only non-converged on-shell elements for faster gemm)
	 * _array_prod: compact, next iteration of Neumann terms (compactified to contain only non-converged on-shell elements for faster gemm) */
	double* re_A_An_row_array 	   = new double [dense_dim * num_on_shell_A_rows];
	double* im_A_An_row_array 	   = new double [dense_dim * num_on_shell_A_rows];
	double* re_A_An_row_array_prev = new double [dense_dim * num_on_shell_A_rows];
	double* im_A_An_row_array_prev = new double [dense_dim * num_on_shell_A_rows];
	double* re_A_An_row_array_comp = new double [dense_dim * num_on_shell_A_rows];
	double* im_A_An_row_array_comp = new double [dense_dim * num_on_shell_A_rows];
	double* re_A_An_row_array_prod = new double [dense_dim * num_on_shell_A_rows];
	double* im_A_An_row_array_prod = new double [dense_dim * num_on_shell_A_rows];
	
	/* Set A_An-arrays to zero */
	for (size_t i=0; i<dense_dim*num_on_shell_A_rows; i++){
		re_A_An_row_array[i] 	  = 0;
		im_A_An_row_array[i] 	  = 0;
		re_A_An_row_array_prev[i] = 0;
		im_A_An_row_array_prev[i] = 0;
		re_A_An_row_array_comp[i] = 0;
		im_A_An_row_array_comp[i] = 0;
		re_A_An_row_array_prod[i] = 0;
		im_A_An_row_array_prod[i] = 0;
	}

	/* Mapping vector from non-converged, compactified rows to full row-storage */
	size_t* A_An_indexing_array = new size_t [num_on_shell_A_rows];

	/* File-paths for storing A_An_row_array and on-shell neumann terms */
	std::string A_An_row_filename      = run_parameters.output_folder + "/An_rows" + file_identification + ".txt";
	std::string neumann_terms_filename = run_parameters.output_folder + "/neumann_terms" + file_identification + ".txt";

	/* Arrays to store Pade-approximants (PA) for each on-shell elastic elements */
	cdouble* pade_approximants_array      = new cdouble [num_EL_A_vals * (NM_max+1)];
	size_t*  pade_approximants_idx_array  = new size_t  [num_EL_A_vals];
	bool*    pade_approximants_conv_array = new bool    [num_EL_A_vals];
	size_t	 num_converged_elements		  = 0;

	/* Arrays to store Pade-approximants (PA) for each on-shell breakup elements */
	cdouble* pade_approximants_BU_array      = NULL;//new cdouble [num_BU_A_vals * (NM_max+1)];
	size_t*  pade_approximants_BU_idx_array  = NULL;//new size_t  [num_BU_A_vals];
	bool*    pade_approximants_BU_conv_array = NULL;//new bool    [num_BU_A_vals];
	size_t	 num_converged_BU_elements		 = 0;

	// [EN] The CPVC kernel is generated in chunks because the full dense object is far larger than the active
	// on-shell row set. Algebraically this still represents the same A=(C^T)PVC action described in the docs.
	// / [CN] CPVC 核按 chunk 生成，是因为完整稠密对象远大于当前活动的 on-shell 行集；但从代数上看，它仍然是文档里
	// 所说的同一个 A=(C^T)PVC 作用。
	/* Define CPVC-chunks size */
	size_t num_Gbytes_per_chunk  = 4;
	size_t num_bytes_per_chunk   = num_Gbytes_per_chunk * std::pow(1024,3);
    size_t num_bytes_in_CPVC_col = (sizeof(cdouble) * dense_dim);
    size_t num_cols_per_chunk    = num_bytes_per_chunk / num_bytes_in_CPVC_col;
    size_t num_chunks            = dense_dim / num_cols_per_chunk + 1;
    size_t block_size 			 = num_cols_per_chunk;

	/* From test-script to program notation */
	size_t    max_num_cols_in_mem = block_size;
	size_t	  num_col_chunks	  = num_chunks;
	double*   CPVC_cols_array     = new double [dense_dim * max_num_cols_in_mem];
	
	/* Allocate row- and column-arrays for (C^T)(P)(VC) */
	int 	 num_threads			    = omp_get_max_threads();
	double*  omp_CPVC_col_array  		= new double [dense_dim * num_threads];
	int*     omp_CPVC_row_to_nnz_array  = new int    [dense_dim * num_threads];
	int*     omp_CPVC_nnz_to_row_array  = new int    [dense_dim * num_threads];

	for (size_t idx_NDOS=0; idx_NDOS<num_EL_A_vals; idx_NDOS++){
		pade_approximants_conv_array[idx_NDOS] = false;
	}

	if (run_parameters.include_breakup_channels){
		num_BU_A_vals 	   				= num_deuteron_states * chn_os_indexing.num_BU_chns;
		a_BU_coeff_array 				= new cdouble [ num_neumann_terms * num_BU_A_vals];
		pade_approximants_BU_array      = new cdouble [num_BU_A_vals * (NM_max+1)];
		pade_approximants_BU_idx_array  = new size_t  [num_BU_A_vals];
		pade_approximants_BU_conv_array = new bool    [num_BU_A_vals];
		
		for (size_t idx_NDOS=0; idx_NDOS<num_BU_A_vals; idx_NDOS++){
			pade_approximants_BU_conv_array[idx_NDOS] = false;
		}
	}

	/* CPVC sparse arrays, used if keep_CPVC_in_mem=true */
	double* CPVC_v_array   = NULL;
	int*    CPVC_c_array   = NULL;
	int*    CPVC_r_array   = NULL;
	size_t* CPVC_csc_array = NULL;
	long long int* CPVC_c_array_LL   = NULL;
	long long int* CPVC_csc_array_LL = NULL;

	
	/* Precalculate kernel CPVC and keep in computer memory, if the option is selected
	 * !!! WARNING: EXTREMELY MEMORY INTENSIVE !!! */
	if (keep_CPVC_in_mem){
		printf("     - Precalculating CPVC-kernel in CSR-sparse format \n"); fflush(stdout);

		size_t CPVC_dim_est = P123_sparse_dim * 100;
		double CPVC_GB_est  = (double) CPVC_dim_est * (2*sizeof(int) + sizeof(double)) / std::pow(1024,3);
		printf("       - Estimated kernel size: %zu non-zero elements (%.2f GB in COO format)\n", CPVC_dim_est, CPVC_GB_est); fflush(stdout);

		/* CPVC-sparse arrays in COO format PER THREAD */
		size_t*  omp_CPVC_dim 	  = new size_t  [num_threads];
		size_t*  omp_CPVC_nnz 	  = new size_t  [num_threads];
		double** omp_CPVC_v_array = new double* [num_threads];
		int**    omp_CPVC_c_array = new int*    [num_threads];
		int**    omp_CPVC_r_array = new int*    [num_threads];
		for (int i=0; i<num_threads; i++){
			size_t omp_dim_est  = P123_sparse_dim * 10 / num_threads; // This is an estimate based on what we see typically, distribute among threads
			omp_CPVC_dim[i]     = omp_dim_est;
			omp_CPVC_nnz[i]     = 0;
			omp_CPVC_v_array[i] = new double [omp_dim_est];
			omp_CPVC_c_array[i] = new int    [omp_dim_est];
			omp_CPVC_r_array[i] = new int    [omp_dim_est];
		}

		/* Fill omp-arrays with nnz elements of CPVC-kernel */
		printf("       - Precalculating ... \n"); fflush(stdout);
		timestamp_start = std::chrono::system_clock::now();
		#pragma omp parallel
		{
			size_t  thread_idx        = omp_get_thread_num();
			double* col_array  	   	  = new double [dense_dim];
			int*    row_to_nnz_array  = new int    [dense_dim];
			int*    nnz_to_row_array  = new int    [dense_dim];

			/* Manually set up parallel for loop, this lets me avoid having to sort columns after matrix construction */
			size_t idx_col_start = (dense_dim/num_threads) *  thread_idx;
			size_t idx_col_end   = (dense_dim/num_threads) * (thread_idx+1) + (thread_idx==num_threads-1)*(dense_dim%num_threads);

			for (size_t idx_col=idx_col_start; idx_col<idx_col_end; idx_col++){
				for (size_t i=0; i<dense_dim; i++){
					col_array[i] = 0;
					row_to_nnz_array[i] = -1;
					nnz_to_row_array[i] = -1;
				}

				size_t idx_alpha_c =  idx_col / (Np_WP*Nq_WP);
				size_t idx_q_c     = (idx_col % (Np_WP*Nq_WP)) /  Np_WP;
				size_t idx_p_c     =  idx_col %  Np_WP;

				/* Calculate CPVC-column */
				size_t num_nnz = 0;
				calculate_CPVC_col(col_array,
								   row_to_nnz_array,
								   nnz_to_row_array,
								   num_nnz,
								   idx_alpha_c, idx_p_c, idx_q_c,
								   Nalpha, Nq_WP, Np_WP,
								   CT_RM_array,
								   VC_CM_array,
								   P123_sparse_val_array,
								   P123_sparse_row_array,
								   P123_sparse_col_array,
								   P123_sparse_dim,
								   tnf_ctx);

				/* Lengthen array by 100*dense_dim if the array cannot fit another dense_dim nnz-elements
				 * (i.e. max nnz elements from next loop-iteration) */
				if ( omp_CPVC_nnz[thread_idx]+num_nnz+dense_dim >= omp_CPVC_dim[thread_idx] ){
					size_t steplength = 100*dense_dim;
					increase_sparse_array_size(&omp_CPVC_v_array[thread_idx], omp_CPVC_dim[thread_idx], steplength);
					increase_sparse_array_size(&omp_CPVC_r_array[thread_idx], omp_CPVC_dim[thread_idx], steplength);
					increase_sparse_array_size(&omp_CPVC_c_array[thread_idx], omp_CPVC_dim[thread_idx], steplength);
					omp_CPVC_dim[thread_idx] += steplength;
				}

				/* Append nnz elements to sparse omp-array */
				size_t curr_nnz = omp_CPVC_nnz[thread_idx];
				for (int nnz=0; nnz<num_nnz; nnz++){
					omp_CPVC_v_array[thread_idx][curr_nnz+nnz] = col_array[nnz_to_row_array[nnz]];
					omp_CPVC_r_array[thread_idx][curr_nnz+nnz] = nnz_to_row_array[nnz];
					omp_CPVC_c_array[thread_idx][curr_nnz+nnz] = idx_col;
				}
				omp_CPVC_nnz[thread_idx] += num_nnz;
			}
		}
		size_t CPVC_num_nnz = 0;
		for (int i=0; i<num_threads; i++){
			CPVC_num_nnz += omp_CPVC_nnz[i];
		}
		double CPVC_GB_true  = (double) CPVC_num_nnz * (2*sizeof(int) + sizeof(double)) / std::pow(1024,3);
		double CPVC_density  = (double) 100.0 * CPVC_num_nnz / std::pow(dense_dim, 2);
		printf("         - Actual kernel size: %zu non-zero elements (%.2f GB in COO format) (%.2f%% density)\n", CPVC_num_nnz, CPVC_GB_true, CPVC_density); fflush(stdout);
		timestamp_end = std::chrono::system_clock::now();
		std::chrono::duration<double>  time_CPVC_construction = timestamp_end - timestamp_start;
		printf("         - Time spent:     %.6f \n", time_CPVC_construction.count()); fflush(stdout);
		printf("         - Done. \n"); fflush(stdout);

		/* Consolidate omp-arrays into a single COO-format array
		 * NOTE: This can be divided into three sections that deallocate e.g. omp-value array
		 * before allocating col-array such that memory is used more efficiently. */
		printf("       - Consolidating distributed arrays ... \n"); fflush(stdout);
		timestamp_start = std::chrono::system_clock::now();

		/* Consolidate values and deallocate parallel memory */
		CPVC_v_array = new double [CPVC_num_nnz];
		size_t idx_nnz = 0;
		for (int i=0; i<num_threads; i++){
			for (size_t j=0; j<omp_CPVC_nnz[i]; j++){
				CPVC_v_array[idx_nnz] = omp_CPVC_v_array[i][j];
				idx_nnz += 1;
			}
		}
		for (int i=0; i<num_threads; i++){
			delete [] omp_CPVC_v_array[i];
		}

		/* Consolidate column-indices and deallocate parallel memory  */
		CPVC_c_array = new int [CPVC_num_nnz];
		idx_nnz = 0;
		for (int i=0; i<num_threads; i++){
			for (size_t j=0; j<omp_CPVC_nnz[i]; j++){
				CPVC_c_array[idx_nnz] = omp_CPVC_c_array[i][j];
				idx_nnz += 1;
			}
		}
		for (int i=0; i<num_threads; i++){
			delete [] omp_CPVC_c_array[i];
		}

		/* Consolidate row-indices and deallocate parallel memory  */
		CPVC_r_array = new int [CPVC_num_nnz];
		idx_nnz = 0;
		for (int i=0; i<num_threads; i++){
			for (size_t j=0; j<omp_CPVC_nnz[i]; j++){
				CPVC_r_array[idx_nnz] = omp_CPVC_r_array[i][j];
				idx_nnz += 1;
			}
		}
		for (int i=0; i<num_threads; i++){
			delete [] omp_CPVC_r_array[i];
		}

		/* Deallocate parallel pointer-arrays */
		delete [] omp_CPVC_dim;
		delete [] omp_CPVC_nnz;
		delete [] omp_CPVC_v_array;
		delete [] omp_CPVC_c_array;
		delete [] omp_CPVC_r_array;

		timestamp_end = std::chrono::system_clock::now();
		std::chrono::duration<double>  time_CPVC_consolidation = timestamp_end - timestamp_start;
		printf("         - Time spent:     %.6f \n", time_CPVC_consolidation.count()); fflush(stdout);
		printf("         - Done. \n"); fflush(stdout);

		std::string filename = "kernel.h5";
		/* Store sparse CPVC-array to file */
		printf("       - Writing kernel to h5 ... \n"); fflush(stdout);
		timestamp_start = std::chrono::system_clock::now();
		store_sparse_matrix_h5(CPVC_v_array,
							   CPVC_r_array,
							   CPVC_c_array,
							   CPVC_num_nnz,
							   dense_dim,
					   		   filename,
							   false);
		timestamp_end = std::chrono::system_clock::now();
		std::chrono::duration<double>  time_CPVC_storage = timestamp_end - timestamp_start;
		printf("         - Successs. \n"); fflush(stdout);
		printf("         - Time spent:     %.6f \n", time_CPVC_storage.count()); fflush(stdout);
		printf("         - Done. \n"); fflush(stdout);

		/* Read sparse CPVC-array from file */
		printf("       - Reading kernel from h5 ... \n"); fflush(stdout);
		timestamp_start = std::chrono::system_clock::now();
		delete [] CPVC_v_array;
		CPVC_v_array = NULL;
		delete [] CPVC_r_array;
		CPVC_r_array = NULL;
		delete [] CPVC_c_array;
		CPVC_c_array = NULL;
		read_sparse_matrix_h5(&CPVC_v_array,
							   &CPVC_r_array,
							   &CPVC_c_array,
							   CPVC_num_nnz,
							   dense_dim,
					   		   filename,
							   false);
		timestamp_end = std::chrono::system_clock::now();
		std::chrono::duration<double>  time_CPVC_reading = timestamp_end - timestamp_start;
		printf("         - Successs. \n"); fflush(stdout);
		printf("         - Time spent:     %.6f \n", time_CPVC_reading.count()); fflush(stdout);
		printf("         - Done. \n"); fflush(stdout);

		printf("       - Done. \n"); fflush(stdout);
	}

	/* Set initial values for A_Kn_row_array, where K^n=1 for n=0 */
	// [EN] The first stored object is a_0 = A, i.e. the driving term before any rescattering. In the notes this is
	// the first term of the Neumann series, representing a single application of the kernel without additional G
	// propagation. / [CN] 这里首先存储的是 a_0 = A，也就是没有任何再散射前的驱动项；按讲稿的说法，它对应 Neumann 级数的
	// 第一项，表示只施加一次核而没有额外的 G 传播。
	printf("     - Working on Pade approximant P[N,M] for N=%d, M=%d \n",0,0); fflush(stdout);
	printf("       - Calculating on-shell rows of A*K^n for n=%d. \n", 0); fflush(stdout);
	timestamp_start = std::chrono::system_clock::now();
	/* Calculate CPVC-row and write to A_Kn_row_array_prev */
	calculate_all_CPVC_rows(re_A_An_row_array_prev,
							q_com_idx_array, num_q_com,
			   				deuteron_idx_array, num_deuteron_states,
							Nalpha, Nq_WP, Np_WP,
							CT_RM_array,
							VC_CM_array,
							P123_sparse_val_array,
							P123_sparse_row_array,
							P123_sparse_col_array,
							P123_sparse_dim,
							tnf_ctx);
	timestamp_end = std::chrono::system_clock::now();
	std::chrono::duration<double> time = timestamp_end - timestamp_start;
	printf("         - Time generating CPVC-rows:     %.6f \n", time.count()); fflush(stdout);
	printf("         - Done \n"); fflush(stdout);

	/* Define Neumann-printout commentary array */
	for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
		for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				const elastic_on_shell_index ndos = make_elastic_on_shell_index(idx_d_row,
																				 idx_d_col,
																				 idx_q_com,
																				 deuteron_idx_array,
																				 q_com_idx_array,
																				 num_deuteron_states,
																				 num_q_com,
																				 Nq_WP,
																				 Np_WP);

				/* Store indices of coefficient */
				comment_array[ndos.value_storage_idx] =   "\t#\t alpha'-idx=" + std::to_string(ndos.alpha_row)
													   + "\t alpha-idx="+ std::to_string(ndos.alpha_col)
													   + "\t q-idx=" + std::to_string(ndos.q_idx);
			}
		}
	}
	
	// [EN] a_0 is the driving term itself: one application of A=(C^T)PVC with no intermediate propagation. Every
	// later coefficient adds one more G propagation and one more rescattering by the same kernel. / [CN] a_0 就是驱动项本身：
	// 只作用一次 A=(C^T)PVC，中间没有传播；之后每一阶都会再多插入一次 G 传播和一次同一核的再散射。
	/* First Neumann-term */
	printf("       - Extracting on-shell Neumann-series terms a_n=A*K^n for n=%d. \n",0); fflush(stdout);
	/* Extract elastic terms */
	for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
		for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				const elastic_on_shell_index ndos = make_elastic_on_shell_index(idx_d_row,
																				 idx_d_col,
																				 idx_q_com,
																				 deuteron_idx_array,
																				 q_com_idx_array,
																				 num_deuteron_states,
																				 num_q_com,
																				 Nq_WP,
																				 Np_WP);

				/* Calculate coefficient */
				cdouble a_coeff = re_A_An_row_array_prev[ndos.row_storage_idx*dense_dim + ndos.col_storage_idx];

				/* Store coefficient */
				a_coeff_array[ndos.value_storage_idx*num_neumann_terms] = a_coeff;

				/* Store coefficient in print-to-file format */
				neumann_store_array[ndos.value_storage_idx] = a_coeff;
				
				if (print_neumann_terms){
					printf("         - Neumann term %d for alpha'=%ld, alpha=%ld, q=%ld: %.16e + %.16ei \n", 0, ndos.alpha_row, ndos.alpha_col, ndos.q_idx, a_coeff.real(), a_coeff.imag());
					fflush(stdout);
				}
			}
		}
	}
	/* Extract breakup terms */
	if (run_parameters.include_breakup_channels){
		for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
			int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[idx_q_com];
			int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[idx_q_com+1];
			for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
				for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){
					
					size_t idx_alpha_NDOS = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 0];
					size_t idx_q_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 1];
					size_t idx_p_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 2];

					size_t idx_row_NDOS   = idx_d_row*num_q_com + idx_q_com;
					size_t idx_col_NDOS   = idx_alpha_NDOS*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;

					/* Calculate coefficient */
					cdouble a_BU_coeff = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx_col_NDOS];
					
					/* Store coefficient */
					size_t idx_NDOS = breakup_value_storage_index(idx_d_row,
																 idx_BU_chn,
																 chn_os_indexing.num_BU_chns);
					a_BU_coeff_array[idx_NDOS*num_neumann_terms] = a_BU_coeff;
				}
			}
		}
	}
	printf("         - Done \n"); fflush(stdout);

	if (store_An_arrays){
		printf("       - Storing matrix A*K^n for n=%d to output-folder. \n", 0); fflush(stdout);
		std::string array_seperator_text = "n = " + std::to_string(0);
		store_sep_complex_matrix(re_A_An_row_array_prev,
								 im_A_An_row_array_prev,
    	                         num_on_shell_A_rows,
						         dense_dim,
						         dense_dim,
						         A_An_row_filename,
						         true,
							     array_seperator_text);
		printf("         - Done \n");
	}
	if (store_neumann_terms){
		printf("       - Storing on-shell Neumann-series terms a_n=A*K^n for n=%d to output-folder. \n", 0); fflush(stdout);
		std::string array_seperator_text = "n = " + std::to_string(0);
		store_complex_vector_with_comments(neumann_store_array,
										   comment_array,
										   num_deuteron_states*num_deuteron_states*num_q_com,
										   neumann_terms_filename,
										   true,
										   array_seperator_text);
		printf("         - Done \n");
	}

	/* Loop over number of Pade-terms we use */
	for (size_t NM=0; NM<NM_max+1; NM++){
		if (num_converged_elements==num_EL_A_vals){
			printf("     - Convergence reached for all on-shell elements! \n"); fflush(stdout);
			break;
		}

		if (NM!=0){
			printf("     - Working on Pade approximant P[N,M] for N=%ld, M=%ld \n",NM,NM); fflush(stdout);
		}
		
		size_t counter_array [100];
		for (size_t i=0; i<100; i++){
			counter_array[i] = 0;
		}

		/* Time-keeper array for parallel environment */
		double*  times_array = new double [3*num_threads];
		
		for (int n=2*NM-1; n<2*NM+1; n++){
			printf("       - Working on Neumann-terms for n=%d. \n", n); fflush(stdout);
			/* We've already done n=0 above */
			if (n<=0){
				continue;
			}
			
			/* Initialise time-profile variables */
			double time_resolvent        = 0;
			double time_CPVC_cols        = 0;
			double time_An_CPVC_multiply = 0;
			double time_neumann          = 0;

			double timestamp_neumann_start = omp_get_wtime();

			// [EN] The Neumann recursion is implemented exactly as in the lecture notes: first multiply the previous
			// term by the diagonal channel resolvent G, then apply the CPVC kernel to generate the next rescattering
			// contribution. / [CN] Neumann 递推严格按照讲稿里的顺序实现：先把上一阶乘上对角的通道分辨算符 G，再施加 CPVC
			// 核，生成下一阶再散射贡献。
			double timestamp_resolvent_start = omp_get_wtime();
			printf("       - Multiplying in resolvent with An. \n"); fflush(stdout);
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
					size_t idx_row_NDOS = idx_d_row*num_q_com + idx_q_com;

					/* Multiply An by G */
					for (size_t idx=0; idx<dense_dim; idx++){
						double re_An = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx];
						double im_An = im_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx];
						double re_G  = G_array[idx_q_com*dense_dim + idx].real();
						double im_G  = G_array[idx_q_com*dense_dim + idx].imag();
						re_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx] = re_An*re_G - im_An*im_G;
						im_A_An_row_array_prev[idx_row_NDOS*dense_dim + idx] = re_An*im_G + im_An*re_G;
					}
				}
			}
			double timestamp_resolvent_end    = omp_get_wtime();
			time_resolvent = timestamp_resolvent_end - timestamp_resolvent_start;

			// [EN] Once some on-shell amplitudes have converged, there is no value in pushing their full dense rows
			// through later rescattering steps. Compacting only the non-converged rows is a pure algebraic shortcut:
			// it preserves the exact iteration on the active rows while reducing GEMM cost. / [CN] 当部分 on-shell
			// 振幅已经收敛后，就没有必要再把它们对应的整条稠密行送入后续再散射步骤；这里只压缩未收敛的行是纯代数层面的加速，在保持活动行迭代完全一致的同时降低了 GEMM 成本。
			size_t num_non_conv_rows = 0;
			for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
				for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
					size_t idx_row_NDOS = elastic_row_storage_index(idx_d_row, idx_q_com, num_q_com);
					if (row_has_only_converged_targets(idx_d_row,
													   idx_q_com,
													   num_deuteron_states,
													   num_q_com,
													   pade_approximants_conv_array,
													   pade_approximants_BU_conv_array,
													   chn_os_indexing,
													   run_parameters.include_breakup_channels)==false){
						for (size_t i=0; i<dense_dim; i++){
							re_A_An_row_array_comp[num_non_conv_rows*dense_dim + i] = re_A_An_row_array_prev[idx_row_NDOS*dense_dim + i];
							im_A_An_row_array_comp[num_non_conv_rows*dense_dim + i] = im_A_An_row_array_prev[idx_row_NDOS*dense_dim + i];
						}
						A_An_indexing_array[num_non_conv_rows] = idx_row_NDOS;
						num_non_conv_rows += 1;
					}
				}
			}
			
			printf("       - Calculating on-shell rows of A*K^n for n=%d. \n", n); fflush(stdout);
			if (keep_CPVC_in_mem==false){
				// [EN] The kernel columns are regenerated in chunks so the dense GEMM path can stream through the
				// active part of CPVC without materializing the full dense matrix. / [CN] 这里按块重建核列，这样稠密 GEMM
				// 路径就能流式处理 CPVC 的活动部分，而不必把整个稠密矩阵完整落在内存中。
				for (size_t idx_col_chunk=0; idx_col_chunk<num_col_chunks; idx_col_chunk++){

					double timestamp_CPVC_chunk_start = omp_get_wtime();

					size_t idx_col_start =  idx_col_chunk    * max_num_cols_in_mem;
					size_t idx_col_end   = (idx_col_chunk+1) * max_num_cols_in_mem;

					if (idx_col_end>dense_dim){
						idx_col_end = dense_dim;
					}

					size_t cols_in_chunk = idx_col_end - idx_col_start;

					/* Reset CPVC-columns when array is filled */
					for (size_t idx=0; idx<dense_dim * max_num_cols_in_mem; idx++){
						CPVC_cols_array[idx] = 0;
					}

					#pragma omp parallel //num_threads(1)
					{

					size_t  thread_idx             = omp_get_thread_num();
					double* CPVC_col_array  	   = &omp_CPVC_col_array	    [thread_idx*dense_dim];
					int*    CPVC_row_to_nnz_array  = &omp_CPVC_row_to_nnz_array [thread_idx*dense_dim];
					int*    CPVC_nnz_to_row_array  = &omp_CPVC_nnz_to_row_array [thread_idx*dense_dim];

					#pragma omp for
					for (size_t idx_col=idx_col_start; idx_col<idx_col_end; idx_col++){
						size_t idx_alpha_c = idx_col / (Np_WP*Nq_WP);
						size_t idx_q_c     = (idx_col % (Np_WP*Nq_WP)) /  Np_WP;
						size_t idx_p_c     = idx_col %  Np_WP;

						/* Calculate CPVC-column */
						size_t CPVC_num_nnz = 0;
						size_t CPVC_col_idx = idx_col - idx_col_start;
						calculate_CPVC_col(&CPVC_cols_array[CPVC_col_idx*dense_dim],//CPVC_col_array,
										   CPVC_row_to_nnz_array,
										   CPVC_nnz_to_row_array,
										   CPVC_num_nnz,
										   idx_alpha_c, idx_p_c, idx_q_c,
										   Nalpha, Nq_WP, Np_WP,
										   CT_RM_array,
										   VC_CM_array,
										   P123_sparse_val_array,
										   P123_sparse_row_array,
										   P123_sparse_col_array,
										   P123_sparse_dim,
										   tnf_ctx);
						counter_array[thread_idx] += CPVC_num_nnz;
					}
					}
					double timestamp_CPVC_chunk_end   = omp_get_wtime();
					time_CPVC_cols += timestamp_CPVC_chunk_end - timestamp_CPVC_chunk_start;

					double beta  = 0;
					double alpha = 1;
					int M    = num_non_conv_rows;// num_on_shell_A_rows
					int N    = cols_in_chunk;// max_num_cols_in_mem;
					int K    = dense_dim;
					int lda  = dense_dim;
					int ldb  = dense_dim;//max_num_cols_in_mem;
					int ldc  = dense_dim;
					double* re_A = &re_A_An_row_array_comp[0];
					double* im_A = &im_A_An_row_array_comp[0];
					double* B 	 = &CPVC_cols_array[0];
					double* re_C = &re_A_An_row_array_prod[idx_col_start];
					double* im_C = &im_A_An_row_array_prod[idx_col_start];

					double timestamp_gemm_start = omp_get_wtime();
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, re_A, lda, B, ldb, beta, re_C, ldc);	// real multiplication
					cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, M, N, K, alpha, im_A, lda, B, ldb, beta, im_C, ldc);	// imag multiplication
					double timestamp_gemm_end   = omp_get_wtime();
					time_An_CPVC_multiply += timestamp_gemm_end - timestamp_gemm_start;
				}
			}
			else{
				double timestamp_gemm_start = omp_get_wtime();
				const char   ordering = 'R';
				const char   trans 	  = 'T';
				const double alpha 	  = 1.0;
				/* Transpose An before sparse multiplication */
				inplace_transpose(re_A_An_row_array_comp, num_non_conv_rows, dense_dim);
				inplace_transpose(im_A_An_row_array_comp, num_non_conv_rows, dense_dim);
				/* Multiply CPVC with An using sparse multiplication */
				//dot_MM_sparse(CPVC_v_array, CPVC_c_array_LL, CPVC_csc_array_LL, re_A_An_row_array_comp, re_A_An_row_array_prod, dense_dim, dense_dim, num_non_conv_rows, true);
				//dot_MM_sparse(CPVC_v_array, CPVC_c_array_LL, CPVC_csc_array_LL, re_A_An_row_array_comp, im_A_An_row_array_prod, dense_dim, dense_dim, num_non_conv_rows, true);
				/* Transpose An+1 after sparse multiplication */
				inplace_transpose(re_A_An_row_array_prod, dense_dim, num_non_conv_rows);
				inplace_transpose(im_A_An_row_array_prod, dense_dim, num_non_conv_rows);
				double timestamp_gemm_end   = omp_get_wtime();
				time_An_CPVC_multiply += timestamp_gemm_end - timestamp_gemm_start;
			}

			/* Write compact format back to full format */
			for (size_t r=0; r<num_non_conv_rows; r++){
				size_t idx_row_NDOS = A_An_indexing_array[r];
				for (size_t i=0; i<dense_dim; i++){
					re_A_An_row_array[idx_row_NDOS*dense_dim + i] = re_A_An_row_array_prod[r*dense_dim + i];
					im_A_An_row_array[idx_row_NDOS*dense_dim + i] = im_A_An_row_array_prod[r*dense_dim + i];
				}
			}
			double timestamp_neumann_end = omp_get_wtime();
			time_neumann = timestamp_neumann_end - timestamp_neumann_start;

			printf("         - Time multiplying An with G:    %.6f s \n", time_resolvent);
			printf("         - Time generating CPVC-cols:     %.6f s \n", time_CPVC_cols);
			printf("         - Time multiplying An with CPVC: %.6f s \n", time_An_CPVC_multiply);
			printf("         - Total time:                    %.6f s \n", time_neumann);
			printf("         - Done \n"); fflush(stdout);

			/* Rewrite previous A_An with current A_An */
			for (size_t i=0; i<num_on_shell_A_rows*dense_dim; i++){
				re_A_An_row_array_prev[i] = re_A_An_row_array[i];
				im_A_An_row_array_prev[i] = im_A_An_row_array[i];
			}

			// [EN] Only the on-shell entries are fed into the Padé build. This follows the Miller workflow: compute the
			// physically required scalar coefficient sequence a_n for each elastic/breakup amplitude, then resum those
			// short sequences rather than the full matrix. / [CN] 只有 on-shell 元素会被送入 Padé 构造。这正对应 Miller
			// 工作流：先为每个弹性/破裂振幅提取物理上需要的标量系数序列 a_n，再对这些短序列做重求和，而不是对整块矩阵重求和。
			printf("       - Extracting on-shell Neumann-series terms a_n=A*K^n for n=%d. \n", n); fflush(stdout);
			/* Extract elastic terms */
			for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
				for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
					for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
						const elastic_on_shell_index ndos = make_elastic_on_shell_index(idx_d_row,
																						 idx_d_col,
																						 idx_q_com,
																						 deuteron_idx_array,
																						 q_com_idx_array,
																						 num_deuteron_states,
																						 num_q_com,
																						 Nq_WP,
																						 Np_WP);

						/* Check if we've already reached convergence for this on-shell element */
						if (pade_approximants_conv_array[ndos.value_storage_idx]==true){
							continue;
						}

						/* Calculate coefficient */
						cdouble a_coeff = {re_A_An_row_array[ndos.row_storage_idx*dense_dim + ndos.col_storage_idx], im_A_An_row_array[ndos.row_storage_idx*dense_dim + ndos.col_storage_idx]};

						/* Store coefficient */
						a_coeff_array[ndos.value_storage_idx*num_neumann_terms + n] = a_coeff;

						/* Store coefficient in print-to-file format */
						neumann_store_array[ndos.value_storage_idx] = a_coeff;

						if (print_neumann_terms){
							printf("         - Neumann term %d for alpha'=%ld, alpha=%ld, q=%ld: %.16e + %.16ei \n", n, ndos.alpha_row, ndos.alpha_col, ndos.q_idx, a_coeff.real(), a_coeff.imag());
						}
					}
				}
			}
			/* Extract breakup terms */
			if (run_parameters.include_breakup_channels){
				for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
					int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[idx_q_com];
					int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[idx_q_com+1];
					for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
						for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){
							
							size_t idx_alpha_NDOS = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 0];
							size_t idx_q_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 1];
							size_t idx_p_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 2];

							size_t idx_row_NDOS   = idx_d_row*num_q_com + idx_q_com;
							size_t idx_col_NDOS   = idx_alpha_NDOS*Nq_WP*Np_WP + idx_q_NDOS*Np_WP + idx_p_NDOS;

							size_t idx_NDOS = breakup_value_storage_index(idx_d_row,
																		 idx_BU_chn,
																		 chn_os_indexing.num_BU_chns);

							/* Check if we've already reached convergence for this on-shell element */
							if (pade_approximants_BU_conv_array[idx_NDOS]==true){
								continue;
							}

							/* Calculate coefficient */
							cdouble a_BU_coeff = {re_A_An_row_array[idx_row_NDOS*dense_dim + idx_col_NDOS], im_A_An_row_array[idx_row_NDOS*dense_dim + idx_col_NDOS]};
							
							/* Store coefficient */
							a_BU_coeff_array[idx_NDOS*num_neumann_terms + n] = a_BU_coeff;
						}
					}
				}
			}
			printf("         - Done \n");

			if (store_An_arrays){
			printf("       - Storing matrix A*K^n for n=%d to output-folder. \n", n); fflush(stdout);
			std::string array_seperator_text = "n = " + std::to_string(n);
		    store_sep_complex_matrix(re_A_An_row_array,
									 im_A_An_row_array,
                                     num_on_shell_A_rows,
		    						 dense_dim,
		    						 dense_dim,
		    						 A_An_row_filename,
		    						 false,
									 array_seperator_text);
			printf("         - Done \n");
			}
			if (store_neumann_terms){
				printf("       - Storing on-shell Neumann-series terms a_n=A*K^n for n=%d to output-folder. \n", n); fflush(stdout);
				std::string array_seperator_text = "n = " + std::to_string(n);
				store_complex_vector_with_comments(neumann_store_array,
										   		   comment_array,
												   num_deuteron_states*num_deuteron_states*num_q_com,
												   neumann_terms_filename,
												   false,
												   array_seperator_text);
				printf("         - Done \n");
			}
		}
		delete [] times_array;

		printf("       - Calculating Pade approximants PA[%ld,%ld]. \n", NM, NM); fflush(stdout);
		// [EN] Padé resummation is what turns a slowly convergent or even divergent Neumann history into stable
		// amplitudes. Each on-shell element is treated independently because different channels can converge at very
		// different rates. / [CN] Padé 重求和是把收敛缓慢甚至发散的 Neumann 历史转化为稳定振幅的关键。这里每个 on-shell
		// 元素独立处理，因为不同通道的收敛速度可能相差很大。
		/* Calculate Pade approximants (PA) for elastic amplitudes */
		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
				for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){

					size_t idx_NDOS = elastic_value_storage_index(idx_d_row,
															   idx_d_col,
															   idx_q_com,
															   num_deuteron_states,
															   num_q_com);

					/* Check and skip if we've already reached convergence for this on-shell element */
					if (pade_approximants_conv_array[idx_NDOS]==true){
						continue;
					}

					/* Calculate and append PA */
					cdouble PA = pade_approximant(&a_coeff_array[idx_NDOS*num_neumann_terms], NM, NM, 1);

					pade_approximants_array[idx_NDOS*(NM_max+1) + NM] = PA;
					
					/* See if we've reached convergence with this iteration */
					size_t idx_best_PA = 0;
					bool convergence_reached = false;

					/* Find minimum PA from previous calculations */
					double min_PA_diff = 1;
					for (int NM_prev=0; NM_prev<NM; NM_prev++){
						cdouble PA_prev = pade_approximants_array[idx_NDOS*(NM_max+1) + NM_prev];

						/* Calculate difference between PAs from previous PA-calculations */
						double PA_diff_prev = std::abs(PA_prev - pade_approximants_array[idx_NDOS*(NM_max+1) + NM_prev-1]);

						/* Ignore PA_diff_prev if numerically equal to the previous PA_diff, overwrite if smaller than min_PA_diff */
						if (PA_diff_prev<min_PA_diff && PA_diff_prev>1e-15){
							idx_best_PA = NM_prev;
							min_PA_diff = PA_diff_prev;
						}
					}
					/* See if current PA is better/worse than previous minimum */
					double PA_diff_curr = std::abs(PA - pade_approximants_array[idx_NDOS*(NM_max+1) + NM - 1]);
					
					/* Ignore PA_diff_curr if numerically equal to the previous PA_diff_curr, overwrite if smaller than min_PA_diff */
					if (PA_diff_curr<min_PA_diff && PA_diff_curr>1e-15){
						idx_best_PA = NM;
						min_PA_diff = PA_diff_curr;
					}

					/* Criterias for convergence¨
					 * If any are fulfilled, we set convergence to true for current on-shell element */
					bool convergence_criteria_0 = (NM==NM_max);													// Cannot go past max NM
					bool convergence_criteria_1 = (NM-idx_best_PA>4);												// If nothing better is found in the last 3 PAs, we assume we found the best
					bool convergence_criteria_2 = (min_PA_diff<1e-6*std::abs(pade_approximants_array[idx_NDOS*(NM_max+1) + idx_best_PA]));	// If the difference is less than the 4th significant digit, we assume "good enough"
					bool convergence_criteria_3 = (min_PA_diff<1e-7);												// If we are below single precision resolution, assume convergence
					
					if (convergence_criteria_0 ||
						convergence_criteria_1 ||
						convergence_criteria_2 ||
						convergence_criteria_3){
						pade_approximants_conv_array[idx_NDOS] = true;
						pade_approximants_idx_array[idx_NDOS]  = idx_best_PA;
						num_converged_elements += 1;
					}
				}
			}
		}
		/* Calculate Pade approximants (PA) for breakup amplitudes */
		if (run_parameters.include_breakup_channels){
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[idx_q_com];
				int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[idx_q_com+1];
				for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
					for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){

						size_t idx_NDOS = breakup_value_storage_index(idx_d_row,
																 idx_BU_chn,
																 chn_os_indexing.num_BU_chns);

						/* Check and skip if we've already reached convergence for this on-shell element */
						if (pade_approximants_BU_conv_array[idx_NDOS]==true){
							continue;
						}

						/* Calculate and append PA */
						cdouble PA = pade_approximant(&a_BU_coeff_array[idx_NDOS*num_neumann_terms], NM, NM, 1);

						pade_approximants_BU_array[idx_NDOS*(NM_max+1) + NM] = PA;
						
						/* See if we've reached convergence with this iteration */
						size_t idx_best_PA = 0;
						bool convergence_reached = false;

						/* Find minimum PA from previous calculations */
						double min_PA_diff = 1;
						for (int NM_prev=0; NM_prev<NM; NM_prev++){
							cdouble PA_prev = pade_approximants_BU_array[idx_NDOS*(NM_max+1) + NM_prev];

							/* Calculate difference between PAs from previous PA-calculations */
							double PA_diff_prev = std::abs(PA_prev - pade_approximants_BU_array[idx_NDOS*(NM_max+1) + NM_prev-1]);

							/* Ignore PA_diff_prev if numerically equal to the previous PA_diff, overwrite if smaller than min_PA_diff */
							if (PA_diff_prev<min_PA_diff && PA_diff_prev>1e-15){
								idx_best_PA = NM_prev;
								min_PA_diff = PA_diff_prev;
							}
						}
						/* See if current PA is better/worse than previous minimum */
						double PA_diff_curr = std::abs(PA - pade_approximants_BU_array[idx_NDOS*(NM_max+1) + NM - 1]);
						
						/* Ignore PA_diff_curr if numerically equal to the previous PA_diff_curr, overwrite if smaller than min_PA_diff */
						if (PA_diff_curr<min_PA_diff && PA_diff_curr>1e-15){
							idx_best_PA = NM;
							min_PA_diff = PA_diff_curr;
						}

						/* Criterias for convergence¨
						* If any are fulfilled, we set convergence to true for current on-shell element */
						bool convergence_criteria_0 = (NM==NM_max);													// Cannot go past max NM
						bool convergence_criteria_1 = (NM-idx_best_PA>4);												// If nothing better is found in the last 3 PAs, we assume we found the best
						bool convergence_criteria_2 = (min_PA_diff<1e-6*std::abs(pade_approximants_BU_array[idx_NDOS*(NM_max+1) + idx_best_PA]));	// If the difference is less than the 4th significant digit, we assume "good enough"
						bool convergence_criteria_3 = (min_PA_diff<1e-7);												// If we are below single precision resolution, assume convergence
						
						if (convergence_criteria_0 ||
							convergence_criteria_1 ||
							convergence_criteria_2 ||
							convergence_criteria_3){
							pade_approximants_BU_conv_array[idx_NDOS] = true;
							pade_approximants_BU_idx_array[idx_NDOS]  = idx_best_PA;
							num_converged_elements += 1;
						}
					}
				}
			}
		}
		printf("         - Done \n"); fflush(stdout);
	}

	printf("     - Extracting on-shell U-matrix elements \n"); fflush(stdout);
	// [EN] After convergence selection, the final U elements are just the best Padé values written back into the
	// elastic and breakup storage layout expected by the rest of the code. / [CN] 在选定收敛阶数后，最终的 U 元素
	// 就是把最佳 Padé 值回写到程序其余部分期望的弹性/破裂存储布局中。
	/* Set on-shell elastic U-matrix elements equal "best" PA */
	for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
		for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
			for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
				const elastic_on_shell_index ndos = make_elastic_on_shell_index(idx_d_row,
																				 idx_d_col,
																				 idx_q_com,
																				 deuteron_idx_array,
																				 q_com_idx_array,
																				 num_deuteron_states,
																				 num_q_com,
																				 Nq_WP,
																				 Np_WP);

				size_t idx_best_PA = pade_approximants_idx_array[ndos.value_storage_idx];

				U_array[ndos.value_storage_idx] = pade_approximants_array[ndos.value_storage_idx*(NM_max+1) + idx_best_PA];
				printf("       - U-matrix element for alpha'=%ld, alpha=%ld, q=%ld: %.10e + %.10ei \n", ndos.alpha_row, ndos.alpha_col, ndos.q_idx, U_array[ndos.value_storage_idx].real(), U_array[ndos.value_storage_idx].imag());
			}
		}
	}
	/* Set on-shell breakup U-matrix elements equal "best" PA */
	if (run_parameters.include_breakup_channels){
		for (size_t idx_q_com=0; idx_q_com<num_q_com; idx_q_com++){
			int BU_chn_start = chn_os_indexing.q_com_BU_idx_array[idx_q_com];
			int BU_chn_end   = chn_os_indexing.q_com_BU_idx_array[idx_q_com+1];
			for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
				for (size_t idx_BU_chn=BU_chn_start; idx_BU_chn<BU_chn_end; idx_BU_chn++){
					/* Nucleon-deuteron breakup on-shell (NDOS) indices */
					size_t idx_alpha_NDOS = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 0];
					size_t idx_q_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 1];
					size_t idx_p_NDOS 	  = chn_os_indexing.alphapq_idx_array[idx_BU_chn*3 + 2];

					size_t idx_NDOS = breakup_value_storage_index(idx_d_row,
																 idx_BU_chn,
																 chn_os_indexing.num_BU_chns);

					size_t idx_best_PA = pade_approximants_BU_idx_array[idx_NDOS];

					U_BU_array[idx_NDOS] = pade_approximants_BU_array[idx_NDOS*(NM_max+1) + idx_best_PA];
				}
			}
		}
	}
	printf("       - Done \n"); fflush(stdout);

	/* Free allocated working space */
	delete [] a_coeff_array;
	delete [] a_BU_coeff_array;
	delete [] neumann_store_array;
	delete [] CPVC_cols_array;
	delete [] omp_CPVC_col_array;
	delete [] omp_CPVC_row_to_nnz_array;
	delete [] omp_CPVC_nnz_to_row_array;
	delete [] re_A_An_row_array;
	delete [] im_A_An_row_array;
	delete [] re_A_An_row_array_prev;
	delete [] im_A_An_row_array_prev;
	delete [] re_A_An_row_array_comp;
	delete [] im_A_An_row_array_comp;
	delete [] re_A_An_row_array_prod;
	delete [] im_A_An_row_array_prod;
	delete [] pade_approximants_array;
	delete [] pade_approximants_idx_array;
	delete [] pade_approximants_conv_array;
	delete [] pade_approximants_BU_array;
	delete [] pade_approximants_BU_idx_array;
	delete [] pade_approximants_BU_conv_array;
}

void solve_faddeev_equations(cdouble*  U_array,
							 cdouble*  U_BU_array,
							 cdouble*  G_array,
							 double*   P123_sparse_val_array,
							 int*      P123_sparse_row_array,
							 size_t*   P123_sparse_col_array_csc,
							 size_t    P123_sparse_dim,
							 double*   V_WP_unco_array,
							 double*   V_WP_coup_array,
							 swp_statespace swp_states,
							 fwp_statespace fwp_states,
							 channel_os_indexing chn_os_indexing,
							 pw_3N_statespace pw_states,
							 std::string file_identification,
					         run_params run_parameters,
							 const three_nucleon_force_model* tnf){
	
	/* Make local pointers & variables for on-shell channel-indexing */
	int*   q_com_idx_array		= chn_os_indexing.q_com_idx_array;
	int*   deuteron_idx_array	= chn_os_indexing.deuteron_idx_array;
	size_t num_q_com			= (size_t) chn_os_indexing.num_T_lab;
	size_t num_deuteron_states	= (size_t) chn_os_indexing.num_deuteron_states;
	/* Make local pointers & variables for SWP-statespace */
	size_t    Nq_WP					= (size_t) swp_states.Nq_WP;
	size_t    Np_WP					= (size_t) swp_states.Np_WP;
	double*   C_WP_unco_array		= swp_states.C_SWP_unco_array;
	double*   C_WP_coup_array		= swp_states.C_SWP_coup_array;
	int  	  num_2N_unco_states	= swp_states.num_2N_unco_states;
	int  	  num_2N_coup_states	= swp_states.num_2N_coup_states;
	/* Make local variable for pw-statespace */
	size_t Nalpha = pw_states.Nalpha;

	/* Test PVC- and CPVC-column multiplication routines with brute-force routines
	 * WARNING: VERY SLOW TEST, ONLY FOR BENCHMARKING */
	bool test_PVC_col_routine   = false;
	bool test_CPVC_col_routine  = false;

	/* Create C^T-product pointer-arrays in row-major format */
	double** CT_RM_array = new double* [Nalpha*Nalpha];
	create_CT_row_maj_3N_pointer_array(CT_RM_array,
									   C_WP_unco_array,
									   C_WP_coup_array,
									   num_2N_unco_states,
									   num_2N_coup_states,
									   Np_WP,
									   pw_states,
									   run_parameters);

	/* Create VC-product pointer-arrays in column-major format */
	double** VC_CM_array  = new double* [Nalpha*Nalpha];
	create_VC_col_maj_3N_pointer_array(VC_CM_array,
									   C_WP_unco_array,
									   C_WP_coup_array,
									   V_WP_unco_array,
									   V_WP_coup_array,
									   num_2N_unco_states,
									   num_2N_coup_states,
									   Np_WP,
									   pw_states,
									   run_parameters);
	
	size_t  dense_dim = Nalpha * Nq_WP * Np_WP;

	// [EN] Bundle 3NF context for the column hot-path. When tnf->enabled()==false (null object) the
	// context still exists but the hot-path 3NF branch is skipped via a single test.
	// [CN] 把 3NF 上下文打包给列计算热路径。当 tnf->enabled()==false 时上下文仍存在但热路径通过一次测试跳过。
	tnf_kernel_context tnf_ctx;
	tnf_ctx.tnf        = tnf;
	tnf_ctx.pw_states  = &pw_states;
	tnf_ctx.p_WP_array = fwp_states.p_WP_array;
	tnf_ctx.q_WP_array = fwp_states.q_WP_array;
	tnf_ctx.CT_RM_array = CT_RM_array;
	tnf_ctx.w1_scale   = run_parameters.w1_scale;
	tnf_ctx.w1_cache   = nullptr;

	// [EN] Pre-evaluate W^(1) on the WP bin-midpoint grid so the CPVC hot path
	// becomes O(1) lookups instead of full PW recoupling + GL quadrature per call.
	// Skipped when 3NF is null/disabled or w1_scale=0 (the hot-path early-exit guards
	// the same conditions, so cache build would just waste memory).
	// / [CN] 提前在 WP bin 中点网格上算好 W^(1)，让 CPVC 热路径变成 O(1) 查表。
	W1_PW_cache w1_cache_storage;
	if (tnf != nullptr && tnf->enabled() && run_parameters.w1_scale != 0.0) {
		printf(" - Pre-building W^(1) PW cache (Np=%zu, Nq=%zu, Nalpha=%d) ... \n",
			   Np_WP, Nq_WP, pw_states.Nalpha); fflush(stdout);
		auto t0 = std::chrono::system_clock::now();
		w1_cache_storage.build(*tnf,
							   fwp_states.p_WP_array, Np_WP,
							   fwp_states.q_WP_array, Nq_WP,
							   pw_states,
							   run_parameters);
		auto t1 = std::chrono::system_clock::now();
		printf("   - Done. blocks=%zu, size=%.2f GB, build=%.1f s \n",
			   w1_cache_storage.num_blocks(),
			   w1_cache_storage.total_bytes() / 1.073741824e9,
			   std::chrono::duration<double>(t1 - t0).count()); fflush(stdout);
		tnf_ctx.w1_cache = &w1_cache_storage;
	}

	/* Test optimized routine for PVC columns */
	if (test_PVC_col_routine){
		printf("   - Testing PVC-column routine ... \n");
		PVC_col_calc_test(Nalpha,
						  Nq_WP,
						  Np_WP,
						  VC_CM_array,
						  P123_sparse_val_array,
						  P123_sparse_row_array,
						  P123_sparse_col_array_csc,
						  P123_sparse_dim);
		printf("     - Done \n");
	}
	/* Test optimized routine for CPVC columns */
	if (test_CPVC_col_routine){
		printf("   - Testing CPVC-column routine ... \n");
		CPVC_col_calc_test(Nalpha,
						   Nq_WP,
						   Np_WP,
						   CT_RM_array,
						   VC_CM_array,
						   P123_sparse_val_array,
						   P123_sparse_row_array,
						   P123_sparse_col_array_csc,
						   P123_sparse_dim);
		printf("     - Done \n");
	}
	
	auto timestamp_solve_start = std::chrono::system_clock::now();
	if (run_parameters.solve_dense==false){
		pade_method_solve(U_array,
						  U_BU_array,
						  G_array,
						  q_com_idx_array,	  num_q_com,
						  deuteron_idx_array, num_deuteron_states,
						  Nalpha,
						  Nq_WP,
						  Np_WP,
						  CT_RM_array,
						  VC_CM_array,
						  P123_sparse_val_array,
						  P123_sparse_row_array,
						  P123_sparse_col_array_csc,
						  P123_sparse_dim,
						  chn_os_indexing,
					      run_parameters,
						  file_identification,
						  tnf_ctx);
	}
	else{
		/* Solve the Faddeev eq. using a dense MKL-solver.
	 	 * Obviously, this only works for small systems and is meant for benchmarking only */
		printf("     - Solving Faddeev equation using a dense direct solver (WARNING: CAN TAKE LONG) ... \n");
		faddeev_dense_solver(U_array,
						     G_array,
						     q_com_idx_array,	 num_q_com,
						     deuteron_idx_array, num_deuteron_states,
						     Nalpha,
						     Nq_WP,
						     Np_WP,
						     CT_RM_array,
						     VC_CM_array,
						     P123_sparse_val_array,
						     P123_sparse_row_array,
						     P123_sparse_col_array_csc,
						     P123_sparse_dim,
						     tnf_ctx);
	}

	auto timestamp_solve_end = std::chrono::system_clock::now();
	std::chrono::duration<double> time_solve = timestamp_solve_end - timestamp_solve_start;
	printf("     - Done. Time used: %.6f\n", time_solve.count());
}

/* Create array of pointers to C^T matrices for product (C^T)PVC in row-major format */
void create_CT_row_maj_3N_pointer_array(double** CT_RM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										int  	 num_2N_unco_states,
										int  	 num_2N_coup_states,
										size_t   Np_WP,
										pw_3N_statespace pw_states,
										run_params run_parameters){
	
	size_t Nalpha		  = pw_states.Nalpha;
	int*   L_2N_array	  = pw_states.L_2N_array;
	int*   S_2N_array	  = pw_states.S_2N_array;
	int*   J_2N_array	  = pw_states.J_2N_array;
	int*   T_2N_array	  = pw_states.T_2N_array;
	int*   L_1N_array	  = pw_states.L_1N_array;
	int*   two_J_1N_array = pw_states.two_J_1N_array;
	int*   two_T_3N_array = pw_states.two_T_3N_array;

	/* This test will be reused several times */
	bool tensor_force_true = (run_parameters.tensor_force==true);
	
	/* Number of uncoupled and coupled 2N-channels */
	int num_unco_chns = num_2N_unco_states;
	int num_coup_chns = num_2N_coup_states;

	double* C_subarray  = NULL;
	double* CT_subarray = NULL;

	double* CT_unco_array = new double [Np_WP*Np_WP   * num_unco_chns];
	double* CT_coup_array = new double [Np_WP*Np_WP*4 * num_coup_chns];
	
	/* Copy and transpose all 2N-uncoupled C-arrays */
	for (int idx_chn_unco=0; idx_chn_unco<num_unco_chns; idx_chn_unco++){
		size_t idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
			
		C_subarray  = &C_WP_unco_array[idx_2N_mat_WP_unco];
		CT_subarray = &CT_unco_array  [idx_2N_mat_WP_unco];

		/* Copy content to avoid rewriting C-arrays */
		std::copy(C_subarray, C_subarray + Np_WP*Np_WP, CT_subarray);
			
		/* Transpose C to get C^T */
		simple_transpose_matrix_routine(CT_subarray, Np_WP);
	}

	/* Copy and transpose all 2N-coupled C-arrays */
	for (int idx_chn_coup=0; idx_chn_coup<num_coup_chns; idx_chn_coup++){
		size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;

		C_subarray  = &C_WP_coup_array[idx_2N_mat_WP_coup];
		CT_subarray = &CT_coup_array  [idx_2N_mat_WP_coup];
		
		/* Copy content to avoid rewriting C-arrays */
		std::copy(C_subarray, C_subarray + 4*Np_WP*Np_WP, CT_subarray);
		
		/* Transpose C to get C^T */
		simple_transpose_matrix_routine(CT_subarray, 2*Np_WP);

		/* Restucture coupled matrix array into 4 seperate arrays */
		restructure_coupled_VC_product(CT_subarray, Np_WP);
	}
	
	/* Row state */
	for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		int L_2N_r 	   = L_2N_array[idx_alpha_r];
		int S_2N_r 	   = S_2N_array[idx_alpha_r];
		int J_2N_r 	   = J_2N_array[idx_alpha_r];
		int T_2N_r 	   = T_2N_array[idx_alpha_r];
		int L_1N_r 	   = L_1N_array[idx_alpha_r];
		int two_J_1N_r = two_J_1N_array[idx_alpha_r];
		int two_T_3N_r = two_T_3N_array[idx_alpha_r];

		/* Column state */
		for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_2N_c 	   = L_2N_array[idx_alpha_c];
			int S_2N_c 	   = S_2N_array[idx_alpha_c];
			int J_2N_c 	   = J_2N_array[idx_alpha_c];
			int T_2N_c 	   = T_2N_array[idx_alpha_c];
			int L_1N_c 	   = L_1N_array[idx_alpha_c];
			int two_J_1N_c = two_J_1N_array[idx_alpha_c];
			int two_T_3N_c = two_T_3N_array[idx_alpha_c];

			/* Check if possible channel through interaction */
			bool check_T = (T_2N_r==T_2N_c);
			bool check_J = (J_2N_r==J_2N_c);
			bool check_S = (S_2N_r==S_2N_c);
			bool check_L = ( (tensor_force_true && abs(L_2N_r-L_2N_c)<=2) || L_2N_r==L_2N_c);
			bool check_l = (L_1N_r==L_1N_c);
			bool check_j = (two_J_1N_r==two_J_1N_c);

			/* Check if possible channel through interaction */
			if (check_T && check_J && check_S && check_L && check_l && check_j){

				/* Detemine if this is a coupled channel.
				 * !!! With isospin symmetry-breaking we count 1S0 as a coupled matrix via T_3N-coupling !!! */
				bool coupled_matrix = false;
				bool state_1S0 = (S_2N_r==0 && J_2N_r==0 && L_2N_r==0);
				bool coupled_via_L_2N = (tensor_force_true && (L_2N_r!=L_2N_c || (L_2N_r==L_2N_c && L_2N_r!=J_2N_r && J_2N_r!=0)));
				bool coupled_via_T_3N = (state_1S0==true && run_parameters.isospin_breaking_1S0==true);
				if (coupled_via_L_2N && coupled_via_T_3N){
					raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
				}
				if (coupled_via_L_2N || coupled_via_T_3N){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}

				/* find which VC-product corresponds to the current coupling */
				if (coupled_matrix){
					size_t idx_chn_coup       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, coupled_matrix, run_parameters);
					size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;
					if ( (coupled_via_L_2N && L_2N_r<L_2N_c) || (coupled_via_T_3N && two_T_3N_r<two_T_3N_c) ){       // L_r=J_r-1, L_c=J_r+1 OR (for 1S0) two_T_3N_r==1/2, two_T_3N_c==3/2
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 1*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r>L_2N_c) || (coupled_via_T_3N && two_T_3N_r>two_T_3N_c) ){  // L_r=J_r+1, L_c=J_r-1 OR (for 1S0) two_T_3N_r==3/2, two_T_3N_c==1/2
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 2*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r<J_2N_c) || (coupled_via_T_3N && two_T_3N_r==1) ){ 		 // L_r=J_r-1, L_c=J_r-1 OR (for 1S0) two_T_3N_r==1/2, two_T_3N_c==1/2
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 0*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r>J_2N_c) || (coupled_via_T_3N && two_T_3N_r==3) ){  		 // L_r=J_r+1, L_c=J_r+1 OR (for 1S0) two_T_3N_r==3/2, two_T_3N_c==3/2
						CT_subarray = &CT_coup_array[idx_2N_mat_WP_coup + 3*Np_WP*Np_WP];
					}
					else{
						raise_error("Unknown 2N coupled-matrix requested in CT-array setup.");
					}
				}
				else{
					size_t idx_chn_unco       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, coupled_matrix, run_parameters);
					size_t idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
					CT_subarray = &CT_unco_array[idx_2N_mat_WP_unco];
				}
			}
			else{
				CT_subarray = NULL;
			}

			/* Unique index for two given states alpha */
			size_t idx_CT_RM = idx_alpha_r*Nalpha + idx_alpha_c;
			CT_RM_array[idx_CT_RM] = CT_subarray;
		}
	}
	
}

/* Restructures NN coupled matrix as 4 seperate matrices
 * !!! WARNING: COLUMN-MAJOR ALGORITHM !!! */
void restructure_coupled_VC_product(double* VC_product, size_t Np_WP){
	double* VC_product_temp = new double [4*Np_WP*Np_WP];

	for (size_t col_block=0; col_block<2; col_block++){
		for (size_t row_block=0; row_block<2; row_block++){

			size_t idx_block = col_block*2*Np_WP*Np_WP + row_block*Np_WP*Np_WP;

			for (size_t col=0; col<Np_WP; col++){
				for (size_t row=0; row<Np_WP; row++){
					size_t prestructure_col_idx = (col + Np_WP*col_block)*2*Np_WP;
					size_t prestructure_row_idx =  row + Np_WP*row_block;

					VC_product_temp[idx_block + col*Np_WP+row] = VC_product[prestructure_col_idx + prestructure_row_idx];
				}
			}

		}
	}

	/* Copy from temp array */
	std::copy(VC_product_temp, VC_product_temp + 4*Np_WP*Np_WP, VC_product);

	/* Delete temporary array */
	delete [] VC_product_temp;
}

/* Create array of pointers to VC-product matrices for product (C^T)PVC in column-major format*/
void create_VC_col_maj_3N_pointer_array(double** VC_CM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										double*  V_WP_unco_array,
										double*  V_WP_coup_array,
										int  	 num_2N_unco_states,
										int  	 num_2N_coup_states,
										size_t   Np_WP,
										pw_3N_statespace pw_states,
										run_params run_parameters){
	
	size_t Nalpha		  = pw_states.Nalpha;
	int*   L_2N_array	  = pw_states.L_2N_array;
	int*   S_2N_array	  = pw_states.S_2N_array;
	int*   J_2N_array	  = pw_states.J_2N_array;
	int*   T_2N_array	  = pw_states.T_2N_array;
	int*   L_1N_array	  = pw_states.L_1N_array;
	int*   two_J_1N_array = pw_states.two_J_1N_array;
	int*   two_T_3N_array = pw_states.two_T_3N_array;
	
	/* This test will be reused several times */
	bool tensor_force_true = (run_parameters.tensor_force==true);

	/* Number of uncoupled and coupled 2N-channels */
	int num_unco_chns = num_2N_unco_states;
	int num_coup_chns = num_2N_coup_states;

	double* V_subarray = NULL;
	double* C_subarray = NULL;

	double* VC_unco_array = new double [Np_WP*Np_WP   * num_unco_chns];
	double* VC_coup_array = new double [Np_WP*Np_WP*4 * num_coup_chns];

	double* VC_product = NULL;

	/* Calculate all 2N-uncoupled VC-products and convert to column-major format */
	for (int idx_chn_unco=0; idx_chn_unco<num_unco_chns; idx_chn_unco++){
		size_t idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
			
		V_subarray = &V_WP_unco_array[idx_2N_mat_WP_unco];
		C_subarray = &C_WP_unco_array[idx_2N_mat_WP_unco];
		VC_product = &VC_unco_array  [idx_2N_mat_WP_unco];

		/* Multiply V and C using BLAS */
		dot_MM(V_subarray, C_subarray, VC_product, Np_WP, Np_WP, Np_WP);

		/* Transpose VC-product to get column-major format */
		simple_transpose_matrix_routine(VC_product, Np_WP);
	}

	/* Calculate all 2N-coupled VC-products and convert to column-major format */
	for (int idx_chn_coup=0; idx_chn_coup<num_coup_chns; idx_chn_coup++){
		size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;

		V_subarray = &V_WP_coup_array[idx_2N_mat_WP_coup];
		C_subarray = &C_WP_coup_array[idx_2N_mat_WP_coup];
		VC_product = &VC_coup_array  [idx_2N_mat_WP_coup];

		/* Multiply V and C using BLAS */
		dot_MM(V_subarray, C_subarray, VC_product, 2*Np_WP, 2*Np_WP, 2*Np_WP);
		/* Transpose VC-product to get column-major format */
		simple_transpose_matrix_routine(VC_product, 2*Np_WP);
		/* Restucture coupled matrix array into 4 seperate arrays */
		restructure_coupled_VC_product(VC_product, Np_WP);
	}

	/* Row state */
	for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		int L_2N_r     = L_2N_array[idx_alpha_r];
		int S_2N_r     = S_2N_array[idx_alpha_r];
		int J_2N_r     = J_2N_array[idx_alpha_r];
		int T_2N_r     = T_2N_array[idx_alpha_r];
		int L_1N_r 	   = L_1N_array[idx_alpha_r];
		int two_J_1N_r = two_J_1N_array[idx_alpha_r];
		int two_T_3N_r = two_T_3N_array[idx_alpha_r];

		/* Column state */
		for (size_t idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_2N_c     = L_2N_array[idx_alpha_c];
			int S_2N_c     = S_2N_array[idx_alpha_c];
			int J_2N_c     = J_2N_array[idx_alpha_c];
			int T_2N_c     = T_2N_array[idx_alpha_c];
			int L_1N_c 	   = L_1N_array[idx_alpha_c];
			int two_J_1N_c = two_J_1N_array[idx_alpha_c];
			int two_T_3N_c = two_T_3N_array[idx_alpha_c];

			/* Check if possible channel through interaction */
			bool check_T = (T_2N_r==T_2N_c);
			bool check_J = (J_2N_r==J_2N_c);
			bool check_S = (S_2N_r==S_2N_c);
			bool check_L = ( (tensor_force_true && abs(L_2N_r-L_2N_c)<=2) || L_2N_r==L_2N_c);
			bool check_l = (L_1N_r==L_1N_c);
			bool check_j = (two_J_1N_r==two_J_1N_c);

			/* Check if possible channel through interaction */
			if (check_T && check_J && check_S && check_L && check_l && check_j){

				/* Detemine if this is a coupled channel.
				 * !!! With isospin symmetry-breaking we count 1S0 as a coupled matrix via T_3N-coupling !!! */
				bool coupled_matrix = false;
				bool state_1S0 = (S_2N_r==0 && J_2N_r==0 && L_2N_r==0);
				bool coupled_via_L_2N = (tensor_force_true && (L_2N_r!=L_2N_c || (L_2N_r==L_2N_c && L_2N_r!=J_2N_r && J_2N_r!=0)));
				bool coupled_via_T_3N = (state_1S0==true && run_parameters.isospin_breaking_1S0==true);
				if (coupled_via_L_2N && coupled_via_T_3N){
					raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
				}
				if (coupled_via_L_2N || coupled_via_T_3N){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}

				/* find which VC-product corresponds to the current coupling */
				if (coupled_matrix){
					size_t idx_chn_coup       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, coupled_matrix, run_parameters);
					size_t idx_2N_mat_WP_coup = idx_chn_coup*4*Np_WP*Np_WP;
					if ( (coupled_via_L_2N && L_2N_r<L_2N_c) || (coupled_via_T_3N && two_T_3N_r<two_T_3N_c) ){       // L_r=J_r-1, L_c=J_r+1 OR (for 1S0) two_T_3N_r==1/2, two_T_3N_c==3/2
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 1*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r>L_2N_c) || (coupled_via_T_3N && two_T_3N_r>two_T_3N_c) ){  // L_r=J_r+1, L_c=J_r-1 OR (for 1S0) two_T_3N_r==3/2, two_T_3N_c==1/2
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 2*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r<J_2N_c) || (coupled_via_T_3N && two_T_3N_r==1) ){ 		 // L_r=J_r-1, L_c=J_r-1 OR (for 1S0) two_T_3N_r==1/2, two_T_3N_c==1/2
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 0*Np_WP*Np_WP];
					}
					else if ( (coupled_via_L_2N && L_2N_r>J_2N_c) || (coupled_via_T_3N && two_T_3N_r==3) ){  		 // L_r=J_r+1, L_c=J_r+1 OR (for 1S0) two_T_3N_r==3/2, two_T_3N_c==3/2
						VC_product = &VC_coup_array[idx_2N_mat_WP_coup + 3*Np_WP*Np_WP];
					}
					else{
						raise_error("Unknown 2N coupled-matrix requested in VC-array setup.");
					}
				}
				else{
					size_t idx_chn_unco       = (size_t) unique_2N_idx(L_2N_r, S_2N_r, J_2N_r, T_2N_r, coupled_matrix, run_parameters);
					size_t idx_2N_mat_WP_unco = idx_chn_unco*Np_WP*Np_WP;
					VC_product = &VC_unco_array[idx_2N_mat_WP_unco];
				}
			}
			else{
				VC_product = NULL;
			}

			/* Unique index for two given states alpha */
			size_t idx_VC_CM = idx_alpha_r*Nalpha + idx_alpha_c;
			VC_CM_array[idx_VC_CM] = VC_product;
		}
	}
}

void PVC_col_brute_force(double*  col_array,
					     size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
					     size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
					     double** VC_CM_array,
					     double*  P123_val_array,
					     int*  	  P123_row_array,
					     size_t*  P123_col_array,
					     size_t   P123_dim){
	
	bool   print_content = false;
	size_t idx1 = idx_alpha_c*Np_WP*Nq_WP + idx_q_c*Np_WP + idx_p_c;
	double* VC_ptr     = NULL;

	for (size_t idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		for (size_t idx_q_r=0; idx_q_r<Nq_WP; idx_q_r++){
			for (size_t idx_p_r=0; idx_p_r<Np_WP; idx_p_r++){
				size_t idx_P123_row = idx_alpha_r*Nq_WP*Np_WP + idx_q_r*Np_WP + idx_p_r;

				double sum = 0; 
				for (size_t idx_alpha_j=0; idx_alpha_j<Nalpha; idx_alpha_j++){
					VC_ptr = VC_CM_array[idx_alpha_c*Nalpha + idx_alpha_j];
					if (VC_ptr!=NULL){
						for (size_t idx_q_j=0; idx_q_j<Nq_WP; idx_q_j++){
							if (idx_q_j == idx_q_c){
								for (size_t idx_p_j=0; idx_p_j<Np_WP; idx_p_j++){
									if (print_content){std::cout << "  - Index j " << idx_alpha_j << " " << idx_q_j << " " << idx_p_j << std::endl;}

									size_t idx_P123_col = idx_alpha_j*Nq_WP*Np_WP + idx_q_j*Np_WP + idx_p_j;

									size_t idx_P123_row_lower = P123_col_array[idx_P123_col];
									size_t idx_P123_row_upper = P123_col_array[idx_P123_col+1];
									bool idx_found = false;
									size_t idx_P123_val = 0;
									for (size_t nnz_idx=idx_P123_row_lower; nnz_idx<idx_P123_row_upper; nnz_idx++){
										if (P123_row_array[nnz_idx] == idx_P123_row){
											idx_found = true;
											idx_P123_val = nnz_idx;
										}
									}

									if (idx_found){
										double P_element  = P123_val_array[idx_P123_val];
										double VC_element = VC_ptr[idx_p_c*Np_WP + idx_p_j];

										sum += P_element * VC_element;
									}
								}
							}
						}
					}
				}

				size_t idx_row = idx_alpha_r*Np_WP*Nq_WP + idx_q_r*Np_WP + idx_p_r;
				/* Write inner product to col_array of PVC-product */
				col_array[idx_row] = sum;
			}
		}
	}
}

void CPVC_col_brute_force(double*  col_array,
						  size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
						  size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
						  double** CT_RM_array,
						  double** VC_CM_array,
						  double*  P123_val_array,
						  int*     P123_row_array,
						  size_t*  P123_col_array,
						  size_t   P123_dim){
	double* CT_ptr = NULL;
	double* VC_ptr = NULL;
}
