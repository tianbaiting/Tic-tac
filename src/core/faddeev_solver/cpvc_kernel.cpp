#include "cpvc_kernel.h"
#include "constants.h"
#include "interactions/three_nucleon_force_model.h"
#include "interactions/w1_pw_cache.h"
#include <cmath>
#include <cstdio>
#include <vector>

// [EN] Self-contained AGS kernel-algebra builders, moved verbatim from
// solve_faddeev.cpp (lines 165-493 pre-extraction) so the Phase 0 oracle
// (tests/cpp/test_3nf_operator_oracle.cpp) can link against the production
// kernel without pulling in the full solver / HDF5 / 2NF-model object graph.
// Behaviour is byte-for-byte identical to the pre-extraction code: the only
// change is the translation unit boundary. / [CN] 从 solve_faddeev.cpp 原样
// 抽出的自包含 AGS 核代数构造函数，便于 Phase 0 oracle 链接而不引入完整求解器依赖。
// 行为与抽取前逐字节一致，仅翻译单元边界变化。

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

	// [EN] 3NF contribution (Born + iteration kernel): add W^(1)·(1+P)·C to PVC_col.
	//
	// OPERATOR ORDERING (locked, see docs/treatise/chapters/15_3nf_physics.tex
	// §operator-ordering and tests/cpp/test_faddeev_operator_order.cpp):
	//
	// The AGS kernel with 3NF (Witała 2008 PRC 77 034004 eq. 3) is
	//     K_AGS = P·V + W^(1)·(1 + P)
	// with W^(1) on the LEFT and (1+P) on the RIGHT. This is NOT the same as
	// (1+P)·W^(1) in the kernel, even though the two coincide as matrix elements
	// between fully antisymmetric states. The kernel acts in spectator-1 Faddeev
	// space where W^(1) is well-defined; left-multiplying by P would move the
	// operator out of its definition frame.
	//
	// Code structure (matches the algebra above):
	//   Identity part:    W^(1)·C  — direct sum over intermediate p_j with C-weight
	//   Permutation part: W^(1)·P·C — sparse P applied to C column first, then W^(1)
	//
	// / [CN] 3NF 贡献（Born + 迭代核）：把 W^(1)·(1+P)·C 加到 PVC_col 缓冲区。
	// 算符顺序锁定：W^(1) 在左，(1+P) 在右（Witała 2008 AGS 形式）。
	//   单位部分：W^(1)·C —— 对中间 p_j 直接做 W1×C 求和
	//   置换部分：W^(1)·P·C —— 用稀疏 P 先作用 C 列，再施加 W^(1)
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
