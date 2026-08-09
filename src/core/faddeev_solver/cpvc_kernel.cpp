#include "cpvc_kernel.h"
#include "constants.h"
#include "interactions/three_nucleon_force_model.h"
#include "interactions/w1_pw_cache.h"
#include <cmath>
#include <cstdio>
#include <vector>

// [EN] Self-contained AGS kernel-algebra builders shared by the production
// solver and the independent finite-dimensional oracle.  Keeping the complete
// W^(1)·(1+P)·C application here prevents the dense-column and on-shell-row
// paths from drifting apart. / [CN] 自包含的 AGS 核代数构造函数，由生产求解器与
// 独立有限维 oracle 共用；完整 W^(1)·(1+P)·C 作用集中在此，避免稠密列路径与
// on-shell 行路径出现代数分叉。

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

void add_W1_one_plus_P_C_col(double*  col_array,
							 size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
							 size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
							 double** CT_RM_array,
							 double*  P123_val_array,
							 int*     P123_row_array,
							 size_t*  P123_col_array,
							 size_t   P123_dim,
							 const tnf_kernel_context& tnf_ctx){

	if (tnf_ctx.tnf == nullptr || !tnf_ctx.tnf->enabled() || tnf_ctx.w1_scale == 0.0){
		return;
	}

	const three_nucleon_force_model* tnf = tnf_ctx.tnf;
	const pw_3N_statespace& pw_st = *tnf_ctx.pw_states;
	const double* p_WP = tnf_ctx.p_WP_array;
	const double* q_WP = tnf_ctx.q_WP_array;

	const double inv_hbarc  = 1.0 / hbarc;
	const double inv_hbarc5 = inv_hbarc * inv_hbarc * inv_hbarc * inv_hbarc * inv_hbarc;
	const double w1_unit    = inv_hbarc5 * tnf_ctx.w1_scale;

	const double q_c_mid = 0.5 * (q_WP[idx_q_c] + q_WP[idx_q_c + 1]);
	const double dq_c    = q_WP[idx_q_c + 1] - q_WP[idx_q_c];
	const double wq_c    = q_c_mid * std::sqrt(dq_c);
	const double q_c_fm  = q_c_mid * inv_hbarc;

	// Identity part: W^(1)·C.  C acts only on the pair line, so q_j=q_c,
	// but a tensor-coupled C has off-diagonal alpha blocks and requires the
	// complete intermediate-alpha contraction.
	for (size_t idx_alpha_j = 0; idx_alpha_j < Nalpha; idx_alpha_j++){
		double* C_block = CT_RM_array[idx_alpha_c * Nalpha + idx_alpha_j];
		if (C_block == nullptr) continue;

		for (size_t idx_alpha_r = 0; idx_alpha_r < Nalpha; idx_alpha_r++){
			// These guards belong to W1(alpha_r,alpha_j), not to the external
			// column alpha_c.  The C block independently enforces allowed pair
			// mixing between alpha_j and alpha_c.
			if (pw_st.two_J_3N_array[idx_alpha_r] != pw_st.two_J_3N_array[idx_alpha_j]) continue;
			if (pw_st.two_T_3N_array[idx_alpha_r] != pw_st.two_T_3N_array[idx_alpha_j]) continue;
			if (pw_st.P_3N_array[idx_alpha_r]     != pw_st.P_3N_array[idx_alpha_j])     continue;

			for (size_t idx_q_r = 0; idx_q_r < Nq_WP; idx_q_r++){
				const double q_r_mid = 0.5 * (q_WP[idx_q_r] + q_WP[idx_q_r + 1]);
				const double dq_r    = q_WP[idx_q_r + 1] - q_WP[idx_q_r];
				const double wq_r    = q_r_mid * std::sqrt(dq_r);
				const double q_r_fm  = q_r_mid * inv_hbarc;

				for (size_t idx_p_r = 0; idx_p_r < Np_WP; idx_p_r++){
					const double p_r_mid = 0.5 * (p_WP[idx_p_r] + p_WP[idx_p_r + 1]);
					const double dp_r    = p_WP[idx_p_r + 1] - p_WP[idx_p_r];
					const double wp_r    = p_r_mid * std::sqrt(dp_r);
					const double p_r_fm  = p_r_mid * inv_hbarc;
					double w1c = 0.0;

					for (size_t idx_p_j = 0; idx_p_j < Np_WP; idx_p_j++){
						// C_(alpha_j p_j,alpha_c p_c)
						//   = (C^T)_(alpha_c p_c,alpha_j p_j).
						const double C_val = C_block[idx_p_c * Np_WP + idx_p_j];
						if (C_val == 0.0) continue;

						const double p_j_mid = 0.5 * (p_WP[idx_p_j] + p_WP[idx_p_j + 1]);
						const double dp_j    = p_WP[idx_p_j + 1] - p_WP[idx_p_j];
						const double wp_j    = p_j_mid * std::sqrt(dp_j);
						const double p_j_fm  = p_j_mid * inv_hbarc;
						double w1_bin;
						if (tnf_ctx.w1_cache != nullptr){
							w1_bin = tnf_ctx.w1_cache->get(idx_alpha_r, idx_alpha_j,
														   idx_p_r, idx_q_r, idx_p_j, idx_q_c)
							       * tnf_ctx.w1_scale;
						}
						else{
							const double w1_raw = tnf->W1_element(idx_alpha_r, idx_alpha_j,
																p_r_fm, q_r_fm,
																p_j_fm, q_c_fm, pw_st);
							w1_bin = (w1_raw * w1_unit) * (wp_r * wq_r * wp_j * wq_c);
						}
						w1c += w1_bin * C_val;
					}

					if (w1c != 0.0){
						const size_t idx_row = idx_alpha_r * Nq_WP * Np_WP
											 + idx_q_r * Np_WP + idx_p_r;
						col_array[idx_row] += w1c;
					}
				}
			}
		}
	}

	// Permutation part: W^(1)·P·C.  calculate_PVC_col treats CT_RM as
	// column-major C, which is the same memory identity used above.
	std::vector<double> PC_col(Nalpha * Nq_WP * Np_WP, 0.0);
	calculate_PVC_col(PC_col.data(),
					  idx_alpha_c, idx_p_c, idx_q_c,
					  Nalpha, Nq_WP, Np_WP,
					  CT_RM_array,
					  P123_val_array,
					  P123_row_array,
					  P123_col_array,
					  P123_dim);

	for (size_t idx_k = 0; idx_k < PC_col.size(); idx_k++){
		const double pc_val = PC_col[idx_k];
		if (pc_val == 0.0) continue;

		const size_t idx_alpha_k = idx_k / (Nq_WP * Np_WP);
		const size_t idx_q_k = (idx_k % (Nq_WP * Np_WP)) / Np_WP;
		const size_t idx_p_k = idx_k % Np_WP;
		const double p_k_mid = 0.5 * (p_WP[idx_p_k] + p_WP[idx_p_k + 1]);
		const double q_k_mid = 0.5 * (q_WP[idx_q_k] + q_WP[idx_q_k + 1]);
		const double dp_k    = p_WP[idx_p_k + 1] - p_WP[idx_p_k];
		const double dq_k    = q_WP[idx_q_k + 1] - q_WP[idx_q_k];
		const double wp_k    = p_k_mid * std::sqrt(dp_k);
		const double wq_k    = q_k_mid * std::sqrt(dq_k);
		const double p_k_fm  = p_k_mid * inv_hbarc;
		const double q_k_fm  = q_k_mid * inv_hbarc;

		for (size_t idx_alpha_r = 0; idx_alpha_r < Nalpha; idx_alpha_r++){
			if (pw_st.two_J_3N_array[idx_alpha_r] != pw_st.two_J_3N_array[idx_alpha_k]) continue;
			if (pw_st.two_T_3N_array[idx_alpha_r] != pw_st.two_T_3N_array[idx_alpha_k]) continue;
			if (pw_st.P_3N_array[idx_alpha_r]     != pw_st.P_3N_array[idx_alpha_k])     continue;

			for (size_t idx_q_r = 0; idx_q_r < Nq_WP; idx_q_r++){
				const double q_r_mid = 0.5 * (q_WP[idx_q_r] + q_WP[idx_q_r + 1]);
				const double dq_r    = q_WP[idx_q_r + 1] - q_WP[idx_q_r];
				const double wq_r    = q_r_mid * std::sqrt(dq_r);
				const double q_r_fm  = q_r_mid * inv_hbarc;

				for (size_t idx_p_r = 0; idx_p_r < Np_WP; idx_p_r++){
					const double p_r_mid = 0.5 * (p_WP[idx_p_r] + p_WP[idx_p_r + 1]);
					const double dp_r    = p_WP[idx_p_r + 1] - p_WP[idx_p_r];
					const double wp_r    = p_r_mid * std::sqrt(dp_r);
					const double p_r_fm  = p_r_mid * inv_hbarc;
					double w1_wp;
					if (tnf_ctx.w1_cache != nullptr){
						w1_wp = tnf_ctx.w1_cache->get(idx_alpha_r, idx_alpha_k,
														  idx_p_r, idx_q_r, idx_p_k, idx_q_k)
						        * tnf_ctx.w1_scale;
					}
					else{
						const double w1_raw = tnf->W1_element(idx_alpha_r, idx_alpha_k,
																p_r_fm, q_r_fm,
																p_k_fm, q_k_fm, pw_st);
						w1_wp = (w1_raw * w1_unit) * (wp_r * wq_r * wp_k * wq_k);
					}
					if (w1_wp != 0.0){
						const size_t idx_row = idx_alpha_r * Nq_WP * Np_WP
											 + idx_q_r * Np_WP + idx_p_r;
						col_array[idx_row] += w1_wp * pc_val;
					}
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

	// The shared helper is also called by calculate_all_CPVC_rows, so the
	// dense-column and Padé/Neumann row paths use identical W1 contractions.
	add_W1_one_plus_P_C_col(PVC_col.data(),
							 idx_alpha_c, idx_p_c, idx_q_c,
							 Nalpha, Nq_WP, Np_WP,
							 CT_RM_array,
							 P123_val_array,
							 P123_row_array,
							 P123_col_array,
							 P123_dim,
							 tnf_ctx);
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
