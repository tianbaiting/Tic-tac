#ifndef TICTAC_CORE_REFERENCE_BRUTE_FORCE_H
#define TICTAC_CORE_REFERENCE_BRUTE_FORCE_H

// [EN] Reference / diagnostic builders for the packet-space AGS kernel.
//
// These O(dense_dim^2) brute-force routines and their *_calc_test harnesses are
// NOT on the production hot path. They exist to cross-check the optimized
// sparse column/row builders of cpvc_kernel against a naive reference, gated by
// the run-time flags `test_PVC_col_routine` / `test_CPVC_col_routine` inside
// solve_faddeev_equations. They are kept in their own header so the production
// API (solve_faddeev.h) reads as the clean entry surface, separate from
// reference/diagnostic utilities (task Phase 5).
//
// The mathematical object under test is unchanged:
//     A = C^T · ( P·V + (1 + P)·W^(1) ) · C   (docs/three_nf_equation_contract.md §4.2)
// / [CN] AGS 核的参考/诊断构造函数：O(dense_dim^2) 暴力例程与 *_calc_test 对照，
// 仅在 test_*_routine 旗标打开时运行，用于核验 cpvc_kernel 的稀疏优化路径。不在生产热路径上。

#include <cstddef>

void PVC_col_brute_force(double*  col_array,
						size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
						size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
						double** VC_CM_array,
						double*  P123_val_array,
						size_t*  P123_row_array,
						int*     P123_col_array,
						size_t   P123_dim);

void CPVC_col_brute_force(double*  col_array,
						  size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
						  size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
						  double** CT_RM_array,
						  double** VC_CM_array,
						  double*  P123_val_array,
						  size_t*  P123_row_array,
						  int*     P123_col_array,
						  size_t   P123_dim);

void PVC_col_calc_test(size_t   Nalpha,
					   size_t 	Nq_WP,
					   size_t 	Np_WP,
					   double** VC_CM_array,
					   double*  P123_sparse_val_array,
					   int*     P123_sparse_row_array,
					   size_t*  P123_sparse_col_array,
					   size_t   P123_sparse_dim);

void CPVC_col_calc_test(size_t   Nalpha,
						size_t 	 Nq_WP,
						size_t 	 Np_WP,
						double** CT_RM_array,
						double** VC_CM_array,
						double*  P123_sparse_val_array,
						int*     P123_sparse_row_array,
						size_t*  P123_sparse_col_array,
						size_t   P123_sparse_dim);

#endif // TICTAC_CORE_REFERENCE_BRUTE_FORCE_H
