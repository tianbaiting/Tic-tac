#ifndef CPVC_KERNEL_H
#define CPVC_KERNEL_H

#include <cstddef>

#include "type_defs.h"
#include "tnf_kernel_context.h"

// [EN] Self-contained AGS kernel-algebra builders, extracted from
// solve_faddeev.cpp so they can be unit-tested in isolation
// (tests/cpp/test_3nf_operator_oracle.cpp, Phase 0 of fix/3nf-physics-contract)
// and so the dense / sparse / row / column / Pade solver paths share ONE
// implementation of the algebra
//
//     A = C^T · ( P·V  +  (1 + P)·W^(1) ) · C
//
// (see docs/three_nf_equation_contract.md §4.2 / §6).
//
// Mathematical aliases (for navigation; historical names kept for WPCD-literature
// correspondence):
//   cpvc_kernel            == the packet-space AGS kernel module  ("AgsKernel")
//   calculate_CPVC_col     == build one driving column of A        ("build_driving_column")
//   add_one_plus_P_W1_C_col== the (1+P)·W^(1)·C 3NF part of A
//   calculate_PVC_col      == the P·V·C 2NF part of A
//
// These functions are deliberately free of HDF5 / 2NF-model / solver-loop
// dependencies: they only consume the pre-built WP arrays (CT_RM, VC_CM,
// P123 sparse) and the three_nucleon_force_model abstract interface.
//
// / [CN] 自包含的 AGS 核代数构造函数，从 solve_faddeev.cpp 抽出以便单独单元测试
// （Phase 0 oracle），并使 dense/sparse/row/column/Pade 路径共享同一实现。
// 仅依赖 WP 数组与 three_nucleon_force_model 抽象接口，无 HDF5/2NF/求解器循环依赖。

// Raw driving column PVC = P·(V·C) that appears in the AGS kernel before the
// left basis rotation by C^T. Applying the sparse permutation first keeps the
// expensive part of the kernel sparse.
void calculate_PVC_col(double*  col_array,
					   size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
					   size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
					   double** VC_CM_array,
					   double*  P123_val_array,
					   int*  	P123_row_array,
					   size_t*  P123_col_array,
					   size_t   P123_dim);

// Add the bare-3NF part (1+P)·W^(1)·C for one external packet column to a
// pre-existing right-kernel column.  This is the single implementation used by
// both calculate_CPVC_col and calculate_all_CPVC_rows.
//
// CT_RM storage contract:
//   CT_RM[a*Nalpha+b][i*Np+j]
//     = (C^T)_(a i,b j) = C_(b j,a i).
// Hence the identity term W1*C must sum alpha_j and use
//   C_(alpha_j p_j,alpha_c p_c)
//     = CT_RM[alpha_c*Nalpha+alpha_j][p_c*Np+p_j].
void add_one_plus_P_W1_C_col(double*  col_array,
							 size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
							 size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
							 double** CT_RM_array,
							 double*  P123_val_array,
							 int*     P123_row_array,
							 size_t*  P123_col_array,
							 size_t   P123_dim,
							 const tnf_kernel_context& tnf_ctx);

// CPVC = C^T·P·V·C is the packet-space kernel that drives both the first Neumann
// term and every later rescattering step. With 3NF enabled (tnf_ctx.tnf != null
// and tnf->enabled()), the column also accumulates (1+P)·W^(1)·C, i.e. the full
// AGS kernel column  A[:, col] = [C^T·(P·V + (1+P)·W^(1))·C][:, col].
//
// Operator ordering (locked, see docs/three_nf_equation_contract.md §3):
// (1+P) on the LEFT of the bare spectator component W^(1).  This follows by
// reducing the symmetrized AGS equation of Deltuva, Phys. Rev. C 80, 064002
// (2009), Eq. (7a), with tG0=vG and G0(1+tG0)=G.
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
					    const tnf_kernel_context& tnf_ctx);

#endif // CPVC_KERNEL_H
