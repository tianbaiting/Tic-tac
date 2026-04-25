// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_permutation_matrix.cpp
// 行号区段：1..39
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================

#include "make_permutation_matrix.h"

namespace {

// [EN] P123 couples two packet states |X_{alpha,qp,pp}> and |X_{alpha',q,p}>. These helpers keep the four-packet
// indexing explicit so the implementation tracks the matrix element <X'|P|X> from chapter 6 and from the Miller
// WPCD formulation rather than burying it in ad hoc stride arithmetic. / [CN] P123 耦合的是两个波包态
// |X_{alpha,qp,pp}> 与 |X_{alpha',q,p}>；这些 helper 把四个波包指标显式写出来，使实现更直接对应第 6 章和 Miller
// 的 WPCD 形式中的矩阵元 <X'|P|X>，而不是把它淹没在零散 stride 运算里。
inline size_t packet_state_index(size_t alpha_idx,
								 size_t q_idx,
								 size_t p_idx,
								 size_t Nq_WP,
								 size_t Np_WP){
	return alpha_idx*Nq_WP*Np_WP + q_idx*Np_WP + p_idx;
}

inline size_t wp_overlap_index(size_t qp_idx_WP,
							   size_t pp_idx_WP,
							   size_t q_idx_WP,
							   size_t p_idx_WP,
							   size_t Nq_WP,
							   size_t Np_WP){
	return qp_idx_WP*(Np_WP*Nq_WP*Np_WP)
		 + pp_idx_WP*(Nq_WP*Np_WP)
		 + q_idx_WP*Np_WP
		 + p_idx_WP;
}

inline size_t phi_packet_index(size_t qp_idx_WP,
							   size_t pp_idx_WP,
							   size_t phi_idx,
							   size_t Np_WP,
							   size_t Nphi){
	return (qp_idx_WP*Np_WP + pp_idx_WP)*Nphi + phi_idx;
}

} // namespace
