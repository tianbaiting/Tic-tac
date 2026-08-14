#ifndef TICTAC_INTERACTIONS_W1_INTEGRATE_H
#define TICTAC_INTERACTIONS_W1_INTEGRATE_H

#include <cstddef>
#include <utility>
#include <vector>

#include "type_defs.h"

class three_nucleon_force_model;

namespace tictac::interactions {

// [EN] Shared exact W^(1) per-cell Gauss-Legendre integration.
//
// This is the single source of truth for the four-dimensional radial-bin
// integration of one W^(1) channel block.  Both the dense W1_PW_cache::build
// (monolithic, batched across the Hermitian triangle of a sector) and the
// distributed W1BlockExecutor (one block at a time, for resumable/distributed
// construction) call integrate_w1_channel_blocks, so the two paths produce
// bitwise-identical blocks by construction -- the per-cell accumulation order,
// the per-bin quadrature nodes/weights, the WP-bin normalisation, and the
// 1/hbarc^5 conversion are all shared, not re-transcribed.
//
// The value semantics are documented in W1_PW_cache::build; nothing here
// changes them.  No low-rank, no compression, no approximation: this is the
// exact factorized N2LO 3NF W^(1) bin matrix element.
//
// / [CN] 精确 W^(1) 单元高斯积分的共享实现。稠密 W1_PW_cache::build（整块批量）
// 与分布式 W1BlockExecutor（逐块）均调用 integrate_w1_channel_blocks，因此两条
// 路径按构造逐比特一致。不含任何低秩/压缩/近似。

// Build per-bin Gauss-Legendre nodes/weights along one axis (MeV convention).
// nodes_out[i*Npts + k] is the k-th Gauss point in bin i (MeV);
// weights_out[i*Npts + k] is the corresponding weight in MeV (already mapped
// from [-1,1] to [bin_lower, bin_upper]).
//
// For Npts=1 the single-point rule is forced to (midpoint, bin_width); this
// keeps the value bit-for-bit identical to the legacy midpoint formula when
// Np_per_WP_W1 = Nq_per_WP_W1 = 1.
void build_per_bin_quadrature(const double* boundary_array,
                              std::size_t   N_bins,
                              int           Npts,
                              std::vector<double>& nodes_out,
                              std::vector<double>& weights_out);

// [EN] Compute one dense W^(1) block per channel pair in `channels`, writing
// the Nq_WP*Nq_WP*Np_WP*Np_WP double payload of channel k into the buffer
// pointed to by out_block_ptrs[k].  Pure function of (tnf, grids, quadrature
// orders, channels).  OpenMP-parallel over the per_block cells.
//
// `channels` and `out_block_ptrs` must have the same size; each entry of
// out_block_ptrs must point to a contiguous buffer of at least
// Nq_WP*Nq_WP*Np_WP*Np_WP doubles.  The caller owns all buffers.
//
// The model's prepare_W1_channel() is called once per channel before the cell
// loop (matching the monolithic path's momentum-independent warm-up).
//
// / [CN] 为 `channels` 中每个通道对计算一个稠密 W^(1) 块，写入 out_block_ptrs[k]
// 指向的缓冲区。纯函数；OpenMP 并行按 cell。
void integrate_w1_channel_blocks(
	const three_nucleon_force_model& tnf,
	const double* p_WP_array, std::size_t Np_WP,
	const double* q_WP_array, std::size_t Nq_WP,
	const pw_3N_statespace& pw_states,
	const run_params& rp,
	const std::vector<std::pair<int, int>>& channels,
	const std::vector<double*>& out_block_ptrs);

// [EN] Thin executor that computes ONE exact W^(1) block. This is the
// distributed-construction work-horse: a worker calls compute_block() for each
// evaluate-unit it owns, then publishes the result via the hash-keyed cache
// (tictac::cache::store_w1), which atomically writes the shard.  The result is
// bitwise identical to the corresponding block produced by
// W1_PW_cache::build, because both call integrate_w1_channel_blocks.
//
// `out` must point to at least Nq_WP*Nq_WP*Np_WP*Np_WP doubles.
// / [CN] 计算单个精确 W^(1) 块的执行器；分布式 worker 对其分到的每个 evaluate
// 单元调用一次，再交给哈希键缓存原子落盘。结果与 W1_PW_cache::build 逐比特一致。
struct W1BlockExecutor {
	static void compute_block(const three_nucleon_force_model& tnf,
	                          const double* p_WP_array, std::size_t Np_WP,
	                          const double* q_WP_array, std::size_t Nq_WP,
	                          const pw_3N_statespace& pw_states,
	                          const run_params& rp,
	                          int a_r, int a_c,
	                          double* out);
};

} // namespace tictac::interactions

#endif // TICTAC_INTERACTIONS_W1_INTEGRATE_H
