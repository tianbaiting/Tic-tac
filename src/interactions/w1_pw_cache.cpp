#include "w1_pw_cache.h"

#include "w1_integrate.h"   // shared per-cell integration (bitwise-identical path)

#include "three_nucleon_force_model.h"
#include "constants.h"
#include "gauss_legendre.h"

#include <cmath>
#include <cstdio>
#include <omp.h>
#include <stdexcept>

#if TICTAC_USE_NEW_CACHE_LAYER
#include "io/cache_layer/cache_layer.h"
#include "cache_schema.h"
#include <cstring>  // memcpy
#endif

using tictac::interactions::integrate_w1_channel_blocks;

tictac::cache::W1Key make_w1_key(const three_nucleon_force_model& tnf,
                                 const run_params&                run_parameters,
                                 const pw_3N_statespace&          pw_states,
                                 int                              Np_WP,
                                 int                              Nq_WP,
                                 int                              Np_quad,
                                 int                              Nq_quad,
                                 int                              a_r,
                                 int                              a_c,
                                 const std::string&               p_grid_hash,
                                 const std::string&               q_grid_hash)
{
#if TICTAC_USE_NEW_CACHE_LAYER
    tictac::cache::W1Key k{};
    k.schema_version  = tictac::cache::W1_SCHEMA_VERSION;
    k.potential_model = run_parameters.potential_model;
    k.tnf_model       = tnf.name();
    k.Np_WP           = Np_WP;
    k.Nq_WP           = Nq_WP;
    k.J_2N_max        = run_parameters.J_2N_max;
    k.two_J_3N_max    = run_parameters.two_J_3N_max;
    k.two_J_3N        = pw_states.two_J_3N_array[a_r];
    k.P_3N            = pw_states.P_3N_array[a_r];
    k.a_r             = a_r;
    k.a_c             = a_c;
    k.c_D             = run_parameters.c_D;
    k.c_E             = run_parameters.c_E;
    k.Lambda_3NF      = run_parameters.Lambda_3NF;
    k.c_1             = tnf.lec_c1_gev();
    k.c_3             = tnf.lec_c3_gev();
    k.c_4             = tnf.lec_c4_gev();
    k.g_A             = tnf.axial_coupling_3nf();
    k.f_pi_MeV        = tnf.pion_decay_constant_mev_3nf();
    k.m_pi_MeV        = tnf.pion_mass_mev_3nf();
    k.Lambda_chi_MeV  = tnf.chiral_scale_mev_3nf();
    k.hbarc_MeV_fm    = tnf.hbarc_mev_fm_3nf();
    k.regulator_kind  = "gaussian";
    k.chebyshev_s        = run_parameters.chebyshev_s;
    k.chebyshev_t        = run_parameters.chebyshev_t;
    k.tensor_force       = run_parameters.tensor_force;
    k.isospin_breaking_1S0 = run_parameters.isospin_breaking_1S0;
    k.Np_per_WP_W1       = Np_quad;
    k.Nq_per_WP_W1       = Nq_quad;
    k.Nangle_3NF         = tnf.angular_order_3nf();
    k.p_grid_hash        = p_grid_hash;
    k.q_grid_hash        = q_grid_hash;
    return k;
#else
    // Legacy no-cache build: the key is unused, return a default.
    return tictac::cache::W1Key{};
#endif
}

void W1_PW_cache::build(const three_nucleon_force_model& tnf,
                        const double*               p_WP_array,
                        std::size_t                 Np_WP,
                        const double*               q_WP_array,
                        std::size_t                 Nq_WP,
                        const pw_3N_statespace&     pw_states,
                        const run_params&           run_parameters)
{
    m_Nalpha = (std::size_t)pw_states.Nalpha;
    m_Np     = Np_WP;
    m_Nq     = Nq_WP;

    // Enumerate (alpha_r, alpha_c) pairs allowed by 3N conservation: J_3N, T_3N, parity.
    m_block_index.assign(m_Nalpha * m_Nalpha, -1);
    m_blocks.clear();
    for (std::size_t a_r = 0; a_r < m_Nalpha; ++a_r) {
        for (std::size_t a_c = 0; a_c < m_Nalpha; ++a_c) {
            if (pw_states.two_J_3N_array[a_r] != pw_states.two_J_3N_array[a_c]) continue;
            if (pw_states.two_T_3N_array[a_r] != pw_states.two_T_3N_array[a_c]) continue;
            if (pw_states.P_3N_array[a_r]     != pw_states.P_3N_array[a_c])     continue;
            m_block_index[a_r * m_Nalpha + a_c] = (int)m_blocks.size();
            m_blocks.emplace_back((int)a_r, (int)a_c);
        }
    }

    const std::size_t per_block = m_Nq * m_Nq * m_Np * m_Np;
    m_data.assign(m_blocks.size() * per_block, 0.0);

    // [EN] Cache value semantics: for each (alpha_r, alpha_c, p_r, q_r, p_c, q_c) cell
    // tuple we store the WP bin matrix element of W^(1) in MeV (the diagnostic w1_scale
    // is applied at consumption time):
    //
    //   W^(1)_WP[i_p_r, i_q_r, i_p_c, i_q_c]
    //     = (1/hbarc^5)
    //       * (1 / sqrt(Δp_r·Δq_r·Δp_c·Δq_c))
    //       * ∫_bin_p_r dp_r ∫_bin_q_r dq_r ∫_bin_p_c dp_c ∫_bin_q_c dq_c
    //         p_r q_r p_c q_c · W^(1)(p_r,q_r,p_c,q_c)
    //
    // The 1/sqrt(Δ) per side comes from WP basis normalization (orthonormalized
    // packets); p·q per side is the radial measure left over from PW reduction;
    // 1/hbarc^5 converts the W^(1) output (fm^5) to MeV^{-5} so the bin matrix
    // element comes out in MeV, matching V_WP's MeV convention.
    //
    // For Np_per_WP_W1 = Nq_per_WP_W1 = 1 (legacy diagnostic) the integral collapses to the
    // single-point rule and the value reduces bit-for-bit to:
    //   p_r·q_r·p_c·q_c · sqrt(Δp_r·Δq_r·Δp_c·Δq_c) · W^(1)(midpoints) / hbarc^5
    // which is what the consumer in solve_faddeev.cpp used to compute inline,
    // so the legacy code path is preserved exactly.
    //
    // / [CN] 缓存值语义：每个 (α_r, α_c, p_r, q_r, p_c, q_c) cell 元组存放 W^(1) 在
    // WP 基下的 bin 矩阵元（MeV，单位与 V_WP 一致；w1_scale 诊断标度在消费侧再乘）。
    // Np_per_WP_W1 = Nq_per_WP_W1 = 1 时退化为现行 midpoint 公式，逐比特复现旧结果。

    // Quadrature orders are part of the cache key (make_w1_key) and select the
    // shared integration helper's per-bin Gauss-Legendre rule.  The actual
    // node/weight/normalisation pre-computation lives in integrate_w1_channel_blocks
    // so the monolithic dense path and the distributed per-block executor share
    // one source of truth (bitwise-identical W^(1) blocks).
    const int Np_quad = std::max(1, run_parameters.Np_per_WP_W1);
    const int Nq_quad = std::max(1, run_parameters.Nq_per_WP_W1);

    const std::size_t num_blocks = m_blocks.size();

    // 3NF audit B5 (2026-06-21): compute SHA-256 of the grid arrays once,
    // outside the parallel block loop. Different grids (different boundaries,
    // different mappings) produce different hashes, so the cache cannot
    // wrongly reuse blocks computed on an older grid.
    std::string p_grid_hash;
    std::string q_grid_hash;
#if TICTAC_USE_NEW_CACHE_LAYER
    {
        // Hash the binary representation: (N+1) doubles covering [0, p_max].
        const std::size_t p_bytes = (m_Np + 1) * sizeof(double);
        const std::size_t q_bytes = (m_Nq + 1) * sizeof(double);
        p_grid_hash = tictac::cache::hash_full_raw(
            reinterpret_cast<const unsigned char*>(p_WP_array), p_bytes);
        q_grid_hash = tictac::cache::hash_full_raw(
            reinterpret_cast<const unsigned char*>(q_WP_array), q_bytes);
    }
#endif

#if TICTAC_USE_NEW_CACHE_LAYER
    // 3NF audit B5 (2026-06-21): added c_1, c_3, c_4 + grid hashes to the key
    // so caches cannot be wrongly reused across different Hamiltonians or
    // momentum grids.
    auto build_w1_key_for_block = [&](int a_r, int a_c) {
        return make_w1_key(tnf, run_parameters, pw_states,
                           (int)m_Np, (int)m_Nq, Np_quad, Nq_quad,
                           a_r, a_c, p_grid_hash, q_grid_hash);
    };
#endif

    // Load complete slabs first, then batch all missing channel pairs at each
    // common Jacobi-momentum tuple.  The complete factorized projector can
    // thereby reuse its momentum-dependent orbital kernels across the whole
    // partial-wave block set.  The default model implementation still calls
    // W1_element independently, preserving behavior for other 3NF models.
    std::vector<std::size_t> missing_blocks;
    missing_blocks.reserve(num_blocks);
    for (std::size_t blk = 0; blk < num_blocks; ++blk) {
        double* slab = &m_data[blk * per_block];
#if TICTAC_USE_NEW_CACHE_LAYER
        const int a_r = m_blocks[blk].first;
        const int a_c = m_blocks[blk].second;
        const tictac::cache::W1Key k = build_w1_key_for_block(a_r, a_c);
        tictac::cache::W1Block out_block{};
        const auto res = tictac::cache::lookup_w1(k, &out_block);
        if (res.hit
            && static_cast<std::size_t>(out_block.Nq) == m_Nq
            && static_cast<std::size_t>(out_block.Np) == m_Np
            && out_block.data.size() == per_block)
        {
            std::memcpy(slab, out_block.data.data(), per_block * sizeof(double));
            continue;
        }
#endif
        missing_blocks.push_back(blk);
    }
	if (missing_blocks.empty()) return;

	// The complete factorized projector is exactly Hermitian at fixed transfer
	// order.  Evaluate only one orientation of a missing channel pair; if its
	// reverse slab was loaded from cache, evaluate neither.  Other models retain
	// the complete ordered-block path.
	std::vector<std::size_t> evaluation_blocks;
	evaluation_blocks.reserve(missing_blocks.size());
	if (tnf.W1_is_exactly_hermitian()) {
		std::vector<unsigned char> is_missing(num_blocks, 0);
		for (const std::size_t blk : missing_blocks) is_missing[blk] = 1;
		for (const std::size_t blk : missing_blocks) {
			const int a_r = m_blocks[blk].first;
			const int a_c = m_blocks[blk].second;
			const int reverse = m_block_index[
				static_cast<std::size_t>(a_c) * m_Nalpha
				+ static_cast<std::size_t>(a_r)];
			if (reverse < 0 || static_cast<std::size_t>(reverse) == blk) {
				evaluation_blocks.push_back(blk);
			} else if (is_missing[static_cast<std::size_t>(reverse)] != 0
			           && blk < static_cast<std::size_t>(reverse)) {
				evaluation_blocks.push_back(blk);
			}
		}
		if (evaluation_blocks.size() < missing_blocks.size()) {
			std::fprintf(stderr,
				"[3NF W1] exact Hermiticity: evaluating %zu of %zu missing "
				"channel blocks and transpose-filling %zu.\n",
				evaluation_blocks.size(), missing_blocks.size(),
				missing_blocks.size() - evaluation_blocks.size());
		}
	} else {
		evaluation_blocks = missing_blocks;
	}

	std::vector<std::pair<int, int>> evaluation_channels;
	evaluation_channels.reserve(evaluation_blocks.size());
	for (const std::size_t blk : evaluation_blocks) {
		evaluation_channels.push_back(m_blocks[blk]);
	}

	// [EN] Exact per-cell integration of every Hermitian-triangle block.
	// integrate_w1_channel_blocks is the shared source of truth: the monolithic
	// dense path here and the distributed per-block W1BlockExecutor call the same
	// function, so cached W^(1) values are bitwise-identical regardless of which
	// path produced them.  Each evaluation block occupies a contiguous per_block
	// region of m_data (block-major layout), so we hand the helper direct pointers
	// into m_data -- no copy, no value change.  The helper also performs the
	// momentum-independent prepare_W1_channel warm-up once per channel.
	// / [CN] 逐 cell 精确积分 Hermitian 三角块。与分布式 W1BlockExecutor 共用同一
	// 函数，保证两条路径逐比特一致；直接写入 m_data，无拷贝、无值变化。
	{
		std::vector<double*> out_ptrs;
		out_ptrs.reserve(evaluation_blocks.size());
		for (const std::size_t blk : evaluation_blocks) {
			out_ptrs.push_back(&m_data[blk * per_block]);
		}
		integrate_w1_channel_blocks(tnf, p_WP_array, m_Np, q_WP_array, m_Nq,
		                            pw_states, run_parameters,
		                            evaluation_channels, out_ptrs);
	}

	if (tnf.W1_is_exactly_hermitian()) {
		std::vector<unsigned char> was_evaluated(num_blocks, 0);
		for (const std::size_t blk : evaluation_blocks) was_evaluated[blk] = 1;
		for (const std::size_t blk : missing_blocks) {
			if (was_evaluated[blk] != 0) continue;
			const int a_r = m_blocks[blk].first;
			const int a_c = m_blocks[blk].second;
			const int reverse = m_block_index[
				static_cast<std::size_t>(a_c) * m_Nalpha
				+ static_cast<std::size_t>(a_r)];
			if (reverse < 0) {
				throw std::runtime_error(
					"exact-Hermitian W1 channel has no reverse block");
			}
			double* target = &m_data[blk * per_block];
			const double* source = &m_data[static_cast<std::size_t>(reverse) * per_block];
			for (std::size_t iqr = 0; iqr < m_Nq; ++iqr) {
				for (std::size_t iqc = 0; iqc < m_Nq; ++iqc) {
					for (std::size_t ipr = 0; ipr < m_Np; ++ipr) {
						for (std::size_t ipc = 0; ipc < m_Np; ++ipc) {
							const std::size_t target_cell =
								((iqr * m_Nq + iqc) * m_Np + ipr) * m_Np + ipc;
							const std::size_t source_cell =
								((iqc * m_Nq + iqr) * m_Np + ipc) * m_Np + ipr;
							target[target_cell] = source[source_cell];
						}
					}
				}
			}
		}
	}

#if TICTAC_USE_NEW_CACHE_LAYER
    for (const std::size_t blk : missing_blocks) {
        const int a_r = m_blocks[blk].first;
        const int a_c = m_blocks[blk].second;
        const double* slab = &m_data[blk * per_block];
        tictac::cache::W1Block in_block{};
        in_block.Nq  = static_cast<int>(m_Nq);
        in_block.Np  = static_cast<int>(m_Np);
        in_block.a_r = a_r;
        in_block.a_c = a_c;
        in_block.data.assign(slab, slab + per_block);
        tictac::cache::store_w1(build_w1_key_for_block(a_r, a_c), in_block);
    }
#endif
}
