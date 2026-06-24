#include "w1_pw_cache.h"

#include "three_nucleon_force_model.h"
#include "constants.h"
#include "gauss_legendre.h"

#include <cmath>
#include <omp.h>

#if TICTAC_USE_NEW_CACHE_LAYER
#include "io/cache_layer/cache_layer.h"
#include "cache_schema.h"
#include <cstring>  // memcpy
#endif

namespace {

// Build per-bin Gauss-Legendre nodes/weights along one axis.
// nodes_out[i*Npts + k] is the k-th Gauss point in bin i (in MeV).
// w_axis_out[i*Npts + k] is the corresponding weight in MeV (already mapped
// from [-1,1] to [bin_lower, bin_upper]).
//
// For Npts=1 the single-point rule is forced to (midpoint, bin_width); this
// keeps the cache value bit-for-bit identical to the legacy midpoint formula
// when Np_per_WP_W1 = Nq_per_WP_W1 = 1.
void build_per_bin_quadrature(const double* boundary_array,
                              std::size_t   N_bins,
                              int           Npts,
                              std::vector<double>& nodes_out,
                              std::vector<double>& weights_out)
{
    nodes_out.assign(N_bins * (std::size_t)Npts, 0.0);
    weights_out.assign(N_bins * (std::size_t)Npts, 0.0);

    if (Npts == 1) {
        for (std::size_t i = 0; i < N_bins; ++i) {
            double lo = boundary_array[i];
            double hi = boundary_array[i + 1];
            nodes_out[i]   = 0.5 * (lo + hi);
            weights_out[i] = hi - lo;          // single-point rule weight = bin width
        }
        return;
    }

    std::vector<double> xref(Npts), wref(Npts);
    gauss(xref.data(), wref.data(), Npts);     // nodes/weights on [-1, 1]

    for (std::size_t i = 0; i < N_bins; ++i) {
        double lo = boundary_array[i];
        double hi = boundary_array[i + 1];
        double a  = 0.5 * (hi - lo);
        double b  = 0.5 * (hi + lo);
        double*  n  = &nodes_out[i * (std::size_t)Npts];
        double*  w  = &weights_out[i * (std::size_t)Npts];
        for (int k = 0; k < Npts; ++k) {
            n[k] = a * xref[k] + b;
            w[k] = a * wref[k];
        }
    }
}

}  // namespace

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
    // For Np_per_WP_W1 = Nq_per_WP_W1 = 1 (default) the integral collapses to the
    // single-point rule and the value reduces bit-for-bit to:
    //   p_r·q_r·p_c·q_c · sqrt(Δp_r·Δq_r·Δp_c·Δq_c) · W^(1)(midpoints) / hbarc^5
    // which is what the consumer in solve_faddeev.cpp used to compute inline,
    // so the legacy code path is preserved exactly.
    //
    // / [CN] 缓存值语义：每个 (α_r, α_c, p_r, q_r, p_c, q_c) cell 元组存放 W^(1) 在
    // WP 基下的 bin 矩阵元（MeV，单位与 V_WP 一致；w1_scale 诊断标度在消费侧再乘）。
    // Np_per_WP_W1 = Nq_per_WP_W1 = 1 时退化为现行 midpoint 公式，逐比特复现旧结果。

    const double inv_hbarc  = 1.0 / hbarc;
    const double inv_hbarc5 = std::pow(inv_hbarc, 5);

    const int Np_quad = std::max(1, run_parameters.Np_per_WP_W1);
    const int Nq_quad = std::max(1, run_parameters.Nq_per_WP_W1);

    // Per-bin Gauss nodes/weights in MeV (consumer convention).
    std::vector<double> p_nodes_MeV, p_w_MeV, q_nodes_MeV, q_w_MeV;
    build_per_bin_quadrature(p_WP_array, m_Np, Np_quad, p_nodes_MeV, p_w_MeV);
    build_per_bin_quadrature(q_WP_array, m_Nq, Nq_quad, q_nodes_MeV, q_w_MeV);

    // Pre-compute WP normalization factor 1/sqrt(Δ) per bin.
    std::vector<double> p_inv_sqrtD(m_Np), q_inv_sqrtD(m_Nq);
    for (std::size_t i = 0; i < m_Np; ++i)
        p_inv_sqrtD[i] = 1.0 / std::sqrt(p_WP_array[i + 1] - p_WP_array[i]);
    for (std::size_t i = 0; i < m_Nq; ++i)
        q_inv_sqrtD[i] = 1.0 / std::sqrt(q_WP_array[i + 1] - q_WP_array[i]);

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
        tictac::cache::W1Key k{};
        k.schema_version  = tictac::cache::W1_SCHEMA_VERSION;
        k.potential_model = run_parameters.potential_model;
        k.tnf_model       = tnf.name();
        k.Np_WP           = (int)m_Np;
        k.Nq_WP           = (int)m_Nq;
        k.J_2N_max        = run_parameters.J_2N_max;
        k.two_J_3N_max    = run_parameters.two_J_3N_max;
        k.two_J_3N        = pw_states.two_J_3N_array[a_r];
        k.P_3N            = pw_states.P_3N_array[a_r];
        k.a_r             = a_r;
        k.a_c             = a_c;
        k.c_D             = run_parameters.c_D;
        k.c_E             = run_parameters.c_E;
        k.Lambda_3NF      = run_parameters.Lambda_3NF;
        // NEW (B5): c_1, c_3, c_4 from the model itself (default 0 for null/stub).
        k.c_1             = tnf.lec_c1_gev();
        k.c_3             = tnf.lec_c3_gev();
        k.c_4             = tnf.lec_c4_gev();
        // Today every supported chiral 3NF uses a Gaussian regulator; if a
        // different regulator family is added later, expose it through
        // three_nucleon_force_model and source it here instead of hardcoding.
        k.regulator_kind  = "gaussian";
        k.chebyshev_s        = run_parameters.chebyshev_s;
        k.chebyshev_t        = run_parameters.chebyshev_t;
        k.tensor_force       = run_parameters.tensor_force;
        k.isospin_breaking_1S0 = run_parameters.isospin_breaking_1S0;
        k.Np_per_WP_W1       = Np_quad;
        k.Nq_per_WP_W1       = Nq_quad;
        // NEW (B5): SHA-256 of the binary representation of the grid arrays.
        // Computed once outside the block loop and captured by reference.
        k.p_grid_hash        = p_grid_hash;
        k.q_grid_hash        = q_grid_hash;
        return k;
    };
#endif

    #pragma omp parallel for schedule(dynamic)
    for (std::size_t blk = 0; blk < num_blocks; ++blk) {
        const int a_r = m_blocks[blk].first;
        const int a_c = m_blocks[blk].second;
        // 3NF audit B6: storage is now double (was float).
        double*   slab = &m_data[blk * per_block];

#if TICTAC_USE_NEW_CACHE_LAYER
        tictac::cache::W1Key k = build_w1_key_for_block(a_r, a_c);
        tictac::cache::W1Block out_block{};
        auto res = tictac::cache::lookup_w1(k, &out_block);
        if (res.hit
            && (size_t)out_block.Nq == m_Nq
            && (size_t)out_block.Np == m_Np
            && out_block.data.size() == per_block)
        {
            std::memcpy(slab, out_block.data.data(), per_block * sizeof(double));
            continue;
        }
#endif

        // Cache miss: integrate W^(1) over each (p_r, q_r, p_c, q_c) bin tuple.
        for (std::size_t iqr = 0; iqr < m_Nq; ++iqr) {
            const double inv_sqrt_dq_r = q_inv_sqrtD[iqr];
            for (std::size_t iqc = 0; iqc < m_Nq; ++iqc) {
                const double inv_sqrt_dq_c = q_inv_sqrtD[iqc];
                for (std::size_t ipr = 0; ipr < m_Np; ++ipr) {
                    const double inv_sqrt_dp_r = p_inv_sqrtD[ipr];
                    for (std::size_t ipc = 0; ipc < m_Np; ++ipc) {
                        const double inv_sqrt_dp_c = p_inv_sqrtD[ipc];

                        double accum = 0.0;
                        for (int kqr = 0; kqr < Nq_quad; ++kqr) {
                            const double q_r_MeV = q_nodes_MeV[iqr * Nq_quad + kqr];
                            const double w_qr    = q_w_MeV   [iqr * Nq_quad + kqr];
                            const double q_r_fm  = q_r_MeV * inv_hbarc;
                            for (int kqc = 0; kqc < Nq_quad; ++kqc) {
                                const double q_c_MeV = q_nodes_MeV[iqc * Nq_quad + kqc];
                                const double w_qc    = q_w_MeV   [iqc * Nq_quad + kqc];
                                const double q_c_fm  = q_c_MeV * inv_hbarc;
                                for (int kpr = 0; kpr < Np_quad; ++kpr) {
                                    const double p_r_MeV = p_nodes_MeV[ipr * Np_quad + kpr];
                                    const double w_pr    = p_w_MeV   [ipr * Np_quad + kpr];
                                    const double p_r_fm  = p_r_MeV * inv_hbarc;
                                    for (int kpc = 0; kpc < Np_quad; ++kpc) {
                                        const double p_c_MeV = p_nodes_MeV[ipc * Np_quad + kpc];
                                        const double w_pc    = p_w_MeV   [ipc * Np_quad + kpc];
                                        const double p_c_fm  = p_c_MeV * inv_hbarc;

                                        const double w1 = tnf.W1_element(a_r, a_c,
                                                                          p_r_fm, q_r_fm,
                                                                          p_c_fm, q_c_fm,
                                                                          pw_states);
                                        accum += (p_r_MeV * q_r_MeV * p_c_MeV * q_c_MeV)
                                               * (w_pr * w_qr * w_pc * w_qc)
                                               * w1;
                                    }
                                }
                            }
                        }

                        const double bin_norm = inv_sqrt_dp_r * inv_sqrt_dq_r * inv_sqrt_dp_c * inv_sqrt_dq_c;
                        slab[((iqr * m_Nq + iqc) * m_Np + ipr) * m_Np + ipc]
                            = accum * bin_norm * inv_hbarc5;
                    }
                }
            }
        }

#if TICTAC_USE_NEW_CACHE_LAYER
        tictac::cache::W1Block in_block{};
        in_block.Nq  = (int)m_Nq;
        in_block.Np  = (int)m_Np;
        in_block.a_r = a_r;
        in_block.a_c = a_c;
        in_block.data.assign(slab, slab + per_block);
        tictac::cache::store_w1(k, in_block);
#endif
    }
}
