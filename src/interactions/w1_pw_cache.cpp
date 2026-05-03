#include "w1_pw_cache.h"

#include "three_nucleon_force_model.h"
#include "constants.h"

#include <omp.h>

#if TICTAC_USE_NEW_CACHE_LAYER
#include "io/cache_layer/cache_layer.h"
#include "cache_schema.h"
#include <cstring>  // memcpy
#endif

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

    // Bin midpoints in fm^-1 (W1_element API expects Jacobi momenta in fm^-1).
    const double inv_hbarc = 1.0 / hbarc;
    std::vector<double> p_mid(m_Np), q_mid(m_Nq);
    for (std::size_t i = 0; i < m_Np; ++i)
        p_mid[i] = 0.5 * (p_WP_array[i] + p_WP_array[i + 1]) * inv_hbarc;
    for (std::size_t i = 0; i < m_Nq; ++i)
        q_mid[i] = 0.5 * (q_WP_array[i] + q_WP_array[i + 1]) * inv_hbarc;

    const std::size_t num_blocks = m_blocks.size();

#if TICTAC_USE_NEW_CACHE_LAYER
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
        // Today every supported chiral 3NF uses a Gaussian regulator; if a
        // different regulator family is added later, expose it through
        // three_nucleon_force_model and source it here instead of hardcoding.
        k.regulator_kind  = "gaussian";
        k.chebyshev_s        = run_parameters.chebyshev_s;
        k.chebyshev_t        = run_parameters.chebyshev_t;
        k.tensor_force       = run_parameters.tensor_force;
        k.isospin_breaking_1S0 = run_parameters.isospin_breaking_1S0;
        return k;
    };
#endif

    #pragma omp parallel for schedule(dynamic)
    for (std::size_t blk = 0; blk < num_blocks; ++blk) {
        const int a_r = m_blocks[blk].first;
        const int a_c = m_blocks[blk].second;
        float*    slab = &m_data[blk * per_block];

#if TICTAC_USE_NEW_CACHE_LAYER
        tictac::cache::W1Key k = build_w1_key_for_block(a_r, a_c);
        tictac::cache::W1Block out_block{};
        auto res = tictac::cache::lookup_w1(k, &out_block);
        if (res.hit
            && (size_t)out_block.Nq == m_Nq
            && (size_t)out_block.Np == m_Np
            && out_block.data.size() == per_block)
        {
            std::memcpy(slab, out_block.data.data(), per_block * sizeof(float));
            continue;
        }
#endif

        // Cache miss: do the original computation.
        for (std::size_t iqr = 0; iqr < m_Nq; ++iqr) {
            for (std::size_t iqc = 0; iqc < m_Nq; ++iqc) {
                for (std::size_t ipr = 0; ipr < m_Np; ++ipr) {
                    for (std::size_t ipc = 0; ipc < m_Np; ++ipc) {
                        slab[((iqr * m_Nq + iqc) * m_Np + ipr) * m_Np + ipc]
                            = (float)tnf.W1_element(a_r, a_c,
                                                    p_mid[ipr], q_mid[iqr],
                                                    p_mid[ipc], q_mid[iqc],
                                                    pw_states);
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
