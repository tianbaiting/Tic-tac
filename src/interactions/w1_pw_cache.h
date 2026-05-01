#ifndef W1_PW_CACHE_H
#define W1_PW_CACHE_H

#include <cstddef>
#include <utility>
#include <vector>

#include "type_defs.h"

class three_nucleon_force_model;

// [EN] Dense W^(1) PW cache.
//
// The chiral N2LO 3NF W^(1) matrix element is independent of the outer column index
// inside the CPVC hot path: the same momentum/quantum-number tuple is recomputed
// O(Np_WP * num_Pade_iter) times per outer (alpha_c, q_c, p_c). Caching the
// midpoint-evaluated W1 once and looking it up in the inner loop eliminates the
// per-Padé-iter overhead.
//
// Storage layout (block-major):
//   blocks ↔ allowed (alpha_r, alpha_c) pairs (J_3N, T_3N, parity match)
//   per-block index = ((idx_q_r * Nq + idx_q_c) * Np + idx_p_r) * Np + idx_p_c
// Memory footprint = num_blocks * Nq^2 * Np^2 * 8 B; see total_bytes().
//
// / [CN] 稠密 W^(1) PW 缓存。在 CPVC 热路径中 W1 与外层列下标无关，但当前实现按
// 每个外层迭代/每次 Padé 重新评估同一动量与量子数组合。该缓存按 bin 中点把 W1 一
// 次性算好，热路径只做 O(1) 查表。
class W1_PW_cache {
public:
    void build(const three_nucleon_force_model& tnf,
               const double*               p_WP_array,
               std::size_t                 Np_WP,
               const double*               q_WP_array,
               std::size_t                 Nq_WP,
               const pw_3N_statespace&     pw_states,
               const run_params&           run_parameters);

    inline double get(int alpha_r, int alpha_c,
                      std::size_t idx_p_r, std::size_t idx_q_r,
                      std::size_t idx_p_c, std::size_t idx_q_c) const
    {
        const int blk = m_block_index[(std::size_t)alpha_r * m_Nalpha + (std::size_t)alpha_c];
        if (blk < 0) return 0.0;
        const std::size_t idx =
            ((((std::size_t)blk * m_Nq + idx_q_r) * m_Nq + idx_q_c) * m_Np + idx_p_r) * m_Np + idx_p_c;
        return (double)m_data[idx];
    }

    bool        valid()       const { return !m_data.empty(); }
    std::size_t Nalpha()      const { return m_Nalpha; }
    std::size_t Nq()          const { return m_Nq; }
    std::size_t Np()          const { return m_Np; }
    std::size_t num_blocks()  const { return m_blocks.size(); }
    std::size_t total_bytes() const { return m_data.size() * sizeof(float); }

private:
    std::vector<float>                   m_data;
    std::vector<int>                     m_block_index;  // (alpha_r, alpha_c) → block id, -1 if forbidden
    std::vector<std::pair<int, int>>     m_blocks;       // block id → (alpha_r, alpha_c)
    std::size_t                          m_Nalpha = 0;
    std::size_t                          m_Nq     = 0;
    std::size_t                          m_Np     = 0;
};

#endif // W1_PW_CACHE_H
