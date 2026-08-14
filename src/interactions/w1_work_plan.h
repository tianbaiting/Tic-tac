#ifndef TICTAC_INTERACTIONS_W1_WORK_PLAN_H
#define TICTAC_INTERACTIONS_W1_WORK_PLAN_H

#include <cstddef>
#include <vector>

#include "type_defs.h"

namespace tictac::interactions {

// [EN] Slice one conserved J^pi channel out of the global basis: a non-owning
// pw_3N_statespace view over [chn_3N_idx_array[chn], chn_3N_idx_array[chn+1]).
// Identical to the slice the solver and w1_worker use; shared here so the work
// plan, manifest, and assembler agree on sector identity.
// / [CN] 从全局基切出一个守恒 J^pi 通道的非拥有视图。
pw_3N_statespace make_channel_view(const pw_3N_statespace& pw, int chn);

// [EN] Role of one W^(1) work unit within a conserved J^pi sector.
//
//   evaluate        -- the worker must directly integrate this block via
//                      W1BlockExecutor::compute_block and publish the shard.
//   transpose_fill  -- the block is the exact Hermitian transpose of its
//                      conjugate; workers MUST NOT integrate it.  It is
//                      produced at solve time by W1_PW_cache::build, which
//                      loads the (cached) conjugate and transpose-fills this
//                      orientation (the existing exact-Hermitian contract).
//
// The classification is deterministic and cache-state-independent: for each
// Hermitian pair {(a_r,a_c),(a_c,a_r)} the member with the SMALLER sector-local
// block id is `evaluate`, the other is `transpose_fill`.  Diagonal blocks
// (a_r==a_c) are always `evaluate`.  This is the absolute form of the rule
// already used inside W1_PW_cache::build, so pre-building all `evaluate`
// units makes the solver a pure cache-hit + cheap transpose-fill pass with
// zero expensive W1 evaluation.  Workers therefore never duplicate work and
// never evaluate both orientations of a pair.
//
// / [CN] 一个 W^(1) 工作单元在守恒 J^pi 扇区内的角色。evaluate 需由 worker 直接
// 积分；transpose_fill 是其共轭的精确 Hermitian 转置，worker 不积分，由求解器在
// load 已缓存的共轭后转置填充。规则确定且与缓存状态无关。
enum class W1UnitRole {
	evaluate = 0,
	transpose_fill = 1,
};

// [EN] One immutable work-unit description. `id` is the sector-local block id
// (position in the deterministic (a_r,a_c) enumeration), stable across runs
// and workers.  For transpose_fill units, `conj_a_r`/`conj_a_c` identify the
// evaluated conjugate whose transpose reproduces this block.
// / [CN] 一个不可变工作单元描述。id 为扇区内确定枚举位置，跨运行/worker 稳定。
struct W1WorkUnit {
	int          id = -1;
	int          a_r = -1;
	int          a_c = -1;
	W1UnitRole   role = W1UnitRole::evaluate;
	int          conj_a_r = -1;   // populated iff role == transpose_fill
	int          conj_a_c = -1;
};

// [EN] Deterministic work plan for one conserved J^pi sector (a channel view
// of pw_3N_statespace).  Pure function of (pw, hermitian): it enumerates the
// allowed (a_r,a_c) blocks in the same order as W1_PW_cache::build and
// classifies each as evaluate / transpose_fill.
//
// The plan is independent of the cache state, the grids, the LECs, and the
// model instance; those enter the per-block cache key (make_w1_key) and the
// integration (integrate_w1_channel_blocks).  This separation keeps the
// decomposition inspectable and testable without constructing a 3NF model.
//
// / [CN] 一个守恒 J^pi 扇区的确定性工作计划。是 (pw, hermitian) 的纯函数：按与
// W1_PW_cache::build 相同的顺序枚举允许的 (a_r,a_c) 块并分类。
class W1WorkPlan {
public:
	// `hermitian` should be tnf.W1_is_exactly_hermitian(). When false, every
	// allowed block is `evaluate` (no transpose reuse), matching the
	// W1_PW_cache::build fallback for non-Hermitian models.
	W1WorkPlan(const pw_3N_statespace& pw, bool hermitian);

	const std::vector<W1WorkUnit>& units() const { return m_units; }
	std::size_t num_blocks()        const { return m_units.size(); }
	std::size_t num_evaluate()      const { return m_num_eval; }
	std::size_t num_transpose_fill() const { return m_units.size() - m_num_eval; }
	int Nalpha() const { return m_Nalpha; }

	// (a_r,a_c) -> position in units(), or -1 if the pair is forbidden by
	// J_3N/T_3N/parity conservation.  Same indexing as W1_PW_cache::m_block_index.
	int block_index(int a_r, int a_c) const {
		return m_block_index[static_cast<std::size_t>(a_r) * m_Nalpha + a_c];
	}
	const W1WorkUnit* unit_for(int a_r, int a_c) const {
		const int idx = block_index(a_r, a_c);
		return idx < 0 ? nullptr : &m_units[idx];
	}

	// Evaluate-only units in deterministic (a_r,a_c) order.  Used by workers
	// for static partitioning: worker w owns units where
	// (global_eval_index + w) % worker_count == 0.
	std::vector<const W1WorkUnit*> evaluate_units() const;

private:
	std::vector<W1WorkUnit> m_units;
	std::vector<int>        m_block_index;   // Nalpha*Nalpha
	std::size_t             m_num_eval = 0;
	int                     m_Nalpha = 0;
};

} // namespace tictac::interactions

#endif // TICTAC_INTERACTIONS_W1_WORK_PLAN_H
