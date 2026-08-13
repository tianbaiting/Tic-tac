#ifndef THREE_NUCLEON_FORCE_MODEL_H
#define THREE_NUCLEON_FORCE_MODEL_H

#include <string>
#include <memory>
#include <utility>
#include <vector>

#include "type_defs.h"

// [EN] Abstract base for three-nucleon force (3NF) models. Mirrors the role of `potential_model` for 2NF, but is
// kept deliberately separate: 2NF and 3NF query signatures, quantum-number structure, and caching needs differ
// enough that a single virtual interface would degrade to a tag-union. The kernel assembler asks `enabled()`
// first; when a null-object (`three_nucleon_force_none`) is installed the hot-path 3NF branch is skipped and
// the solver reproduces the pre-3NF 2NF-only code path exactly. / [CN] 三体核力 (3NF) 模型的抽象基类。功能上
// 与 2NF 的 `potential_model` 对称，但刻意保持分离：2NF 与 3NF 的查询签名、量子数结构、缓存需求差异足够大，
// 塞进同一个虚接口会退化成联合 tag。内核组装器先问 `enabled()`；当装入 null 对象 (`three_nucleon_force_none`)
// 时，热路径上的 3NF 分支被跳过，求解器与 3NF 引入前的纯 2NF 路径完全一致。
class three_nucleon_force_model
{
public:
	three_nucleon_force_model() = default;
	virtual ~three_nucleon_force_model() = default;

	// Factory: dispatches on run_parameters.three_nucleon_force. Unknown strings raise an error; "none" (the
	// default) returns a three_nucleon_force_none null object so callers never need to null-check the pointer.
	// [EN] Modern factory returning an owning unique_ptr. Production code should
	// prefer create(); fetch() below is retained as a thin compatibility wrapper
	// that releases the pointer for legacy raw-pointer call sites. / [CN] 新工厂
	// 返回持有所有权的 unique_ptr，生产代码应优先使用；fetch() 作为兼容包装保留，
	// 向旧式裸指针调用点释放所有权。
	static std::unique_ptr<three_nucleon_force_model> create(run_params run_parameters);
	static three_nucleon_force_model* fetch(run_params run_parameters);

	// True only for concrete (real) 3NF implementations. The null object returns false so one test in the
	// kernel assembler skips the entire 3NF contribution path.
	virtual bool enabled() const = 0;

	// Display name used by logs and the run-parameter printout; must match the factory key.
	virtual std::string name() const = 0;

	// [EN] Optional hook for parameter walks (same role as potential_model::update_parameters). The base
	// implementation is a no-op so null and simple models don't need to override it. / [CN] 参数扫描的可选钩子
	// (与 potential_model::update_parameters 同样的角色)；基类给 no-op，null 与简单模型不必覆盖。
	virtual void update_parameters(const double* parameters) {}

	// [EN] Evaluate a single W^(1) matrix element in the 3N partial-wave & Jacobi-momentum basis:
	//   ⟨alpha_r, p_r, q_r | W^(1) | alpha_c, p_c, q_c⟩
	// where p/q are Jacobi momenta in fm^{-1} (NOT WP bin indices). alpha_r/alpha_c are indices into
	// pw_states arrays. W^(1) is the 3NF decomposition where particle 1 is the spectator: it is
	// symmetric under exchange of pair particles 2 and 3.
	//
	// OPERATOR ORDERING CONVENTION (locked, see docs/treatise/chapters/15_3nf_physics.tex
	// §operator-ordering and tests/cpp/test_faddeev_operator_order.cpp for the dense test):
	//
	// The full physical 3NF is the sum of the three cyclic spectator components,
	//     W = W^(1) + W^(2) + W^(3).
	// This decomposition must not be used to commute a spectator component
	// through P inside an unsaturated integral equation.
	//
	// Tic-tac solves the symmetrised elastic AGS equation of Deltuva,
	// Phys. Rev. C 80, 064002 (2009), Eq. (7a), after the exact channel-resolvent
	// reduction t G0 = V G and G0(1+t G0) = G.  Its kernel is
	//     K_AGS = P·V + (1 + P)·W^(1)
	// i.e. (1+P) is on the LEFT of W^(1), with P = P_{123} + P_{132}.
	//
	// W^(1)(1+P) belongs to the distinct Faddeev breakup-component equation.  It
	// may NOT replace (1+P)W^(1) here when [P,W^(1)] != 0.
	//
	// / [CN] 在 3N partial-wave 与 Jacobi 动量基下计算单个 W^(1) 矩阵元。p/q 为 Jacobi 动量
	// (fm^{-1})，而非 WP bin 索引。
	// 算符顺序约定（锁定）：弹性 AGS 核使用 (1+P)·W^(1)，(1+P) 在左；
	// W^(1)·(1+P) 属于另一条 Faddeev 分量方程，二者不可在核内混用。
	virtual double W1_element(int alpha_r, int alpha_c,
								  double p_r, double q_r,
								  double p_c, double q_c,
								  const pw_3N_statespace& pw_states) const { return 0.0; }

	// Evaluate several channel pairs at one common Jacobi-momentum tuple.  The
	// default preserves the scalar interface exactly; implementations with a
	// factorized angular kernel may override it to share momentum-dependent
	// orbital work across channel pairs.
	virtual void W1_elements_for_channels(
		const std::vector<std::pair<int, int>>& channels,
		double p_r, double q_r, double p_c, double q_c,
		const pw_3N_statespace& pw_states,
		std::vector<double>& values) const
	{
		values.resize(channels.size());
		for (std::size_t index = 0; index < channels.size(); ++index) {
			values[index] = W1_element(
				channels[index].first, channels[index].second,
				p_r, q_r, p_c, q_c, pw_states);
		}
	}

	// Models may opt in only when their finite-order implementation obeys
	// W(ar,ac; pr,qr,pc,qc) = W(ac,ar; pc,qc,pr,qr) exactly.  The W1 packet
	// builder can then integrate one orientation and transpose the reverse
	// channel block.  The conservative default protects diagnostic projectors
	// whose finite quadrature is not Hermitian at a fixed order.
	virtual bool W1_is_exactly_hermitian() const { return false; }

	// Optional metadata-only warm-up before a parallel W1 packet build.  Models
	// with shared momentum-independent recoupling tables may populate them here
	// without evaluating and discarding a full momentum-space matrix element.
	virtual void prepare_W1_channel(
		int alpha_r, int alpha_c,
		const pw_3N_statespace& pw_states) const {}

	// [EN] LEC accessors used by W1_PW_cache to build a complete cache key
	// (3NF audit B5: previously the cache key omitted c_1, c_3, c_4, allowing
	// two runs with different 2NF (e.g. N2LOopt vs Idaho_N3LO) but identical
	// c_D, c_E, Λ to wrongly reuse the cached W^(1) matrix elements).
	// Default implementations return 0 (null/stub models). Real chiral models
	// override these to expose their internal LECs in **GeV⁻¹** units
	// (matching the convention used in run_params and constants.h).
	// / [CN] LEC 访问器，供 W1_PW_cache 构造完整 cache key（3NF 审计 B5）。
	virtual double lec_c1_gev() const { return 0.0; }
	virtual double lec_c3_gev() const { return 0.0; }
	virtual double lec_c4_gev() const { return 0.0; }
	// Angular quadrature order that changes W1 matrix elements.  Zero denotes a
	// model with no runtime angular-order control.
	virtual int angular_order_3nf() const { return 0; }
	// Constants that define a chiral W1 kernel.  They are explicit cache-key
	// inputs so recompiling with a different convention cannot reuse old blocks.
	virtual double axial_coupling_3nf() const { return 0.0; }
	virtual double pion_decay_constant_mev_3nf() const { return 0.0; }
	virtual double pion_mass_mev_3nf() const { return 0.0; }
	virtual double chiral_scale_mev_3nf() const { return 0.0; }
	virtual double hbarc_mev_fm_3nf() const { return 0.0; }
};

#endif // THREE_NUCLEON_FORCE_MODEL_H
