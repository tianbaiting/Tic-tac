#ifndef THREE_NUCLEON_FORCE_MODEL_H
#define THREE_NUCLEON_FORCE_MODEL_H

#include <string>

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
	// The full physical 3NF on the antisymmetrised 3N Hilbert space is
	//     W = W^(1) + W^(2) + W^(3) = (1 + P_{123} + P_{132})·W^(1)
	// which for antisymmetric bras and kets equals 3·W^(1) as a matrix element.
	//
	// HOWEVER, the Faddeev/AGS KERNEL acts in spectator-1 Faddeev space where the
	// ket is a single Faddeev component. The kernel is (Witała 2008 PRC 77 034004
	// eq. 3, Golak 2010 EPJA 43 241 §2):
	//     K_AGS = P·V + W^(1)·(1 + P)
	// i.e. W^(1) on the LEFT, (1+P) on the RIGHT, with P = P_{123} + P_{132}.
	//
	// One may NOT replace W^(1)(1+P) by (1+P)W^(1) in the kernel even though the
	// two coincide for fully antisymmetric states — W^(1) is defined ONLY in the
	// spectator-1 frame and is not the same object after a cyclic permutation.
	//
	// / [CN] 在 3N partial-wave 与 Jacobi 动量基下计算单个 W^(1) 矩阵元。p/q 为 Jacobi 动量
	// (fm^{-1})，而非 WP bin 索引。
	// 算符顺序约定（锁定）：AGS 核使用 W^(1)·(1+P)，W^(1) 在左、(1+P) 在右。完整的物理
	// 3NF 在反对称化 3N 希尔伯特空间中为 W = (1+P+P²)·W^(1)，但核内不能将两者混用。
	virtual double W1_element(int alpha_r, int alpha_c,
							  double p_r, double q_r,
							  double p_c, double q_c,
							  const pw_3N_statespace& pw_states) const { return 0.0; }
};

#endif // THREE_NUCLEON_FORCE_MODEL_H
