#ifndef TNF_KERNEL_CONTEXT_H
#define TNF_KERNEL_CONTEXT_H

#include "type_defs.h"

// [EN] Narrow header carrying only the tnf_kernel_context bundle shared between
// solve_faddeev.h and cpvc_kernel.h, so cpvc_kernel.h no longer needs to pull in
// the whole solve_faddeev.h (and its transitive HDF5 / state-space baggage).
// Both the production Faddeev driver and the isolated kernel-algebra builders
// depend on this single small header. / [CN] 窄头文件，仅承载 solve_faddeev.h 与
// cpvc_kernel.h 共享的 tnf_kernel_context，使 cpvc_kernel.h 不再需要包含整个
// solve_faddeev.h（及其 HDF5 / 态空间依赖）。

class three_nucleon_force_model;
class W1_PW_cache;

// [EN] Bundles all data the column-computation hot path needs for the 3NF contribution so that
// existing function signatures stay manageable. When tnf->enabled()==false (null object) the entire
// 3NF branch is skipped via a single test. / [CN] 把列计算热路径所需的全部 3NF 数据打包到一个结构体中，
// 避免已有函数签名膨胀。当 tnf->enabled()==false（null 对象）时，整个 3NF 分支通过一次判断跳过。
struct tnf_kernel_context {
	const three_nucleon_force_model* tnf;
	const pw_3N_statespace*          pw_states;
	const double*                    p_WP_array;   // WP boundaries, size Np_WP+1
	const double*                    q_WP_array;   // WP boundaries, size Nq_WP+1
	double**                         CT_RM_array;  // C^T row-major = C column-major, [Nalpha*Nalpha]
	double                           w1_scale;     // Overall scale factor applied to W^(1) output (diagnostic knob)
	const W1_PW_cache*               w1_cache;     // optional: radial-cell-integrated W1 lookup table (nullptr -> midpoint diagnostic fallback)
};

#endif // TNF_KERNEL_CONTEXT_H
