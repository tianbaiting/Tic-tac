
#ifndef TYPE_DEFS_H
#define TYPE_DEFS_H

#include <complex>
#include <string>

// [EN] Keep this header as the canonical solver-wide type catalog; `src/type_defs.h` is only a compatibility
// redirect so the refactored tree and the legacy include layout cannot silently drift apart. / [CN] 保持此头文件作为全求解器共享类型的唯一真相源；`src/type_defs.h` 只是兼容性转发层，避免重构目录与旧目录的类型定义悄悄漂移。

typedef unsigned int uint;
//~ typedef float floatType;
typedef double floatType;
typedef std::complex<floatType> cfloatType;
typedef std::complex<double> cdouble;
typedef unsigned long long int ull_int;

// [EN] The three-body basis is stored channel-by-channel in contiguous alpha ranges, so each conserved J^pi block
// can be solved independently after the global basis is generated. This mirrors the block structure used in the
// WPCD Faddeev formulation discussed in the project docs and in Miller et al., Phys. Rev. C 106, 024001 (2022).
// [CN] 三体基以连续 alpha 区间按守恒的 J^pi 通道分块存储，因此在生成全局基之后，每个块都可以独立求解；这对应项目文档和 Miller 等人在 Phys. Rev. C 106, 024001 (2022) 中讨论的 WPCD Faddeev 分块结构。
typedef struct pw_3N_statespace{
	int  Nalpha;				// Number of partial waves, set in dynamical state-space construction (by construct_symmetric_pw_states)
	int  J_2N_max;				// Maximum pair-system total angular momentum
	int* L_2N_array;			// pair-nucleon state angular momentum
	int* S_2N_array;			// pair-nucleon state total spin
	int* J_2N_array;			// pair-nucleon state total angular momentum
	int* T_2N_array;			// pair-nucleon state total isospin
	int* L_1N_array;			// orbital-nucleon state angular momentum
	int* two_J_1N_array;		// orbital-nucleon state total angular momentum x2
	int* two_J_3N_array;		// three-nucleon state total angular momentum x2
	int* two_T_3N_array;		// three-nucleon state total isospin x2
	int* P_3N_array;			// three-nucleon state parity
	int* chn_3N_idx_array;		// Indices of each 3N JP-channel
	int  N_chn_3N;				// Number of distinct 3N JP-channels
} pw_3N_statespace;

// [EN] fWP stores the free Jacobi-momentum packet mesh that replaces continuum integrals by finite bins plus
// quadrature weights. All later operators are projected onto this discretized p-q lattice first. / [CN] fWP 保存自由 Jacobi 动量波包网格，用有限 bin 和求积权重取代连续积分；后续所有算符都先投影到这套离散的 p-q 格点上。
typedef struct fwp_statespace{
	int 	Np_WP;				// Number of p-momentum WPs
	int 	Nq_WP;				// Number of q-momentum WPs
	int 	Np_per_WP;			// Number of p-momentum quadrature points in each p-momentum WP
	int 	Nq_per_WP;			// Number of q-momentum quadrature points in each q-momentum WP
	double* p_WP_array;			// Boundaries of p-momentum WPs
	double* q_WP_array;			// Boundaries of q-momentum WPs
	double* p_array;			// p-momentum quadrature points/average momenta, for all bins
	double* wp_array;			// p-momentum quadrature weights, for all bins
	double* fp_array;			// p-momentum weight-function at p_array-points, for all bins
	double* q_array;			// q-momentum quadrature points/average momenta, for all bins
	double* wq_array;			// q-momentum quadrature weights, for all bins
	double* fq_array;			// q-momentum weight-function at q_array-points, for all bins
	double* norm_p_array;		// p-momentum normaliszation, for all bins
	double* norm_q_array;		// q-momentum normaliszation, for all bins
} fwp_statespace;

// [EN] SWP stores the interacting packet basis obtained by diagonalizing the pair Hamiltonian inside each 2N
// channel. In this basis the channel resolvent is diagonal apart from the bound/continuum branch choice, which is
// the key WPCD simplification exploited by the solver. / [CN] SWP 保存每个两体通道内对角化相互作用哈密顿量后得到的散射波包基；在这套基中，通道分辨算符除束缚/连续分支选择外近似对角，这正是求解器利用的 WPCD 核心简化。
typedef struct swp_statespace{
	int 	Np_WP;				// Number of p-momentum WPs
	int 	Nq_WP;				// Number of q-momentum WPs
	int 	num_2N_unco_states;	// Number of uncoupled NN states
	int 	num_2N_coup_states;	// Number of coupled NN states
	double  E_bound;			// Deuteron bound state energy (with negative sign)
	double* e_SWP_unco_array;	// Boundaries of uncoupled scattering WPs
	double* e_SWP_coup_array;	// Boundaries of coupled scattering WPs
	double* C_SWP_unco_array;	// Basis transformation matrices for uncoupled WPs
	double* C_SWP_coup_array;	// Basis transformation matrices for coupled WPs
	double* q_WP_array;			// Boundaries of q-momentum WPs
} swp_statespace;

// [EN] This is the global energy-selection record: user-requested laboratory energies are first mapped to discrete
// q-shell bins here, then each J^pi channel receives a local view of the relevant entries. This matches the Sean B.
// Miller workflow where observables are extracted only from the finite set of on-shell packet amplitudes.
// / [CN] 这是全局的能量选择记录：用户请求的实验室能量先在这里映射到离散 q-shell bin，然后每个 J^pi 通道再拿到
// 相关条目的局部视图。这对应 Sean B. Miller 工作流里“只从有限个 on-shell 波包振幅提取可观测量”的做法。
// [EN] These arrays translate the large packet lattice back to the physical elastic asymptotic states. Only the
// deuteron rows/columns at the on-shell q bins are finally reported as nd/pd observables. / [CN] 这些数组负责把巨大的波包格点重新映射回物理弹性渐近态；最终真正输出为 nd/pd 可观测量的，只是 on-shell q bin 上的氘核行/列。
typedef struct solution_configuration{
	size_t  num_T_lab;				// Number of on-shell bins/energies to calculate
	double* T_lab_array;			// On-shell lab  energies  (T_lab)
	double* q_com_array;			// On-shell c.m. q-momenta (q_com)
	double* E_com_array;			// On-shell c.m. energies  (E_com)
	int*    q_com_idx_array;		// Index-array for on-shell q-bin (used for elastic scattering)
	int**   deuteron_idx_arrays;	// Index-arrays deuteron-channels in all 3N-channels
	int*    deuteron_num_array;		// Contains number of deuteron-channels in all 3N-channels
} solution_configuration;

// [EN] Channel-local on-shell bookkeeping lets the solver work blockwise: a single J^pi channel receives its own
// elastic rows, breakup rows, and q-bin map without carrying unrelated channels through the hot path. / [CN] 通道局部的 on-shell 索引使求解器可以逐块运行：单个 J^pi 通道只携带自身的弹性行、裂变行和 q-bin 映射，不把无关通道带入热点路径。
typedef struct channel_os_indexing{	// This should be renamed solution_configuration_subchannel
	size_t  num_T_lab;				// Number of elastic on-shell bins/energies to calculate
	int*    q_com_idx_array;		// Index-array for elastic on-shell q-bins
	int*    deuteron_idx_array;		// Index-array for deuteron-channels in given 3N-channel
	int     num_deuteron_states;	// Number of deuteron-channels in given 3N-channel
	int*    alphapq_idx_array;		// Flattened triples (alpha, q, p) for breakup on-shell packet states
	int*    q_com_BU_idx_array;		// Prefix offsets that group breakup triples by elastic q-shell bin
	size_t  num_BU_chns;			// Number of breakup channels (combinations of p,q, and alpha)
} channel_os_indexing;

// [EN] Runtime controls mix physics truncations, discretization knobs, and IO locations because the current solver
// still exposes a legacy key=value interface. Keeping them in one record avoids long argument chains while the
// build is migrated from `CPP/` to `src/`. / [CN] 运行参数同时包含物理截断、离散化旋钮和 IO 路径，因为当前求解器仍然对外暴露旧式 key=value 接口；在从 `CPP/` 迁移到 `src/` 的过程中，把它们放在同一记录里能避免冗长的参数链。
typedef struct run_params{
	int         two_J_3N_max;
	int         Np_WP;
	int         Nq_WP;
	int         J_2N_max;
	int         Nphi;
	int         Nx;
	int         Np_per_WP;
	int         Nq_per_WP;
	// [EN] Quadrature points per WP cell when building the W^(1) bin matrix elements
	// (3NF cache). Default 1 reproduces the legacy midpoint approximation bit-for-bit.
	// Set to >=2 to upgrade W^(1)_WP to a proper cell-averaged matrix element via
	// Gauss-Legendre quadrature, matching the convention V_WP already uses with Np_per_WP.
	// / [CN] 构建 W^(1) WP bin 矩阵元（3NF 缓存）时每个 cell 的 Gauss 积分点数。默认 1
	// 等价于现行 midpoint 近似，逐比特复现旧结果；>=2 时升级为 cell-averaged 矩阵元，
	// 与 V_WP 已使用 Np_per_WP 的约定保持一致。
	int         Np_per_WP_W1;
	int         Nq_per_WP_W1;
	int			channel_idx;
	int			P123_omp_num_threads;
	double 		chebyshev_t;
	double 		chebyshev_s;
	bool        parallel_run;
	bool		P123_recovery;
	bool 		tensor_force;
	bool		isospin_breaking_1S0;
	bool 		midpoint_approx;
	bool		calculate_and_store_P123;
	bool		include_breakup_channels;
	bool		trace_im_path;             // [EN] when true, dump per-stage Re/Im norms
                                         // to <output_folder>/im_path_trace.txt during a
                                         // single-shot diagnostic run. Default false. /
                                         // [CN] true 时单次诊断跑里输出每个阶段的 Re/Im
                                         // 范数到 im_path_trace.txt。默认 false。
	bool		solve_faddeev;
	bool		solve_dense;
	bool		production_run;
	bool	    parameter_walk;
	int			PSI_start;
	int			PSI_end;
	std::string potential_model;
	// [EN] Three-nucleon force selection and LECs. "none" disables 3NF entirely and leaves the kernel assembly
	// path bit-for-bit identical to the pre-3NF 2NF-only solver. / [CN] 三体力选择与低能常数。"none" 时完全关闭
	// 3NF，核组装路径与引入 3NF 之前的纯 2NF 求解器逐比特一致。
	std::string three_nucleon_force;
	double      c_D;
	double      c_E;
	double      Lambda_3NF;
	// Diagnostic knob: scale the entire W^(1) contribution by this factor before adding to the kernel.
	// Default 1.0 = physical strength. Used to empirically locate the normalization that makes the
	// Padé-accelerated Neumann iteration converge while preserving the relative structure of the
	// 3NF contribution across LECs. / [CN] 诊断用标度：在加入核之前对整个 W^(1) 乘以此系数。
	double      w1_scale;
	std::string subfolder;
	std::string p_grid_type;
	std::string p_grid_filename;
	std::string q_grid_type;
	std::string q_grid_filename;
	std::string parameter_file;
	std::string energy_input_file;
    std::string output_folder;
	std::string P123_folder;
	std::string cache_root;     // hash-keyed cache root (P123, W1)
} run_params;

/* Structs for compact storage */
typedef struct pw_table{
	int alpha_pwtable;
	int L_pwtable;
	int S_pwtable;
	int J_pwtable;
	int T_pwtable;
	int l_pwtable;
	int twoj_pwtable;
	int twoJtotal_pwtable;
	int PARtotal_pwtable;
	int twoTtotal_pwtable;
} pw_table;

//typedef struct alpha_table{
//	int alpha_alphaprime_index;
//	int alpha_index;
//	int alphaprime_index;
//} alpha_table;
//
//
//typedef struct pmesh_table{
//	int index_ptable;
//	double mesh_ptable;
//} pmesh_table;
//
//
//typedef struct qmesh_table{
//	int index_qtable;
//	double mesh_qtable;
//} qmesh_table;
//
typedef struct bounds_mesh_table{
	int    index_table;
	double mesh_table;
} bounds_mesh_table;

//typedef struct Psparse_table{
//	int index_row_Ptable;
//	int index_col_Ptable;
//	double value_Ptable;
//} Psparse_table;

#endif // TYPE_DEFS_H
