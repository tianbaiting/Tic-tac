// ===============================================================
// 抽取自仓库 [current]: include/type_defs.h
// 行号区段：105..153
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
typedef struct run_params{
	int         two_J_3N_max;
	int         Np_WP;
	int         Nq_WP;
	int         J_2N_max;
	int         Nphi;
	int         Nx;
	int         Np_per_WP;
	int         Nq_per_WP;
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
} run_params;
