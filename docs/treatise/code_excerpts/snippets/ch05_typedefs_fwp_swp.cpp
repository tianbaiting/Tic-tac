// ===============================================================
// 抽取自仓库 [current]: include/type_defs.h
// 行号区段：40..71
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
