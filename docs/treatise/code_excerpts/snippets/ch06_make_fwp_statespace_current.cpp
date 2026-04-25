// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_wp_states.cpp
// 行号区段：179..245
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void make_fwp_statespace(fwp_statespace& fwp_states, run_params run_parameters){
	printf("Constructing wave-packet (WP) state space ... \n");

	// [EN] The WPCD projection space is fixed once the p and q packet counts are chosen. Everything below prepares
	// the finite tensor-product basis |x_i> \otimes |\bar{x}_j> used later for V, P123 and G. / [CN] 一旦选定 p 与 q
	// 波包个数，WPCD 的投影空间就被固定下来。下面所有步骤都是在准备之后用于 V、P123 和 G 的有限维张量积基
	// |x_i> \otimes |\bar{x}_j>。

	/* Copy space dimensions from input struct */
	fwp_states.Np_WP 	 = run_parameters.Np_WP;
	fwp_states.Nq_WP 	 = run_parameters.Nq_WP;
	fwp_states.Np_per_WP = run_parameters.Np_per_WP;
	fwp_states.Nq_per_WP = run_parameters.Nq_per_WP;

	/* Temporary dimensional integer for p_array, wp_array, fp_array (and same for q_...) */
	const int Np = num_packet_points(fwp_states.Np_WP,
									 fwp_states.Np_per_WP,
									 run_parameters.midpoint_approx);
	const int Nq = num_packet_points(fwp_states.Nq_WP,
									 fwp_states.Nq_per_WP,
									 run_parameters.midpoint_approx);

	// [EN] The packet boundaries define the projection cells in momentum space. The default Chebyshev map clusters
	// cells near threshold where the kernels vary fastest, but custom boundaries can be injected without changing
	// the downstream solver. / [CN] 这些边界定义了动量空间中的投影单元。默认的 Chebyshev 映射会在阈值附近加密，
	// 因为核函数在那里变化最快；也可以注入自定义边界，而不必改动后续求解器。
	printf(" - Constructing wave-packet (WP) p-momentum bin boundaries ... \n");
	fwp_states.p_WP_array = new double [fwp_states.Np_WP+1];
	make_p_bin_grid(fwp_states, run_parameters);
	printf("   - Done \n");
	printf(" - Constructing wave-packet (WP) q-momentum bin boundaries ... \n");
	fwp_states.q_WP_array = new double [fwp_states.Nq_WP+1];
	make_q_bin_grid(fwp_states, run_parameters);
	printf("   - Done \n");

	// [EN] After the packet cells are known, operator averages are evaluated either by one midpoint per cell or by
	// an internal Gauss-Legendre quadrature mesh. This mirrors the notes: the basis is defined by cell averages, and
	// the numerical integration rule only controls how those averages are approximated. / [CN] 在确定波包单元后，
	// 算符平均要么由每个单元的一个中点近似，要么由单元内部的 Gauss-Legendre 求积网格近似。这和讲稿一致：
	// 基底由单元平均定义，而数值积分规则只决定这些平均如何近似计算。
	if (run_parameters.midpoint_approx==true){
		printf(" - Constructing p bin-midpoint mesh ... \n");
		fwp_states.p_array  = new double [Np];
		build_midpoint_mesh(fwp_states.p_array,
						    fwp_states.p_WP_array,
						    fwp_states.Np_WP);
		printf("   - Done \n");
		printf(" - Constructing q bin-midpoint mesh ... \n");
		fwp_states.q_array  = new double [Nq];
		build_midpoint_mesh(fwp_states.q_array,
						    fwp_states.q_WP_array,
						    fwp_states.Nq_WP);
		printf("   - Done \n");
	}
	else{
		printf(" - Constructing p quadrature mesh per WP, for all WPs ... \n");
		fwp_states.p_array  = new double [Np];
		fwp_states.wp_array = new double [Np];
		make_p_bin_quadrature_grids(fwp_states);
		printf("   - Done \n");
		printf(" - Constructing q quadrature mesh per WP, for all WPs ... \n");
		fwp_states.q_array  = new double [Nq];
		fwp_states.wq_array = new double [Nq];
		make_q_bin_quadrature_grids(fwp_states);
		printf("   - Done \n");
	}

