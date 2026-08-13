#include "make_wp_states.h"

namespace {

// [EN] In midpoint mode every packet contributes one representative momentum. Otherwise each packet contributes its
// internal Gauss-Legendre subgrid. / [CN] 在中点模式下，每个波包只贡献一个代表动量；否则每个波包贡献其内部的
// Gauss-Legendre 子网格。
int num_packet_points(int num_packets,
					  int num_points_per_packet,
					  bool midpoint_approx){
	if (midpoint_approx){
		return num_packets;
	}
	return num_packets * num_points_per_packet;
}

// [EN] Midpoints are the cheapest realization of the packet averages described in the WPCD notes. / [CN] 中点是
// WPCD 讲稿中“波包平均态”最廉价的数值实现。
void build_midpoint_mesh(double* momentum_array,
						 const double* packet_boundaries,
						 int num_packets){
	for (int idx_packet=0; idx_packet<num_packets; idx_packet++){
		momentum_array[idx_packet] = 0.5 * (packet_boundaries[idx_packet] + packet_boundaries[idx_packet+1]);
	}
}

// [EN] The legacy call pattern repeatedly invokes gauss/updateRange on the full local packet mesh. We keep that
// exact execution path to avoid changing any downstream numerical behavior. / [CN] 历史实现会对同一局部波包网格重复调用
// `gauss`/`updateRange`；这里保留这条完全相同的执行路径，避免改变任何下游数值行为。
void fill_packet_quadrature_mesh(double* momentum_array,
								 double* quadrature_weights,
								 int num_points,
								 double lower_bound,
								 double upper_bound){
	for (int idx_point=0; idx_point<num_points; idx_point++){
		gauss(momentum_array, quadrature_weights, num_points);
		updateRange_a_b(momentum_array, quadrature_weights, lower_bound, upper_bound, num_points);
	}
}

// [EN] Weight arrays encode the packet measure reused later in V and P123 matrix elements. / [CN] 权函数数组编码了
// 后续 V 与 P123 矩阵元会复用的波包测度。
void fill_weight_array(double* weight_array,
					   const double* momentum_array,
					   int num_points,
					   double (*weight_function)(double)){
	for (int idx_point=0; idx_point<num_points; idx_point++){
		weight_array[idx_point] = weight_function(momentum_array[idx_point]);
	}
}

void fill_normalization_array(double* normalization_array,
							  const double* packet_boundaries,
							  int num_packets,
							  double (*normalization_function)(double, double)){
	for (int idx_packet=0; idx_packet<num_packets; idx_packet++){
		normalization_array[idx_packet] = normalization_function(packet_boundaries[idx_packet],
																 packet_boundaries[idx_packet+1]);
	}
}

} // namespace

double p_normalization(double p0, double p1){
	//double Ep0 = p0*p0/MN;
	//double Ep1 = p1*p1/MN;
	//return sqrt(Ep1-Ep0);	// energy WPs
	return sqrt(p1-p0);	// momentum WPs
}
double q_normalization(double q0, double q1){
	//double Eq0 = com_q_momentum_to_com_energy(q0);
	//double Eq1 = com_q_momentum_to_com_energy(q1);
	//return sqrt(Eq1-Eq0);	// energy WPs
	//return sqrt( (q1*q1-q0*q0)/MN );	// energy WPs (old)
	return sqrt(q1-q0);	// momentum WPs
}

double p_weight_function(double p){
	//return sqrt(2*p/MN);	// energy WPs
	return 1;				// momentum WPs
}
double q_weight_function(double q){
	//return sqrt(2*q/MN);	// energy WPs
	return 1;				// momentum WPs
}

void make_chebyshev_distribution(int N_WP,
								 double* boundary_array,
								 double scale,
								 double	sparseness_degree){
	double tan_term = 0;
	double boundary = 0;
	
	for (int i=1; i<N_WP+1; i++){
		tan_term = tan( (2*i-1)*M_PI/(4*N_WP) ); 

		boundary = scale*pow(tan_term, sparseness_degree);
		
		/* Momentum distribution */
		boundary_array[i] = boundary;
		/* Energy distribution */
		//boundary_array[i] = com_energy_to_com_q_momentum(boundary);
	}
	boundary_array[0] = 0.0;
}

void make_p_bin_grid(fwp_statespace& fwp_states, run_params run_parameters){
	if (run_parameters.p_grid_type=="chebyshev"){
		double scale = run_parameters.p_chebyshev_s > 0.0
			? run_parameters.p_chebyshev_s : run_parameters.chebyshev_s;
		double sparseness_degree = run_parameters.p_chebyshev_t > 0.0
			? run_parameters.p_chebyshev_t : run_parameters.chebyshev_t;
		make_chebyshev_distribution(fwp_states.Np_WP, fwp_states.p_WP_array,
									scale,
									sparseness_degree);
	}
	else if (run_parameters.p_grid_type=="custom"){
		read_WP_boundaries_from_txt(fwp_states.p_WP_array, fwp_states.Np_WP, run_parameters.p_grid_filename);
	}
	else{
		raise_error("Unknown p-momentum gridtype specified.");
	}
}
void make_q_bin_grid(fwp_statespace& fwp_states, run_params run_parameters){
	if (run_parameters.q_grid_type=="chebyshev"){
		double scale = run_parameters.q_chebyshev_s > 0.0
			? run_parameters.q_chebyshev_s : run_parameters.chebyshev_s;
		double sparseness_degree = run_parameters.q_chebyshev_t > 0.0
			? run_parameters.q_chebyshev_t : run_parameters.chebyshev_t;
		make_chebyshev_distribution(fwp_states.Nq_WP, fwp_states.q_WP_array,
									scale,
									sparseness_degree);
	}
	else if (run_parameters.q_grid_type=="custom"){
		read_WP_boundaries_from_txt(fwp_states.q_WP_array, fwp_states.Nq_WP, run_parameters.q_grid_filename);
	}
	else{
		raise_error("Unknown q-momentum gridtype specified.");
	}
}

void make_p_bin_quadrature_grids(fwp_statespace& fwp_states){
	int 	Np_WP		= fwp_states.Np_WP;
	double* p_WP_array	= fwp_states.p_WP_array;
	int 	Np_per_WP	= fwp_states.Np_per_WP;
	double* p_array		= fwp_states.p_array;
	double* wp_array	= fwp_states.wp_array;

	for (int idx_bin=0; idx_bin<Np_WP; idx_bin++){
		double bin_lower_bound = p_WP_array[idx_bin];
		double bin_upper_bound = p_WP_array[idx_bin+1];

		double*  p_array_ptr =  &p_array[idx_bin*Np_per_WP];
		double* wp_array_ptr = &wp_array[idx_bin*Np_per_WP];
		fill_packet_quadrature_mesh(p_array_ptr,
									wp_array_ptr,
									Np_per_WP,
									bin_lower_bound,
									bin_upper_bound);
	}
}
void make_q_bin_quadrature_grids(fwp_statespace& fwp_states){
	int 	Nq_WP		= fwp_states.Nq_WP;
	double* q_WP_array	= fwp_states.q_WP_array;
	int 	Nq_per_WP	= fwp_states.Nq_per_WP;
	double* q_array		= fwp_states.q_array;
	double* wq_array	= fwp_states.wq_array;

	for (int idx_bin=0; idx_bin<Nq_WP; idx_bin++){
		double bin_lower_bound = q_WP_array[idx_bin];
		double bin_upper_bound = q_WP_array[idx_bin+1];

		double*  q_array_ptr =  &q_array[idx_bin*Nq_per_WP];
		double* wq_array_ptr = &wq_array[idx_bin*Nq_per_WP];
		fill_packet_quadrature_mesh(q_array_ptr,
									wq_array_ptr,
									Nq_per_WP,
									bin_lower_bound,
									bin_upper_bound);
	}
}

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

	// [EN] The packet weight and normalization arrays encode the chosen WP convention. In the present momentum-packet
	// implementation the weights are unity and the normalization reduces to the square root of the bin width, but we
	// keep the steps explicit because the energy-packet formulas live in the same interface. / [CN] 波包的权函数和归一化
	// 数组编码了当前采用的 WP 约定。对当前的动量波包实现，权函数为 1，归一化退化为区间宽度的平方根；之所以保留显式步骤，
	// 是因为同一接口也承载了注释中的能量波包公式。
	printf(" - Calculating weight-functions for p-momentum WPs ... \n");
	fwp_states.fp_array  = new double [Np];
	fill_weight_array(fwp_states.fp_array,
					  fwp_states.p_array,
					  Np,
					  p_weight_function);
	printf("   - Done \n");
	printf(" - Calculating weight-functions for q-momentum WPs ... \n");
	fwp_states.fq_array  = new double [Nq];
	fill_weight_array(fwp_states.fq_array,
					  fwp_states.q_array,
					  Nq,
					  q_weight_function);
	printf("   - Done \n");

	/* Calculate normalization constants */
	printf(" - Normalizing p-momentum WPs ... \n");
	fwp_states.norm_p_array  = new double [fwp_states.Np_WP];
	fill_normalization_array(fwp_states.norm_p_array,
							 fwp_states.p_WP_array,
							 fwp_states.Np_WP,
							 p_normalization);
	printf("   - Done \n");
	printf(" - Normalizing q-momentum WPs ... \n");
	fwp_states.norm_q_array  = new double [fwp_states.Nq_WP];
	fill_normalization_array(fwp_states.norm_q_array,
							 fwp_states.q_WP_array,
							 fwp_states.Nq_WP,
							 q_normalization);
	printf("   - Done \n");

	// [EN] The q-bin boundaries are persisted because later scripts map discrete on-shell bins back to lab energies
	// and observables. Keeping the exported grid synchronized with the solve-time grid is essential for reproducible
	// post-processing. / [CN] q-bin 边界会被写盘，因为后处理脚本需要把离散 on-shell 单元重新映射回实验室系能量
	// 和可观测量。保证导出的网格与求解时使用的网格完全一致，是结果可复现的前提。
	printf(" - Storing q boundaries to CSV-file ... \n");
	std::string q_boundaries_filename = run_parameters.output_folder + "/" + "q_boundaries_Nq_" + std::to_string(fwp_states.Nq_WP) + ".csv";
	store_q_WP_boundaries_csv(fwp_states, q_boundaries_filename);
	printf("   - Done \n");

	printf(" - Done \n");
}
