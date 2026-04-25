// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_wp_states.cpp
// 行号区段：3..62
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
