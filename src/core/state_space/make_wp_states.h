#ifndef MAKE_WP_STATES_H
#define MAKE_WP_STATES_H

#include <iostream>
#include <math.h>

#include "constants.h"
#include "type_defs.h"
#include "error_management.h"
#include "disk_io_routines.h"
#include "utils/kinetic_conversion.h"
#include "utils/gauss_legendre.h"

// [EN] These helpers construct the free wave-packet (fWP) discretization used as the finite-dimensional
// projection space for the WPCD solver. The p-grid resolves pair motion and the q-grid resolves spectator
// motion, matching the tensor-product basis discussed in the notes. / [CN] 这些辅助例程构建 WPCD 求解器使用的自由波包
// (fWP) 离散化；其中 p 网格解析成对粒子内部运动，q 网格解析旁观粒子运动，对应讲稿里讨论的张量积基底。
double p_normalization(double p0, double p1);
double q_normalization(double q0, double q1);

// [EN] Weight functions remain explicit even though the active momentum-packet scheme uses unity weights, because
// the later operator formulas are written as cell averages over these measures. / [CN] 即使当前动量波包方案使用单位权重，这些权函数仍然显式保留，
// 因为后续算符公式本质上都是对这些测度上的单元平均。
double p_weight_function(double p);
double q_weight_function(double q);

// [EN] Build the outer packet boundaries for the free pair/spectator grids. / [CN] 为自由 pair / spectator 网格构造外层波包边界。
void make_chebyshev_distribution(int N_WP,
								 double* boundary_array,
								 double scale,
								 double	sparseness_degree);

// [EN] Build the pair and spectator packet boundaries from either analytic rules or user-supplied tables.
// / [CN] 通过解析规则或用户提供的表格构造 pair 与 spectator 的波包边界。
void make_p_bin_grid(fwp_statespace& fwp_states, run_params run_parameters);
void make_q_bin_grid(fwp_statespace& fwp_states, run_params run_parameters);

// [EN] Fill the internal quadrature nodes of each packet cell when midpoint mode is disabled. / [CN] 当关闭中点模式时，为每个波包单元填充内部求积节点。
void make_p_bin_quadrature_grids(fwp_statespace& fwp_states);
void make_q_bin_quadrature_grids(fwp_statespace& fwp_states);

// [EN] Assemble the full free wave-packet basis bookkeeping reused by potential, permutation, and resolvent
// construction. / [CN] 组装势矩阵、置换矩阵和 resolvent 都会复用的自由波包基 bookkeeping。
void make_fwp_statespace(fwp_statespace& fwp_states, run_params run_parameters);

#endif // MAKE_WP_STATES_H
