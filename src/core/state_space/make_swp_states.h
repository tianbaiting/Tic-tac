#ifndef MAKE_SWP_STATES_H
#define MAKE_SWP_STATES_H

#include <iostream>
#include <math.h>
#include <complex>
//#include "mkl.h"

#include "constants.h"
#include "type_defs.h"
#include "error_management.h"
#include <lapacke.h>
#include "disk_io_routines.h"
#include "make_pw_symm_states.h"
#include "utils/kinetic_conversion.h"

// [EN] These routines perform the WPCD basis transformation from free packets to scattering packets by
// diagonalizing H=H0+V in each distinct two-body partial wave. The resulting eigenvectors become the C matrices
// used later in the AGS kernel. / [CN] 这些例程通过在每个不同的两体分波上对角化 H=H0+V，把自由波包基变换为散射波包基；
// 得到的本征向量随后会作为 AGS 核中使用的 C 矩阵。
void diagonalize_real_symm_matrix(float *A, float *w, float *z, int N);
void diagonalize_real_symm_matrix(double *A, double *w, double *z, int N);

void construct_free_hamiltonian(double* H0_WP_array,
								int Np_WP, double* p_WP_array);

// [EN] The Hamiltonian is stored in packed upper-triangular form because the LAPACK eigensolver expects a real
// symmetric matrix in that layout. / [CN] Hamiltonian 以压缩上三角格式存储，因为 LAPACK 的实对称本征求解器就是按
// 这种布局读取矩阵的。
void construct_full_hamiltonian(double* mat_ptr_H,
								double* mat_ptr_H0,
								double* mat_ptr_V,
								int     mat_dim);

// [EN] Coupled channels are diagonalized as one 2Np x 2Np block, but later code expects the two physical branches
// to occupy contiguous packet ranges. / [CN] 耦合通道是按一个 2Np x 2Np 的整体块对角化的，但后续代码期望两条物理分支
// 分别占据连续的波包区间。
void reorder_coupled_eigenspectrum(double* eigenvalues,
								   double* eigenvectors,
								   int     Np_WP);

// [EN] Only the deuteron channel should contain one negative-energy state; other channels must stay purely in the
// continuum. / [CN] 只有氘核通道应当包含一个负能本征态，其余通道必须保持纯连续谱。
void look_for_unphysical_bound_states(double* eigenvalues,
						   			  int     mat_dim,
						   			  bool    chn_3S1,
									  double& E_bound);

// [EN] Convert the diagonalized spectrum into the SWP cell boundaries later consumed by the resolvent. / [CN]
// 把对角化谱转换成后续 resolvent 会消费的 SWP 单元边界。
void make_swp_bin_boundaries(double* eigenvalues,
							 double* e_SWP_array,
							 int	 Np_WP,
							 bool    coupled,
							 bool    chn_3S1);

// [EN] Main WPCD basis-change step for the two-body subsystems. / [CN] 两体子系统的主 WPCD 基变换步骤。
void make_swp_states(double* e_SWP_unco_array,
					 double* e_SWP_coup_array,
					 double* C_WP_unco_array,
					 double* C_WP_coup_array,
					 double* V_WP_unco_array,
					 double* V_WP_coup_array,
					 int num_2N_unco_states,
					 int num_2N_coup_states,
					 double& E_bound,
					 fwp_statespace   fwp_states,
					 swp_statespace&  swp_states,
					 pw_3N_statespace pw_states,
					 run_params run_parameters);

#endif // MAKE_SWP_STATES_H
