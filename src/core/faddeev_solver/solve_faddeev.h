#ifndef SOLVE_FADDEEV_H
#define SOLVE_FADDEEV_H

#include <iostream>
#include <stdlib.h>     /* div, div_t */
#include <math.h>
#include <complex>
#include <fstream>
#include <iomanip>
#include <omp.h>
//#include "mkl.h"

/* Time-keeping modules */
#include <chrono>
#include <ctime>

#include "constants.h"
#include "type_defs.h"
#include "disk_io_routines.h"
#include "error_management.h"
#include "make_pw_symm_states.h"
#include "make_permutation_matrix.h"
#include "utils/matrix_routines.h"

class three_nucleon_force_model;

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
};

void solve_faddeev_equations(cdouble*  U_array,
					   		 cdouble*  U_BU_array,
							 cdouble*  G_array,
							 double*   P123_sparse_val_array,
							 int*      P123_sparse_row_array,
							 size_t*   P123_sparse_col_array_csc,
							 size_t    P123_sparse_dim,
							 double*   V_WP_unco_array,
							 double*   V_WP_coup_array,
							 swp_statespace swp_states,
							 fwp_statespace fwp_states,
							 channel_os_indexing chn_os_indexing,
							 pw_3N_statespace pw_states,
							 std::string file_identification,
					         run_params run_parameters,
							 const three_nucleon_force_model* tnf = nullptr);

/* Create array of pointers to C^T matrices for product (C^T)PVC in row-major format */
void create_CT_row_maj_3N_pointer_array(double** CT_RM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										int  	 num_2N_unco_states,
										int  	 num_2N_coup_states,
										size_t   Np_WP,
										pw_3N_statespace pw_states,
					         			run_params run_parameters);

/* Restructures NN coupled matrix as 4 seperate matrices
 * !!! WARNING: COLUMN-MAJOR ALGORITHM !!! */
void restructure_coupled_VC_product(double* VC_product, size_t Np_WP);

/* Create array of pointers to VC-product matrices for product (C^T)PVC in column-major format*/
void create_VC_col_maj_3N_pointer_array(double** VC_CM_array,
										double*  C_WP_unco_array,
										double*  C_WP_coup_array,
										double*  V_WP_unco_array,
										double*  V_WP_coup_array,
										int  	 num_2N_unco_states,
										int  	 num_2N_coup_states,
										size_t   Np_WP,
										pw_3N_statespace pw_states,
					         			run_params run_parameters);

void PVC_col_brute_force(double*  col_array,
					     size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
					     size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
					     double** VC_CM_array,
					     double*  P123_val_array,
					     size_t*  P123_row_array,
					     int*     P123_col_array,
					     size_t   P123_dim);

void CPVC_col_brute_force(double*  col_array,
						  size_t   idx_alpha_c, size_t idx_p_c, size_t idx_q_c,
						  size_t   Nalpha,      size_t Nq_WP,   size_t Np_WP,
						  double** CT_RM_array,
						  double** VC_CM_array,
						  double*  P123_val_array,
						  size_t*  P123_row_array,
						  int*     P123_col_array,
						  size_t   P123_dim);

void PVC_col_calc_test(size_t   Nalpha,
					   size_t 	Nq_WP,
					   size_t 	Np_WP,
					   double** VC_CM_array,
					   double*  P123_sparse_val_array,
					   int*     P123_sparse_row_array,
					   size_t*  P123_sparse_col_array,
					   size_t   P123_sparse_dim);

void CPVC_col_calc_test(size_t   Nalpha,
						size_t 	 Nq_WP,
						size_t 	 Np_WP,
						double** CT_RM_array,
						double** VC_CM_array,
						double*  P123_sparse_val_array,
						int*     P123_sparse_row_array,
						size_t*  P123_sparse_col_array,
						size_t   P123_sparse_dim);

#endif // SOLVE_FADDEEV_H