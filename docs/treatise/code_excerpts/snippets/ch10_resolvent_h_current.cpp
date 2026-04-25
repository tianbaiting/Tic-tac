// ===============================================================
// 抽取自仓库 [current]: src/core/resolvent/make_resolvent.h
// 行号区段：11..41
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
/* COMMENTARY A:
 * Taken from Kukulin 2007 (https://doi.org/10.1007/s11232-007-0030-3), Eqs. (32), (34-35).
 * In their notation (left: program variables, right: LaTeX symbols from paper):
 * Eq_upper = \Epsilon_{j}
 * Eq_lower = \Epsilon_{j-1}
 * Ep_upper = E_{i}
 * Ep_lower = E_{i-1}
 * Dq       = \Delta_{j}
 * Dp       = \delta_{i}
 * DP       = \Delta_{+}
 * DM       = \Delta_{-}
 * D        = \Delta */

double heaviside_step_function(double val);

cdouble resolvent_bound_continuum(double E, double Eb,
								  double q_bin_upper,
								  double q_bin_lower);

cdouble resolvent_continuum_continuum(double E, double Eb,
									  double q_bin_upper,
									  double q_bin_lower,
									  double e_bin_upper,
									  double e_bin_lower);

void calculate_resolvent_array_in_SWP_basis(cdouble* G_array,
											double   E,
											swp_statespace swp_states,
											pw_3N_statespace pw_states,
											run_params run_parameters);

