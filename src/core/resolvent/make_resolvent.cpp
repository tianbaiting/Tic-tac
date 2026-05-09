
#include "make_resolvent.h"

#include <cmath>
#include <cstdio>
#include <string>

namespace {

bool is_singlet_s_wave_channel(int L_2N, int S_2N, int J_2N){
	return L_2N==0 && S_2N==0 && J_2N==0;
}

double* select_swp_energy_branch(int L_2N,
								 int S_2N,
								 int J_2N,
								 int T_2N,
								 int two_T_3N,
								 int Np_WP,
								 double* e_SWP_unco_array,
								 double* e_SWP_coup_array,
								 const run_params& run_parameters){
	const bool coupled_via_L_2N = run_parameters.tensor_force && L_2N!=J_2N && J_2N!=0;
	const bool coupled_via_T_3N = is_singlet_s_wave_channel(L_2N, S_2N, J_2N) && run_parameters.isospin_breaking_1S0;
	if (coupled_via_L_2N && coupled_via_T_3N){
		raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
	}

	if (!coupled_via_L_2N && !coupled_via_T_3N){
		const int chn_2N_idx = unique_2N_idx(L_2N, S_2N, J_2N, T_2N, false, run_parameters);
		return &e_SWP_unco_array[chn_2N_idx * (Np_WP+1)];
	}

	const int chn_2N_idx = unique_2N_idx(L_2N, S_2N, J_2N, T_2N, true, run_parameters);
	double* channel_branch_ptr = &e_SWP_coup_array[chn_2N_idx * 2*(Np_WP+1)];

	// [EN] The coupled SWP storage contains two contiguous branches. Tensor-force coupling uses the lower/upper
	// orbital branch (L<J or L>J), while isospin-breaking in 1S0 uses the T_3N branch label. / [CN] 耦合 SWP 存储里
	// 连续放着两条分支：张量力耦合时按轨道角动量分支选取（L<J 或 L>J），而 1S0 的同位旋破缺则按 T_3N 分支选取。
	if (coupled_via_L_2N){
		if (L_2N>J_2N){
			channel_branch_ptr += Np_WP + 1;
		}
		return channel_branch_ptr;
	}

	if (two_T_3N!=1){
		channel_branch_ptr += Np_WP + 1;
	}
	return channel_branch_ptr;
}

} // namespace

double heaviside_step_function(double val){
	if (val<0){
		return 0;
	}
	else{
		return 1;
	}
}

// [EN] In the SWP basis the channel resolvent becomes diagonal, and the only remaining work is to integrate the
// resolvent kernel over the packet energy cells. For a deuteron-like bound channel this yields a bound-continuum
// term with the expected logarithmic real part and step-function imaginary part. / [CN] 在 SWP 基中，通道分辨算符变为对角形式；剩余工作只是把分辨核在波包能量元上积分。对于氘核这类束缚分支，会得到具有对数实部和阶跃函数虚部的束缚-连续项。
cdouble resolvent_bound_continuum(double E, double Eb,
								  double q_bin_upper,
								  double q_bin_lower){

	double mu1 = Mn*(Mn+Mp+Eb)/(Mn + Mn + Mp + Eb);

	/* Bin boundaries in energy */
	double Eq_lower = q_bin_lower*q_bin_lower/(2*mu1);
	double Eq_upper = q_bin_upper*q_bin_upper/(2*mu1);

	/* Energy width of SWPs (D is short for Delta) for q momenta */
	double Dq = Eq_upper - Eq_lower;

	/* Real part of the BC resolvent */
	double Re_R = std::log( std::abs( (Eq_lower + Eb - E)/(Eq_upper + Eb - E) ) ) / (Dq); 
	
	/* Imaginary part of the BC resolvent */
	double Im_R = (  heaviside_step_function( Eq_upper + Eb - E )
				   - heaviside_step_function( Eq_lower + Eb - E ) ) * (-M_PI)/ (Dq); 
	
	/* Return the complex BC resolvent term Q */
	return {Re_R, Im_R};
}

// [EN] For continuum-continuum packets we average the free three-body propagator over one p-cell and one q-cell.
// This regularized cell average is one reason WPCD avoids the singular kernels of the continuous formulation.
// [CN] 对连续-连续波包，我们在一个 p 单元和一个 q 单元上对自由三体传播子取平均；这种规则化的单元平均正是 WPCD 能避开连续表述中奇异核的重要原因之一。
cdouble resolvent_continuum_continuum(double E, double Eb,
									  double q_bin_upper,
									  double q_bin_lower,
									  double e_bin_upper,
									  double e_bin_lower){

	double mu1 = Mn*(Mn+Mp+Eb)/(Mn + Mn + Mp + Eb);

	/* Bin boundaries in energy */
	double Eq_lower = q_bin_lower*q_bin_lower/(2*mu1);
	double Eq_upper = q_bin_upper*q_bin_upper/(2*mu1);

	double Ep_lower = e_bin_lower;
	double Ep_upper = e_bin_upper;

	double Ep = 0.5*(Ep_lower + Ep_upper);
	double Eq = 0.5*(Eq_lower + Eq_upper);

	/* Energy widths of SWPs (D is short for Delta) for p and q momenta */
	double Dq = Eq_upper - Eq_lower;
	double Dp = Ep_upper - Ep_lower;

	/* Temporary variables */
	double DM = 0.5*(Dp - Dq);
	double DP = 0.5*(Dp + Dq);
	double D  = Ep + Eq - E;
	
	/* Real part of the CC resolvent */
	double Re_Q = (  (D+DM) * std::log( std::abs(D+DM) )
				   + (D-DM) * std::log( std::abs(D-DM) )
				   - (D+DP) * std::log( std::abs(D+DP) )
				   - (D-DP) * std::log( std::abs(D-DP) ) ) / (Dp*Dq); 
	
	/* Imaginary part of the CC resolvent */
	double Im_Q = (  (D+DM) * heaviside_step_function( D+DM )
				   + (D-DM) * heaviside_step_function( D-DM )
				   - (D+DP) * heaviside_step_function( D+DP )
				   - (D-DP) * heaviside_step_function( D-DP ) ) * M_PI/ (Dp*Dq); 
	
	/* Return the complex CC resolvent term Q */
	return {Re_Q, Im_Q};
}

void calculate_resolvent_array_in_SWP_basis(cdouble* G_array,
											double   E,
											swp_statespace swp_states,
											pw_3N_statespace pw_states,
											run_params run_parameters){

	// [EN] The full AGS/Faddeev iteration only needs the diagonal entries of G in the SWP basis, so this routine
	// fills one complex number per (alpha, q-bin, p-bin). The later Neumann/Padé solver multiplies by these entries
	// as a cheap elementwise step between rescattering iterations. / [CN] 完整的 AGS/Faddeev 迭代在 SWP 基中只需要 G 的对角元，因此这里为每个 (alpha, q-bin, p-bin) 填充一个复数；后续 Neumann/Padé 求解器在每次再散射迭代之间只需做一次廉价的逐元素乘法。
	
	int 	Np_WP			 = swp_states.Np_WP;
	int     Nq_WP			 = swp_states.Nq_WP;
	double* e_SWP_unco_array = swp_states.e_SWP_unco_array;
	double* e_SWP_coup_array = swp_states.e_SWP_coup_array;
	double* q_WP_array		 = swp_states.q_WP_array;

	int  Nalpha			= pw_states.Nalpha;
	int* L_2N_array		= pw_states.L_2N_array;
	int* S_2N_array		= pw_states.S_2N_array;
	int* J_2N_array		= pw_states.J_2N_array;
	int* T_2N_array		= pw_states.T_2N_array;
	int* two_T_3N_array = pw_states.two_T_3N_array;

	/* Loop over states along resolvent diagonal */
	for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
	
		int L = L_2N_array[idx_alpha];
		int S = S_2N_array[idx_alpha];
		int J = J_2N_array[idx_alpha];
		int T = T_2N_array[idx_alpha];

		int two_T_3N = two_T_3N_array[idx_alpha];

		// [EN] The channel resolvent is diagonal in the SWP basis, but each alpha state must still select the correct
		// interacting p-packet branch before the analytic cell average can be evaluated. / [CN] 通道分辨算符在 SWP 基中
		// 虽然是对角的，但每个 alpha 态仍需先选中正确的相互作用 p 波包分支，才能计算解析的单元平均。
		double* e_SWP_array_ptr = select_swp_energy_branch(L,
															   S,
															   J,
															   T,
															   two_T_3N,
															   Np_WP,
															   e_SWP_unco_array,
															   e_SWP_coup_array,
															   run_parameters);

		//printf("%d %d %d %d \n", L, S, J, T);

		/* p-momentum index loop */
		for (int idx_p_bin=0; idx_p_bin<Np_WP; idx_p_bin++){

			/* Upper and lower boundaries of current p-bin (expressed in energy) */
			double e_bin_lower = e_SWP_array_ptr[idx_p_bin    ];
			double e_bin_upper = e_SWP_array_ptr[idx_p_bin + 1];

			//printf("%.4f\n", e_bin_lower);

			/* Bound state check, given by p_bin_lower if bound state exists */
			bool bound_state_exists = false;
			double Eb = 0;
			if (e_bin_lower<0){
				Eb = e_bin_lower;
				bound_state_exists = true;
			}

			/* q-momentum index loop */            
			for (int idx_q_bin=0; idx_q_bin<Nq_WP; idx_q_bin++){

				/* Upper and lower boundaries of current q-bin */
				double q_bin_lower = q_WP_array[idx_q_bin    ];
				double q_bin_upper = q_WP_array[idx_q_bin + 1];

				cdouble R = {0, 0};
				cdouble Q = {0, 0};
				if (bound_state_exists){    // Calculate bound-continuum (BC) resolvent part R
					R = resolvent_bound_continuum(E, Eb,
												  q_bin_upper,
												  q_bin_lower);
				}
				else{                       // Calculate continuum-continuum (CC) resolvent part Q
					Q = resolvent_continuum_continuum(E, Eb,
													  q_bin_upper,
													  q_bin_lower,
													  e_bin_upper,
													  e_bin_lower);
					//std::cout << Q << std::endl;
				}

				/* Use identical indexing as used in permutation matrix */
				int G_idx = idx_alpha*Nq_WP*Np_WP + idx_q_bin*Np_WP + idx_p_bin;
				G_array[G_idx] = R + Q;

				//if (e_bin_lower<0){
				//	std::cout << R << " " << Q << std::endl;
				//}
				//std::cout << R << " " << Q << std::endl;
			}
		}
	}

	// [EN] Trace hook (additive). When trace_im_path=true, dump the on-shell q-bin
	// row of (Re G, Im G) aggregated over all (alpha, p) pairs for this E. Cheap:
	// one extra full sweep, only on the diagnostic run. /
	// [CN] 追踪钩子（新增，可选）。trace_im_path=true 时，把当前 E 的 on-shell q-bin
	// 行（在所有 (alpha, p) 上聚合）的 Re/Im G 写入诊断文件。开销：1 次额外扫描，仅诊断跑生效。
	if (run_parameters.trace_im_path){
		double sum_re_bc_sq = 0.0, sum_im_bc_sq = 0.0;
		double sum_re_cc_sq = 0.0, sum_im_cc_sq = 0.0;
		for (int idx_alpha=0; idx_alpha<Nalpha; idx_alpha++){
			int L = L_2N_array[idx_alpha];
			int S = S_2N_array[idx_alpha];
			int J = J_2N_array[idx_alpha];
			int T = T_2N_array[idx_alpha];
			int two_T_3N = two_T_3N_array[idx_alpha];
			double* e_SWP_array_ptr = select_swp_energy_branch(L, S, J, T, two_T_3N,
															   Np_WP, e_SWP_unco_array,
															   e_SWP_coup_array, run_parameters);
			for (int idx_p=0; idx_p<Np_WP; idx_p++){
				bool bs = (e_SWP_array_ptr[idx_p] < 0);
				for (int idx_q=0; idx_q<Nq_WP; idx_q++){
					int G_idx = idx_alpha*Nq_WP*Np_WP + idx_q*Np_WP + idx_p;
					cdouble g = G_array[G_idx];
					if (bs){
						sum_re_bc_sq += g.real()*g.real();
						sum_im_bc_sq += g.imag()*g.imag();
					} else {
						// [EN] Aggregate all CC cells; cells whose Im is exactly zero
						// (E outside cell support) contribute nothing to Im sum but do
						// contribute to Re sum. The aggregate ratio Im/Re reveals
						// global straddle weight. /
						// [CN] 聚合所有 CC 单元；E 落在单元外时 Im 为零，仅对 Re 求和有贡献，
						// 总比例 Im/Re 反映 straddle 权重。
						sum_re_cc_sq += g.real()*g.real();
						sum_im_cc_sq += g.imag()*g.imag();
					}
				}
			}
		}
		std::string trace_path = run_parameters.output_folder + "/im_path_trace.txt";
		FILE* fp = std::fopen(trace_path.c_str(), "a");
		if (fp){
			std::fprintf(fp, "G0_BC_on_shell_q\t%.6e\t%.6e\t%.6e\n",
						 std::sqrt(sum_re_bc_sq), std::sqrt(sum_im_bc_sq),
						 std::sqrt(sum_im_bc_sq) / (std::sqrt(sum_re_bc_sq) + 1e-30));
			std::fprintf(fp, "G0_CC_straddle_aggregate\t%.6e\t%.6e\t%.6e\n",
						 std::sqrt(sum_re_cc_sq), std::sqrt(sum_im_cc_sq),
						 std::sqrt(sum_im_cc_sq) / (std::sqrt(sum_re_cc_sq) + 1e-30));
			std::fclose(fp);
		}
	}
}
