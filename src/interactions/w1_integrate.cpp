// Shared exact W^(1) per-cell Gauss-Legendre integration.
//
// Extracted verbatim (same node/weight/accumulation order) from the former
// inline loop in W1_PW_cache::build, so the monolithic dense builder and the
// distributed per-block executor produce bitwise-identical W^(1) blocks.
// No physics, normalisation, or summation order is changed.

#include "w1_integrate.h"

#include "three_nucleon_force_model.h"
#include "constants.h"
#include "gauss_legendre.h"

#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <utility>
#include <vector>

#include <omp.h>

namespace tictac::interactions {

void build_per_bin_quadrature(const double* boundary_array,
                              std::size_t   N_bins,
                              int           Npts,
                              std::vector<double>& nodes_out,
                              std::vector<double>& weights_out)
{
    nodes_out.assign(N_bins * (std::size_t)Npts, 0.0);
    weights_out.assign(N_bins * (std::size_t)Npts, 0.0);

    if (Npts == 1) {
        for (std::size_t i = 0; i < N_bins; ++i) {
            double lo = boundary_array[i];
            double hi = boundary_array[i + 1];
            nodes_out[i]   = 0.5 * (lo + hi);
            weights_out[i] = hi - lo;          // single-point rule weight = bin width
        }
        return;
    }

    std::vector<double> xref(Npts), wref(Npts);
    gauss(xref.data(), wref.data(), Npts);     // nodes/weights on [-1, 1]

    for (std::size_t i = 0; i < N_bins; ++i) {
        double lo = boundary_array[i];
        double hi = boundary_array[i + 1];
        double a  = 0.5 * (hi - lo);
        double b  = 0.5 * (hi + lo);
        double*  n  = &nodes_out[i * (std::size_t)Npts];
        double*  w  = &weights_out[i * (std::size_t)Npts];
        for (int k = 0; k < Npts; ++k) {
            n[k] = a * xref[k] + b;
            w[k] = a * wref[k];
        }
    }
}

void integrate_w1_channel_blocks(
	const three_nucleon_force_model& tnf,
	const double* p_WP_array, std::size_t Np_WP,
	const double* q_WP_array, std::size_t Nq_WP,
	const pw_3N_statespace& pw_states,
	const run_params& rp,
	const std::vector<std::pair<int, int>>& channels,
	const std::vector<double*>& out_block_ptrs)
{
	if (channels.empty()) return;
	if (channels.size() != out_block_ptrs.size()) {
		throw std::runtime_error("integrate_w1_channel_blocks: channels/out ptr size mismatch");
	}

	const std::size_t m_Np = Np_WP;
	const std::size_t m_Nq = Nq_WP;

	const double inv_hbarc  = 1.0 / hbarc;
	const double inv_hbarc5 = std::pow(inv_hbarc, 5);

	const int Np_quad = std::max(1, rp.Np_per_WP_W1);
	const int Nq_quad = std::max(1, rp.Nq_per_WP_W1);
	{
		// One-time midpoint warning (process-local). Keeps the helper
		// self-contained without spamming per-block / per-worker.
		static bool warned = false;
		if (!warned && (Np_quad == 1 || Nq_quad == 1)) {
			std::fprintf(stderr,
				"[3NF quadrature] WARNING: Np_per_WP_W1=%d, Nq_per_WP_W1=%d uses "
				"the legacy midpoint rule in at least one dimension. This is a "
				"diagnostic/reproducibility setting, not evidence of a converged "
				"physical 3NF result. Compare at least N=2 and N=4.\n",
				Np_quad, Nq_quad);
			warned = true;
		}
	}

	// Per-bin Gauss nodes/weights in MeV (consumer convention).
	std::vector<double> p_nodes_MeV, p_w_MeV, q_nodes_MeV, q_w_MeV;
	build_per_bin_quadrature(p_WP_array, m_Np, Np_quad, p_nodes_MeV, p_w_MeV);
	build_per_bin_quadrature(q_WP_array, m_Nq, Nq_quad, q_nodes_MeV, q_w_MeV);

	// Pre-compute WP normalization factor 1/sqrt(Δ) per bin.
	std::vector<double> p_inv_sqrtD(m_Np), q_inv_sqrtD(m_Nq);
	for (std::size_t i = 0; i < m_Np; ++i)
		p_inv_sqrtD[i] = 1.0 / std::sqrt(p_WP_array[i + 1] - p_WP_array[i]);
	for (std::size_t i = 0; i < m_Nq; ++i)
		q_inv_sqrtD[i] = 1.0 / std::sqrt(q_WP_array[i + 1] - q_WP_array[i]);

	const std::size_t per_block = m_Nq * m_Nq * m_Np * m_Np;

	// Momentum-independent warm-up (matches the monolithic path: populate LS
	// expansions / spin tables once per channel without evaluating a full
	// transfer integral).  Done sequentially per channel before the cell loop.
	for (const auto& ch : channels) {
		tnf.prepare_W1_channel(ch.first, ch.second, pw_states);
	}

	#pragma omp parallel for schedule(dynamic)
	for (std::size_t cell = 0; cell < per_block; ++cell) {
		std::size_t remaining = cell;
		const std::size_t ipc = remaining % m_Np;
		remaining /= m_Np;
		const std::size_t ipr = remaining % m_Np;
		remaining /= m_Np;
		const std::size_t iqc = remaining % m_Nq;
		const std::size_t iqr = remaining / m_Nq;

		std::vector<double> accum(channels.size(), 0.0);
		std::vector<double> w1_values;
		for (int kqr = 0; kqr < Nq_quad; ++kqr) {
			const double q_r_MeV = q_nodes_MeV[iqr * Nq_quad + kqr];
			const double w_qr    = q_w_MeV   [iqr * Nq_quad + kqr];
			const double q_r_fm  = q_r_MeV * inv_hbarc;
			for (int kqc = 0; kqc < Nq_quad; ++kqc) {
				const double q_c_MeV = q_nodes_MeV[iqc * Nq_quad + kqc];
				const double w_qc    = q_w_MeV   [iqc * Nq_quad + kqc];
				const double q_c_fm  = q_c_MeV * inv_hbarc;
				for (int kpr = 0; kpr < Np_quad; ++kpr) {
					const double p_r_MeV = p_nodes_MeV[ipr * Np_quad + kpr];
					const double w_pr    = p_w_MeV   [ipr * Np_quad + kpr];
					const double p_r_fm  = p_r_MeV * inv_hbarc;
					for (int kpc = 0; kpc < Np_quad; ++kpc) {
						const double p_c_MeV = p_nodes_MeV[ipc * Np_quad + kpc];
						const double w_pc    = p_w_MeV   [ipc * Np_quad + kpc];
						const double p_c_fm  = p_c_MeV * inv_hbarc;

						tnf.W1_elements_for_channels(
							channels, p_r_fm, q_r_fm, p_c_fm, q_c_fm,
							pw_states, w1_values);
						const double radial_weight =
							(p_r_MeV * q_r_MeV * p_c_MeV * q_c_MeV)
							* (w_pr * w_qr * w_pc * w_qc);
						for (std::size_t index = 0; index < accum.size(); ++index) {
							accum[index] += radial_weight * w1_values[index];
						}
					}
				}
			}
		}

		const double bin_norm = p_inv_sqrtD[ipr] * q_inv_sqrtD[iqr]
		                      * p_inv_sqrtD[ipc] * q_inv_sqrtD[iqc];
		for (std::size_t index = 0; index < channels.size(); ++index) {
			out_block_ptrs[index][cell] = accum[index] * bin_norm * inv_hbarc5;
		}
	}
}

void W1BlockExecutor::compute_block(const three_nucleon_force_model& tnf,
                                    const double* p_WP_array, std::size_t Np_WP,
                                    const double* q_WP_array, std::size_t Nq_WP,
                                    const pw_3N_statespace& pw_states,
                                    const run_params& rp,
                                    int a_r, int a_c,
                                    double* out)
{
	std::vector<std::pair<int, int>> channels{{a_r, a_c}};
	std::vector<double*> ptrs{out};
	integrate_w1_channel_blocks(tnf, p_WP_array, Np_WP, q_WP_array, Nq_WP,
	                            pw_states, rp, channels, ptrs);
}

} // namespace tictac::interactions
