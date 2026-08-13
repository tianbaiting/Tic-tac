// FaddeevProblem / FaddeevResult facade implementation (Phase 4).
//
// This is a thin structural facade: it assembles the operands of one elastic
// AGS/Faddeev block solve and forwards to the existing solve_faddeev_equations
// with the SAME argument order the pipeline used. No floating-point path is
// changed; only the call is re-expressed through the semantic types.

#include "faddeev_solver_facade.h"

#include "solve_faddeev.h"
#include "three_nucleon_force_model.h"

namespace tictac::core {

FaddeevResult solve_faddeev_block(const FaddeevProblem&  problem,
                                  const FaddeevSolveOptions& options)
{
	const channel_os_indexing& on_shell = *problem.on_shell;

	FaddeevResult result;
	result.elastic_U.resize(static_cast<std::size_t>(on_shell.num_T_lab)
	                         * static_cast<std::size_t>(on_shell.num_deuteron_states)
	                         * static_cast<std::size_t>(on_shell.num_deuteron_states));

	cdouble* breakup_ptr = nullptr;
	if (options.run_parameters.include_breakup_channels) {
		result.breakup_U.resize(static_cast<std::size_t>(on_shell.num_T_lab)
		                        * static_cast<std::size_t>(on_shell.num_BU_chns));
		breakup_ptr = result.breakup_U.data();
	}

	// Forward with argument order identical to the legacy direct call.
	solve_faddeev_equations(result.elastic_U.data(),
	                        breakup_ptr,
	                        problem.G_array,
	                        problem.permutation_values,
	                        problem.permutation_row_index,
	                        problem.permutation_col_ptr,
	                        problem.permutation_dimension,
	                        problem.V_WP_unco,
	                        problem.V_WP_coup,
	                        *problem.scattering_basis,
	                        *problem.packet_grid,
	                        *problem.on_shell,
	                        *problem.basis,
	                        options.file_identification,
	                        options.run_parameters,
	                        problem.tnf);
	return result;
}

} // namespace tictac::core
