#ifndef TICTAC_CORE_FADDEEV_SOLVER_FACADE_H
#define TICTAC_CORE_FADDEEV_SOLVER_FACADE_H

#include <string>
#include <vector>

#include "type_defs.h"

class three_nucleon_force_model;  // global-namespace forward decl (src/interactions)

namespace tictac::core {

// [EN] The elastic AGS / Faddeev scattering problem on one conserved J^pi block:
//
//     U = K + K G U ,   with   K = P V + (1 + P) W^(1) ,   A = C^T K C ,
//
// (docs/three_nf_equation_contract.md §3 / §4). This bundles the mathematical
// operands of that equation -- the permutation P, the 2NF pair potential V, the
// 3NF spectator component W^(1), the channel resolvent G, and the three nested
// bases (PW / free-WP / scattering-WP) plus the on-shell selection -- so a solve
// reads as "solve this scattering problem" rather than "consume 15 positional
// arrays".
//
// All members are non-owning aliases of storage owned by the caller (the solver
// pipeline), matching the mutability the underlying solver expects. Storage and
// evaluation order are UNCHANGED.
// / [CN] 一个守恒 J^pi 块上的弹性 AGS/Faddeev 散射问题：把方程 U=K+KGU 的算符与基打包。
struct FaddeevProblem {
	// --- Kernel operands of K = P V + (1+P) W^(1) ---
	// P = P_123 + P_132, the cyclic particle permutation (CSC sparse on the packet lattice).
	double*  permutation_values        = nullptr;  // P123_sparse_val_array
	int*     permutation_row_index     = nullptr;  // P123_sparse_row_array
	std::size_t* permutation_col_ptr   = nullptr;  // P123_sparse_col_array_csc
	std::size_t  permutation_dimension = 0;         // P123_sparse_dim

	// V: 2NF pair potential in the WP basis (uncoupled + coupled blocks; coupled may be null).
	double* V_WP_unco = nullptr;
	double* V_WP_coup = nullptr;

	// W^(1): bare spectator-1 3NF component (nullable -> 2NF-only solve).
	const ::three_nucleon_force_model* tnf = nullptr;

	// G: channel resolvent, diagonal in the SWP basis (precomputed per energy).
	cdouble* G_array = nullptr;

	// --- Bases ---
	const swp_statespace*            scattering_basis = nullptr;  // C, e_SWP
	const fwp_statespace*            packet_grid      = nullptr;  // p/q mesh
	const pw_3N_statespace*          basis            = nullptr;  // channel-block PW basis
	const channel_os_indexing*       on_shell         = nullptr;  // elastic/breakup on-shell rows
};

// [EN] Per-solve controls that are not part of the scattering operator itself.
struct FaddeevSolveOptions {
	run_params   run_parameters;
	std::string  file_identification;  // grid/Jpi handshake suffix for output files
};

// [EN] On-shell elastic (and optional breakup) amplitudes produced by the solve.
struct FaddeevResult {
	std::vector<cdouble> elastic_U;  // U^{Jpi}_{alpha' alpha} at the on-shell q bins
	std::vector<cdouble> breakup_U;  // breakup on-shell amplitudes (empty if not requested)
};

// [EN] Solve one elastic AGS/Faddeev block. Allocates the on-shell U buffers and
// forwards -- with argument order preserved verbatim -- to the existing
// solve_faddeev_equations, so the numerical path is identical to the direct call.
FaddeevResult solve_faddeev_block(const FaddeevProblem&  problem,
                                  const FaddeevSolveOptions& options);

} // namespace tictac::core

#endif // TICTAC_CORE_FADDEEV_SOLVER_FACADE_H