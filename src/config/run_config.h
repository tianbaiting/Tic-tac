#ifndef TICTAC_CONFIG_RUN_CONFIG_H
#define TICTAC_CONFIG_RUN_CONFIG_H

#include <string>

#include "type_defs.h"

// [EN] Semantic partition of the legacy flat `run_params` record.
//
// `run_params` (type_defs.h) mixes every control knob of the run -- physics
// model, basis truncation, packet grid, quadrature, Padé, threading, cache
// policy, and IO paths -- into one flat struct. This header reorganizes exactly
// those fields into typed sub-configurations, one per concern, so a reader can
// understand a run as
//
//     RunConfig { physics, truncation, grid, quadrature, solver, runtime, io }
//
// instead of a wall of unrelated keys.
//
// STORAGE IS UNCHANGED. `run_params` remains the authoritative record that the
// solver consumes; the key=value parser (set_run_parameters.cpp) still fills it.
// The adapters below bridge the two representations, so this is a semantic
// separation, not a flag-day rewrite. See task Phase 2.
// / [CN] 把扁平的 run_params 按物理含义重新分组为语义子配置。存储不变：run_params
// 仍是求解器消费的权威记录；下方适配器在两种表示间转换。

namespace tictac::config {

// --- Physics model ---------------------------------------------------------
struct TwoBodyForceConfig {
	std::string potential_model;
	bool        tensor_force            = true;
	bool        isospin_breaking_1S0    = true;
	bool        midpoint_approx         = false;
};

struct ThreeBodyForceConfig {
	std::string three_nucleon_force = "none";  // "none" = 2NF-only, bit-identical
	double      c_D         = 0.0;
	double      c_E         = 0.0;
	double      Lambda_3NF  = 500.0;
	// Diagnostic/fault-injection only; MUST be 1.0 in physics runs
	// (docs/three_nf_equation_contract.md §8).
	double      w1_scale    = 1.0;
	// Independent 3NF-active J cutoff: 3NF active iff two_J_3N <= this.
	// Default -1 = active in ALL solved blocks (bit-identical to pre-cutoff).
	// See docs/j3nf_truncation_design.md.
	int         two_J_3NF_force_max = -1;
};

struct PhysicsConfig {
	TwoBodyForceConfig   two_body;
	ThreeBodyForceConfig three_body;
};

// --- Basis truncation ------------------------------------------------------
struct BasisTruncation {
	int two_J_3N_max = 9;   // 2*J_3N_max (conserved 3N angular momentum cutoff)
	int J_2N_max     = 2;   // pair-system total angular momentum cutoff
};

// --- Packet (p,q) grid -----------------------------------------------------
struct PacketAxisConfig {
	std::string grid_type = "chebyshev";   // "chebyshev" or "file"
	std::string grid_filename;             // populated iff grid_type=="file"
	// Legacy common Chebyshev controls (sparseness t, scale s).
	double chebyshev_t = 1.0;
	double chebyshev_s = 200.0;
	// Optional independent per-axis controls; 0 means "inherit the common value".
	double chebyshev_t_axis = 0.0;
	double chebyshev_s_axis = 0.0;
};

struct PacketGridConfig {
	int Np_WP = 10;
	int Nq_WP = 10;
	PacketAxisConfig p;
	PacketAxisConfig q;
};

// --- Quadrature ------------------------------------------------------------
struct QuadratureConfig {
	int Nphi          = 36;   // azimuthal Gauss-Legendre nodes (permutation build)
	int Nx            = 36;   // polar Gauss-Legendre nodes   (permutation build)
	int Np_per_WP     = 6;    // p quadrature nodes per WP cell (2NF V_WP)
	int Nq_per_WP     = 6;    // q quadrature nodes per WP cell (2NF V_WP)
	int Np_per_WP_W1  = 2;    // p quadrature nodes per cell for W^(1) (3NF cache)
	int Nq_per_WP_W1  = 2;    // q quadrature nodes per cell for W^(1) (3NF cache)
	int Nangle_3NF    = 2;    // 5D angular order of the complete 3NF projector
};

// --- Solver ----------------------------------------------------------------
struct SolverConfig {
	int  pade_max_order              = 14;   // Padé [N/N] diagonal cutoff
	bool solve_faddeev               = true; // Neumann + Padé path
	bool solve_dense                 = false;// LAPACK dense reference path
	bool include_breakup_channels    = false;
	bool n0_neumann_complex_born     = true; // keep Im part of the n=0 Born term
	bool deuteron_binding_only       = false;// stop after SWP, write deuteron BE
	bool calculate_and_store_P123    = true; // build/recover the permutation cache
};

// --- Runtime / execution ---------------------------------------------------
struct RuntimeConfig {
	bool parallel_run         = false;
	int  channel_idx          = -1;   // -1 = all J^pi blocks; >=0 = single block
	int  P123_omp_num_threads = 4;
	bool P123_recovery        = false;
	bool production_run       = false;
	bool trace_im_path        = false; // off-by-default Re/Im diagnostic dump
	bool parameter_walk       = false;
	int  PSI_start            = -1;
	int  PSI_end              = -1;
	std::string parameter_file = "none";
};

// --- IO --------------------------------------------------------------------
struct IoConfig {
	std::string energy_input_file;
	std::string output_folder;
	std::string P123_folder;
	std::string cache_root;   // hash-keyed cache root (P123, W1 caches)
	std::string subfolder;
};

// --- Aggregate -------------------------------------------------------------
struct RunConfig {
	PhysicsConfig      physics;
	BasisTruncation    truncation;
	PacketGridConfig   grid;
	QuadratureConfig   quadrature;
	SolverConfig       solver;
	RuntimeConfig      runtime;
	IoConfig           io;
};

// --- Adapters between RunConfig and the legacy run_params record -----------
// `run_params` stays the authoritative storage consumed by the solver; these
// adapters are the only coupling and keep the two representations in lock-step.

RunConfig make_run_config(const run_params& rp);
void      apply_to(run_params& rp, const RunConfig& cfg);

} // namespace tictac::config

#endif // TICTAC_CONFIG_RUN_CONFIG_H
