// RunConfig adapter round-trip test (Phase 2).
//
// Verifies that make_run_config / apply_to is a faithful, lossless 1:1 mapping
// of every run_params field: rp -> RunConfig -> rp' must reproduce rp exactly.
// This guards the semantic partition against drift from the legacy record.
//
// Returns 0 on success.

#include <cassert>
#include <cmath>
#include <cstdio>
#include <string>

#include "run_config.h"
#include "type_defs.h"

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { ++failures; \
	std::printf("FAIL: %s (line %d)\n", #cond, __LINE__); } } while (0)

static bool same_str(const std::string& a, const std::string& b) { return a == b; }

// Populate every run_params field with distinctive, non-default sentinel values.
static run_params sentinel_run_params()
{
	run_params rp{};
	rp.two_J_3N_max = 7;   rp.Np_WP = 11;   rp.Nq_WP = 13;   rp.J_2N_max = 3;
	rp.Nphi = 41;  rp.Nx = 43;  rp.Np_per_WP = 5;  rp.Nq_per_WP = 7;
	rp.Np_per_WP_W1 = 2;  rp.Nq_per_WP_W1 = 4;  rp.Nangle_3NF = 6;
	rp.pade_max_order = 24;  rp.channel_idx = 2;  rp.P123_omp_num_threads = 8;
	rp.chebyshev_t = 0.88;  rp.chebyshev_s = 325.0;
	rp.p_chebyshev_t = 1.0;  rp.p_chebyshev_s = 90.0;
	rp.q_chebyshev_t = 2.0;  rp.q_chebyshev_s = 0.0;
	rp.parallel_run = true;  rp.P123_recovery = true;
	rp.tensor_force = false;  rp.isospin_breaking_1S0 = false;  rp.midpoint_approx = true;
	rp.calculate_and_store_P123 = false;  rp.include_breakup_channels = true;
	rp.trace_im_path = true;  rp.n0_neumann_complex_born = false;
	rp.solve_faddeev = false;  rp.solve_dense = true;
	rp.deuteron_binding_only = true;  rp.production_run = true;
	rp.parameter_walk = true;  rp.PSI_start = 0;  rp.PSI_end = 9;
	rp.potential_model = "Idaho_N3LO";
	rp.three_nucleon_force = "chiral_N2LO_full_factorized";
	rp.c_D = -0.2;  rp.c_E = -0.205;  rp.Lambda_3NF = 500.0;  rp.w1_scale = 1.0;
	rp.subfolder = "subA";
	rp.p_grid_type = "file";  rp.p_grid_filename = "p.txt";
	rp.q_grid_type = "chebyshev";  rp.q_grid_filename = "";
	rp.parameter_file = "lec.csv";
	rp.energy_input_file = "e.txt";  rp.output_folder = "out";
	rp.P123_folder = "p123";  rp.cache_root = "cache";
	return rp;
}

static bool run_params_equal(const run_params& a, const run_params& b)
{
#define EQF(x) (std::fabs((a.x) - (b.x)) < 1e-15)
#define EQI(x) ((a.x) == (b.x))
#define EQB(x) ((a.x) == (b.x))
#define EQS(x) (same_str(a.x, b.x))
	return EQI(two_J_3N_max) && EQI(Np_WP) && EQI(Nq_WP) && EQI(J_2N_max)
	    && EQI(Nphi) && EQI(Nx) && EQI(Np_per_WP) && EQI(Nq_per_WP)
	    && EQI(Np_per_WP_W1) && EQI(Nq_per_WP_W1) && EQI(Nangle_3NF)
	    && EQI(pade_max_order) && EQI(channel_idx) && EQI(P123_omp_num_threads)
	    && EQF(chebyshev_t) && EQF(chebyshev_s)
	    && EQF(p_chebyshev_t) && EQF(p_chebyshev_s)
	    && EQF(q_chebyshev_t) && EQF(q_chebyshev_s)
	    && EQB(parallel_run) && EQB(P123_recovery)
	    && EQB(tensor_force) && EQB(isospin_breaking_1S0) && EQB(midpoint_approx)
	    && EQB(calculate_and_store_P123) && EQB(include_breakup_channels)
	    && EQB(trace_im_path) && EQB(n0_neumann_complex_born)
	    && EQB(solve_faddeev) && EQB(solve_dense)
	    && EQB(deuteron_binding_only) && EQB(production_run)
	    && EQB(parameter_walk) && EQI(PSI_start) && EQI(PSI_end)
	    && EQS(potential_model) && EQS(three_nucleon_force)
	    && EQF(c_D) && EQF(c_E) && EQF(Lambda_3NF) && EQF(w1_scale)
	    && EQS(subfolder) && EQS(p_grid_type) && EQS(p_grid_filename)
	    && EQS(q_grid_type) && EQS(q_grid_filename) && EQS(parameter_file)
	    && EQS(energy_input_file) && EQS(output_folder) && EQS(P123_folder)
	    && EQS(cache_root);
#undef EQF
#undef EQI
#undef EQB
#undef EQS
}

int main()
{
	using namespace tictac::config;

	const run_params rp0 = sentinel_run_params();

	// rp -> RunConfig -> rp' must be the identity.
	RunConfig cfg = make_run_config(rp0);
	run_params rp1{};
	apply_to(rp1, cfg);
	CHECK(run_params_equal(rp0, rp1));

	// Spot-check the semantic grouping is sensible.
	CHECK(cfg.physics.three_body.three_nucleon_force == "chiral_N2LO_full_factorized");
	CHECK(cfg.physics.two_body.potential_model == "Idaho_N3LO");
	CHECK(cfg.truncation.two_J_3N_max == 7 && cfg.truncation.J_2N_max == 3);
	CHECK(cfg.grid.Np_WP == 11 && cfg.grid.Nq_WP == 13);
	CHECK(cfg.grid.p.grid_type == "file" && cfg.grid.q.grid_type == "chebyshev");
	CHECK(cfg.quadrature.Nangle_3NF == 6);
	CHECK(cfg.solver.pade_max_order == 24);
	CHECK(cfg.runtime.channel_idx == 2);
	CHECK(cfg.io.cache_root == "cache");

	if (failures == 0) {
		std::printf("run_config: round-trip adapter faithful (all fields)\n");
		return 0;
	}
	std::printf("run_config: %d CHECK(s) FAILED\n", failures);
	return 1;
}
