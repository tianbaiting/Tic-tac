// Semantic RunConfig adapters (Phase 2).
//
// Bridges the legacy flat `run_params` record and the semantic `RunConfig`
// partition. `run_params` remains the authoritative storage consumed by the
// solver; these adapters are the only coupling, so the field mapping below must
// stay in exact 1:1 lock-step with type_defs.h.

#include "run_config.h"

namespace tictac::config {

RunConfig make_run_config(const run_params& rp)
{
	RunConfig cfg{};

	cfg.physics.two_body.potential_model       = rp.potential_model;
	cfg.physics.two_body.tensor_force           = rp.tensor_force;
	cfg.physics.two_body.isospin_breaking_1S0   = rp.isospin_breaking_1S0;
	cfg.physics.two_body.midpoint_approx        = rp.midpoint_approx;

	cfg.physics.three_body.three_nucleon_force  = rp.three_nucleon_force;
	cfg.physics.three_body.c_D                  = rp.c_D;
	cfg.physics.three_body.c_E                  = rp.c_E;
	cfg.physics.three_body.Lambda_3NF           = rp.Lambda_3NF;
	cfg.physics.three_body.w1_scale             = rp.w1_scale;
	cfg.physics.three_body.two_J_3NF_force_max  = rp.two_J_3NF_force_max;

	cfg.truncation.two_J_3N_max = rp.two_J_3N_max;
	cfg.truncation.J_2N_max     = rp.J_2N_max;

	cfg.grid.Np_WP = rp.Np_WP;
	cfg.grid.Nq_WP = rp.Nq_WP;
	// chebyshev_t/s are the legacy COMMON controls; the optional per-axis
	// overrides are carried separately (0 == "inherit the common value").
	cfg.grid.p.grid_type         = rp.p_grid_type;
	cfg.grid.p.grid_filename     = rp.p_grid_filename;
	cfg.grid.p.chebyshev_t       = rp.chebyshev_t;
	cfg.grid.p.chebyshev_s       = rp.chebyshev_s;
	cfg.grid.p.chebyshev_t_axis  = rp.p_chebyshev_t;
	cfg.grid.p.chebyshev_s_axis  = rp.p_chebyshev_s;
	cfg.grid.q.grid_type         = rp.q_grid_type;
	cfg.grid.q.grid_filename     = rp.q_grid_filename;
	cfg.grid.q.chebyshev_t       = rp.chebyshev_t;
	cfg.grid.q.chebyshev_s       = rp.chebyshev_s;
	cfg.grid.q.chebyshev_t_axis  = rp.q_chebyshev_t;
	cfg.grid.q.chebyshev_s_axis  = rp.q_chebyshev_s;

	cfg.quadrature.Nphi         = rp.Nphi;
	cfg.quadrature.Nx           = rp.Nx;
	cfg.quadrature.Np_per_WP    = rp.Np_per_WP;
	cfg.quadrature.Nq_per_WP    = rp.Nq_per_WP;
	cfg.quadrature.Np_per_WP_W1 = rp.Np_per_WP_W1;
	cfg.quadrature.Nq_per_WP_W1 = rp.Nq_per_WP_W1;
	cfg.quadrature.Nangle_3NF   = rp.Nangle_3NF;

	cfg.solver.pade_max_order           = rp.pade_max_order;
	cfg.solver.solve_faddeev            = rp.solve_faddeev;
	cfg.solver.solve_dense              = rp.solve_dense;
	cfg.solver.include_breakup_channels = rp.include_breakup_channels;
	cfg.solver.n0_neumann_complex_born  = rp.n0_neumann_complex_born;
	cfg.solver.deuteron_binding_only    = rp.deuteron_binding_only;
	cfg.solver.calculate_and_store_P123 = rp.calculate_and_store_P123;

	cfg.runtime.parallel_run         = rp.parallel_run;
	cfg.runtime.channel_idx          = rp.channel_idx;
	cfg.runtime.P123_omp_num_threads = rp.P123_omp_num_threads;
	cfg.runtime.P123_recovery        = rp.P123_recovery;
	cfg.runtime.production_run       = rp.production_run;
	cfg.runtime.trace_im_path        = rp.trace_im_path;
	cfg.runtime.parameter_walk       = rp.parameter_walk;
	cfg.runtime.PSI_start            = rp.PSI_start;
	cfg.runtime.PSI_end              = rp.PSI_end;
	cfg.runtime.parameter_file       = rp.parameter_file;

	cfg.io.energy_input_file = rp.energy_input_file;
	cfg.io.output_folder     = rp.output_folder;
	cfg.io.P123_folder       = rp.P123_folder;
	cfg.io.cache_root        = rp.cache_root;
	cfg.io.subfolder         = rp.subfolder;

	return cfg;
}

void apply_to(run_params& rp, const RunConfig& cfg)
{
	rp.potential_model       = cfg.physics.two_body.potential_model;
	rp.tensor_force          = cfg.physics.two_body.tensor_force;
	rp.isospin_breaking_1S0  = cfg.physics.two_body.isospin_breaking_1S0;
	rp.midpoint_approx       = cfg.physics.two_body.midpoint_approx;

	rp.three_nucleon_force = cfg.physics.three_body.three_nucleon_force;
	rp.c_D                 = cfg.physics.three_body.c_D;
	rp.c_E                 = cfg.physics.three_body.c_E;
	rp.Lambda_3NF          = cfg.physics.three_body.Lambda_3NF;
	rp.w1_scale            = cfg.physics.three_body.w1_scale;
	rp.two_J_3NF_force_max = cfg.physics.three_body.two_J_3NF_force_max;

	rp.two_J_3N_max = cfg.truncation.two_J_3N_max;
	rp.J_2N_max     = cfg.truncation.J_2N_max;

	rp.Np_WP = cfg.grid.Np_WP;
	rp.Nq_WP = cfg.grid.Nq_WP;
	// The COMMON controls are written back from the p-axis; the per-axis
	// overrides are written back unchanged (0 keeps "inherit" semantics).
	rp.chebyshev_t = cfg.grid.p.chebyshev_t;
	rp.chebyshev_s = cfg.grid.p.chebyshev_s;
	rp.p_grid_type     = cfg.grid.p.grid_type;
	rp.p_grid_filename = cfg.grid.p.grid_filename;
	rp.p_chebyshev_t   = cfg.grid.p.chebyshev_t_axis;
	rp.p_chebyshev_s   = cfg.grid.p.chebyshev_s_axis;
	rp.q_grid_type     = cfg.grid.q.grid_type;
	rp.q_grid_filename = cfg.grid.q.grid_filename;
	rp.q_chebyshev_t   = cfg.grid.q.chebyshev_t_axis;
	rp.q_chebyshev_s   = cfg.grid.q.chebyshev_s_axis;

	rp.Nphi         = cfg.quadrature.Nphi;
	rp.Nx           = cfg.quadrature.Nx;
	rp.Np_per_WP    = cfg.quadrature.Np_per_WP;
	rp.Nq_per_WP    = cfg.quadrature.Nq_per_WP;
	rp.Np_per_WP_W1 = cfg.quadrature.Np_per_WP_W1;
	rp.Nq_per_WP_W1 = cfg.quadrature.Nq_per_WP_W1;
	rp.Nangle_3NF   = cfg.quadrature.Nangle_3NF;

	rp.pade_max_order           = cfg.solver.pade_max_order;
	rp.solve_faddeev            = cfg.solver.solve_faddeev;
	rp.solve_dense              = cfg.solver.solve_dense;
	rp.include_breakup_channels = cfg.solver.include_breakup_channels;
	rp.n0_neumann_complex_born  = cfg.solver.n0_neumann_complex_born;
	rp.deuteron_binding_only    = cfg.solver.deuteron_binding_only;
	rp.calculate_and_store_P123 = cfg.solver.calculate_and_store_P123;

	rp.parallel_run         = cfg.runtime.parallel_run;
	rp.channel_idx          = cfg.runtime.channel_idx;
	rp.P123_omp_num_threads = cfg.runtime.P123_omp_num_threads;
	rp.P123_recovery        = cfg.runtime.P123_recovery;
	rp.production_run       = cfg.runtime.production_run;
	rp.trace_im_path        = cfg.runtime.trace_im_path;
	rp.parameter_walk       = cfg.runtime.parameter_walk;
	rp.PSI_start            = cfg.runtime.PSI_start;
	rp.PSI_end              = cfg.runtime.PSI_end;
	rp.parameter_file       = cfg.runtime.parameter_file;

	rp.energy_input_file = cfg.io.energy_input_file;
	rp.output_folder     = cfg.io.output_folder;
	rp.P123_folder       = cfg.io.P123_folder;
	rp.cache_root        = cfg.io.cache_root;
	rp.subfolder         = cfg.io.subfolder;
}

} // namespace tictac::config
