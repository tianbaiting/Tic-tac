#include "app/solver_pipeline.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <memory>
#include <fstream>

/* Used for printf */
#include <cstdio>
#include <cstring>

/* Time-keeping modules */
#include <chrono>
#include <ctime>

#include "set_run_parameters.h"
#include "make_pw_symm_states.h"
#include "make_permutation_matrix.h"
#include "utils/gauss_legendre.h"
#include "utils/kinetic_conversion.h"
#include "interactions/potential_model.h"
#include "interactions/three_nucleon_force_model.h"
#if TICTAC_USE_NEW_CACHE_LAYER
#include "io/cache_layer/cache_layer.h"
#endif
#include "make_potential_matrix.h"
#include "make_wp_states.h"
#include "make_swp_states.h"
#include "make_resolvent.h"
#include "solve_faddeev.h"
#include "disk_io_routines.h"
#include "run_organizer.h"

namespace tictac::app {
namespace {

// [EN] RAII owner for the three COO arrays built per J^pi channel. The
// underlying builders (fill_P123_arrays, unsorted_sparse_to_coo_col_major_sorter)
// hand ownership back through double/triple-pointer out-params and may reallocate
// the buffers in place; this guard therefore always owns whatever the pointers
// currently reference and frees each buffer exactly once at channel exit, matching
// the original manual delete[] while making early-return / exception paths safe.
// / [CN] 每个 J^pi 通道构造的三个 COO 数组的 RAII 持有者。底层构造函数通过
// 二重/三重指针出参返回所有权并可能就地重新分配缓冲区；本守卫始终持有指针当前
// 所指内存，在通道退出时各释放一次，等价于原来的手动 delete[]，并使提前返回 /
// 异常路径同样安全。
struct P123CooArrays {
	double* val = nullptr;
	int*    row = nullptr;
	int*    col = nullptr;

	P123CooArrays() = default;
	~P123CooArrays() {
		delete[] val;
		delete[] row;
		delete[] col;
	}
	P123CooArrays(const P123CooArrays&) = delete;
	P123CooArrays& operator=(const P123CooArrays&) = delete;
};

struct parameter_loop_bounds {
	int idx_param_min;
	int idx_param_max;
};

struct two_body_channel_layout {
	int num_2N_unco_states;
	int num_2N_coup_states;
	int V_unco_array_size;
	int V_coup_array_size;
};

// [EN] Each conserved J^pi block occupies one contiguous alpha range, so a channel solve only needs a shallow view
// into the global basis arrays rather than copying state tables. / [CN] 每个守恒的 J^pi 块都占据一个连续的 alpha 区间，因此单通道求解只需要对全局基数组做浅视图，而不必复制状态表。
pw_3N_statespace make_channel_view(const pw_3N_statespace& pw_states, int chn_3N) {
	const int idx_alpha_lower = pw_states.chn_3N_idx_array[chn_3N];
	const int idx_alpha_upper = pw_states.chn_3N_idx_array[chn_3N + 1];

	pw_3N_statespace channel_view = {};
	channel_view.Nalpha          = idx_alpha_upper - idx_alpha_lower;
	channel_view.J_2N_max        = pw_states.J_2N_max;
	channel_view.L_2N_array      = &pw_states.L_2N_array[idx_alpha_lower];
	channel_view.S_2N_array      = &pw_states.S_2N_array[idx_alpha_lower];
	channel_view.J_2N_array      = &pw_states.J_2N_array[idx_alpha_lower];
	channel_view.T_2N_array      = &pw_states.T_2N_array[idx_alpha_lower];
	channel_view.L_1N_array      = &pw_states.L_1N_array[idx_alpha_lower];
	channel_view.two_J_1N_array  = &pw_states.two_J_1N_array[idx_alpha_lower];
	channel_view.two_T_3N_array  = &pw_states.two_T_3N_array[idx_alpha_lower];
	channel_view.two_J_3N_array  = &pw_states.two_J_3N_array[idx_alpha_lower];
	channel_view.P_3N_array      = &pw_states.P_3N_array[idx_alpha_lower];
	return channel_view;
}

// [EN] Parameter walks reuse the same solve pipeline but sweep over different interaction samples. Packaging the
// loop bounds separately keeps the channel workflow identical for the single-run and multi-sample cases. / [CN]
// 参数扫描与单次计算复用同一条求解链，只是势模型样本不同；把循环边界单独打包后，单次运行和多样本运行可以共享同一套通道工作流。
parameter_loop_bounds make_parameter_loop_bounds(const run_params& run_parameters){
	parameter_loop_bounds bounds = {0, 1};
	if (run_parameters.parameter_walk){
		bounds.idx_param_min = run_parameters.PSI_start;
		bounds.idx_param_max = run_parameters.PSI_end;
	}
	return bounds;
}

// [EN] The WPCD potential/SWP storage is block-structured by distinct 2N partial waves, not by 3N channels. These
// counts determine how many uncoupled and coupled Hamiltonian blocks must be built for one parameter set. / [CN]
// WPCD 中势矩阵和 SWP 存储按不同的两体分波块组织，而不是按三体通道组织；这里的计数决定了每个参数集需要构建多少个
// 非耦合和耦合 Hamiltonian 块。
two_body_channel_layout make_two_body_channel_layout(int J_2N_max,
													 const run_params& run_parameters,
													 int Np_WP){
	two_body_channel_layout layout = {};
	if (run_parameters.tensor_force){
		layout.num_2N_unco_states = 2*(J_2N_max+1);
		layout.num_2N_coup_states = J_2N_max;
		if (run_parameters.isospin_breaking_1S0){
			layout.num_2N_unco_states -= 1;
			layout.num_2N_coup_states += 1;
		}
	}
	else{
		layout.num_2N_unco_states = 4*J_2N_max + 2;
		layout.num_2N_coup_states = 0;
		if (run_parameters.isospin_breaking_1S0){
			layout.num_2N_unco_states -= 1;
			layout.num_2N_coup_states  = 1;
		}
	}

	layout.V_unco_array_size = Np_WP * Np_WP * layout.num_2N_unco_states;
	layout.V_coup_array_size = 0;
	if (layout.num_2N_coup_states>0){
		layout.V_coup_array_size = Np_WP * Np_WP * 4 * layout.num_2N_coup_states;
	}
	return layout;
}

// [EN] The file suffix is the handshake between the core solver and the Python validation scripts: it records the
// packet grid and conserved quantum numbers needed to regroup amplitudes later. / [CN] 这个文件后缀是核心求解器与
// Python 验证脚本之间的“握手格式”，记录了后续重组振幅所需的波包网格和守恒量子数。
std::string make_file_identification(const swp_statespace& swp_states,
									 int two_J_3N,
									 int P_3N,
									 int J_2N_max){
	return   "_Np_"   + std::to_string(swp_states.Np_WP)
		   + "_Nq_"   + std::to_string(swp_states.Nq_WP)
		   + "_JP_"   + std::to_string(two_J_3N)
		   + "_"      + std::to_string(P_3N)
		   + "_Jmax_" + std::to_string(J_2N_max);
}

// [EN] Internal orchestrator that sequences the documented WPCD workflow. Split into three private stage methods
// that mirror the three responsibility layers: program-level initialization, one J^pi channel, and one parameter
// set within that channel. The call order, loop bounds, array layouts, defaults, log text, output file names and
// numerical formulas are preserved verbatim from the original main(). / [CN] 内部编排器，按文档化的 WPCD 工作流
// 顺序串接。拆成三个私有阶段方法，对应三层职责：程序级初始化、单个 J^pi 通道、通道内单个参数集。调用顺序、循环
// 边界、数组布局、默认值、日志文本、输出文件名与数值公式均与原 main() 逐字保持一致。
class SolverPipeline {
public:
	int run(int argc, char* argv[]) {
		program_start_ = std::chrono::system_clock::now();

		initialize(argc, argv);

		for (int chn_3N=0; chn_3N<pw_states_.N_chn_3N; chn_3N++){
			process_channel(chn_3N);
		}

		auto program_end = std::chrono::system_clock::now();

		std::chrono::duration<double> total_time = program_end - program_start_;
		printf("Total run-time: %.6f\n", total_time.count());

		printf("END OF RUN");

		return 0;
	}

private:
	// ---- Stage 1: program-level initialization -------------------------------------

	void initialize(int argc, char* argv[]) {
		/* Interpret command-line input and save in run_parameters */
		set_run_parameters(argc, argv, run_parameters_);

		/* Wave-packet 3N momenta */
		Np_WP_ = run_parameters_.Np_WP;
		Nq_WP_ = run_parameters_.Nq_WP;

		/* PWE truncation */
		/* Maximum (max) values for J_2N and J_3N (minimum is set to 0 and 1, respectively)*/
		J_2N_max_ = run_parameters_.J_2N_max;

		/* Potential model class pointers */
		pot_ptr_ = potential_model::create(run_parameters_);

		/* Three-nucleon force model (null object when disabled — keeps the kernel path unchanged) */
		tnf_ptr_ = three_nucleon_force_model::create(run_parameters_);
		std::cout << "Three-nucleon force model: " << tnf_ptr_->name()
				  << " (enabled=" << (tnf_ptr_->enabled() ? "true" : "false") << ")" << std::endl;
#if TICTAC_USE_NEW_CACHE_LAYER
		tictac::cache::initialize(run_parameters_.cache_root);
		auto cache_stats = tictac::cache::summary();
		std::cout << "Cache layer initialized at " << run_parameters_.cache_root
				  << " (existing entries: p123=" << cache_stats.n_p123
				  << ", w1=" << cache_stats.n_w1 << ")" << std::endl;
#endif

		/* Start of code segment for reading input model parameters (if enabled) */
		if (run_parameters_.parameter_walk==true){
			num_params_to_loop_ = run_parameters_.PSI_end - run_parameters_.PSI_start;
			read_parameter_sample_list(run_parameters_, parameter_vector_, num_model_params_, num_params_in_file_);
			if (num_params_to_loop_>num_params_in_file_){
				raise_error("PSI-range given is larger than parameter-set input file.");
			}
		}

		/* Start of code segment for state space construction */
		printf("Constructing 3N partial-wave basis ... \n");
		construct_symmetric_pw_states(pw_states_,
									  run_parameters_);

		/* Allocate deuteron-channel index-lookup arrays */
		solve_config_.deuteron_idx_arrays = new int* [pw_states_.N_chn_3N];
		solve_config_.deuteron_num_array  = new int  [pw_states_.N_chn_3N];

		make_fwp_statespace(fwp_states_, run_parameters_);

		/* Verify either P123 or U is being calculated, or both */
		if (run_parameters_.solve_faddeev==false && run_parameters_.calculate_and_store_P123==false){
			raise_error("Both solve_faddeev=false and calculate_and_store_P123==false so code has nothing to do!");
		}
	}

	// ---- Stage 2: one J^pi channel -------------------------------------------------

	void process_channel(int chn_3N) {
		// [EN] Production runs can shard the conserved J^pi blocks across processes; in serial mode every block is
		// still visited in the same order. / [CN] 生产运行可以把守恒的 J^pi 块分片到不同进程；串行模式下仍按相同顺序遍历全部通道块。
		if (run_parameters_.parallel_run==true && chn_3N!=run_parameters_.channel_idx){
			return;
		}

		// [EN] The solver follows the WPCD workflow described in the docs and in Miller et al.: build one global
		// partial-wave basis, then solve each conserved J^pi block with its own sparse permutation/resolvent data.
		// [CN] 求解器遵循文档和 Miller 等人描述的 WPCD 流程：先构建全局分波基，再为每个守恒的 J^pi 块单独生成并求解其稀疏置换/分辨算符数据。
		pw_3N_statespace pw_substates = make_channel_view(pw_states_, chn_3N);

		int two_J_3N = pw_substates.two_J_3N_array[0];
		int P_3N 	 = pw_substates.P_3N_array[0];

		printf("Working on 3N-channel J_3N=%.d/2, PAR=%.d (channel %.d of %.d) with %.d partial-wave states \n",
				two_J_3N, P_3N, chn_3N+1, pw_states_.N_chn_3N, pw_substates.Nalpha);

		/* Start of code segment for permutation matrix construction */
		P123CooArrays p123;
		size_t P123_sparse_dim = 0;

		fill_P123_arrays(&p123.val,
						 &p123.row,
						 &p123.col,
						 P123_sparse_dim,
						 run_parameters_.production_run,
						 fwp_states_,
						 pw_substates,
						 run_parameters_,
						 run_parameters_.P123_folder);

		if (run_parameters_.solve_faddeev){
			/* Convert row-major sparse format to column-major */
			printf(" - Converting P123 from COO to CSC format ... \n");
			printf("   - Converting row- to column-major format ... \n");
			size_t  dense_dim = pw_substates.Nalpha * fwp_states_.Nq_WP * fwp_states_.Np_WP;
			unsorted_sparse_to_coo_col_major_sorter(&p123.val,
													&p123.row,
													&p123.col,
													P123_sparse_dim,
													dense_dim);
			printf("     - Done \n");

			/* Convert from COO format to CSC format */
			printf("   - Converting COO to CSC format ... \n");
			std::vector<size_t> P123_sparse_col_array_csc(dense_dim+1);
			coo_to_csc_format_converter(p123.col,
										P123_sparse_col_array_csc.data(),
										P123_sparse_dim,
										dense_dim);
			printf("     - Done \n");

			printf(" - Looping through input parameter sets ... \n");
			const parameter_loop_bounds param_bounds = make_parameter_loop_bounds(run_parameters_);

			for (int idx_param_set=param_bounds.idx_param_min; idx_param_set<param_bounds.idx_param_max; idx_param_set++){
				process_parameter_set(idx_param_set, chn_3N, pw_substates, two_J_3N, P_3N,
									  p123.val, p123.row, P123_sparse_col_array_csc.data(),
									  P123_sparse_dim);
			}
		}

		/* Delete 3N-channel, model-independent arrays */
#if TICTAC_USE_NEW_CACHE_LAYER
		tictac::cache::shutdown();
#endif
	}

	// ---- Stage 3: one parameter set within a channel ------------------------------

	void process_parameter_set(int idx_param_set,
							   int chn_3N,
							   const pw_3N_statespace& pw_substates,
							   int two_J_3N,
							   int P_3N,
							   double* P123_sparse_val_array,
							   int* P123_sparse_row_array,
							   size_t* P123_sparse_col_array_csc,
							   size_t P123_sparse_dim) {
		if (run_parameters_.parameter_walk==true){
			printf("   - Setting model parameters for PSI %d ... \n", idx_param_set);
			if (idx_param_set==run_parameters_.PSI_start){
				pot_ptr_->first_parameter_sampling(true);
			}
			else{
				pot_ptr_->first_parameter_sampling(false);
			}
			pot_ptr_->update_parameters(&parameter_vector_[idx_param_set*num_model_params_]);
			printf("     - Done \n");
		}

		/* Start of code segment for potential matrix construction */
		const two_body_channel_layout two_body_layout = make_two_body_channel_layout(J_2N_max_,
																					 run_parameters_,
																					 Np_WP_);

		/* Potential matrices */
		std::vector<double> V_WP_unco_array(two_body_layout.V_unco_array_size);
		std::vector<double> V_WP_coup_array(two_body_layout.V_coup_array_size);
		double* V_WP_coup_ptr = two_body_layout.V_coup_array_size>0 ? V_WP_coup_array.data() : nullptr;

		// [EN] Steps 5 and 6 of the documented workflow are parameter-set dependent: build V in the fWP basis,
		// then diagonalize H=H0+V to obtain the interacting SWP basis and the deuteron bound packet. / [CN] 文档
		// 工作流中的第 5、6 步依赖当前势参数集：先在 fWP 基上构建 V，再对角化 H=H0+V 得到相互作用 SWP 基和氘核束缚波包。
		printf(" - Constructing 2N-potential matrices in WP basis ... \n");
		calculate_potential_matrices_array_in_WP_basis(V_WP_unco_array.data(), two_body_layout.num_2N_unco_states,
													V_WP_coup_ptr, two_body_layout.num_2N_coup_states,
												fwp_states_,
												pw_states_,
												pot_ptr_.get(),
												run_parameters_);
		printf("   - Done \n");

		/* Start of code segment for scattering wave-packets construction */
		double  E_bound = 0;

		std::vector<double> e_SWP_unco_array((Np_WP_+1)   * two_body_layout.num_2N_unco_states);
		std::vector<double> e_SWP_coup_array(2*(Np_WP_+1) * two_body_layout.num_2N_coup_states);

		std::vector<double> C_WP_unco_array(two_body_layout.V_unco_array_size);
		std::vector<double> C_WP_coup_array(two_body_layout.V_coup_array_size);
		double* e_SWP_coup_ptr = two_body_layout.num_2N_coup_states>0 ? e_SWP_coup_array.data() : nullptr;
		double* C_WP_coup_ptr = two_body_layout.V_coup_array_size>0 ? C_WP_coup_array.data() : nullptr;

		swp_statespace swp_states;

		printf(" - Constructing 2N SWPs ... \n");
		make_swp_states(e_SWP_unco_array.data(),
						e_SWP_coup_ptr,
						C_WP_unco_array.data(),
						C_WP_coup_ptr,
						V_WP_unco_array.data(),
						V_WP_coup_ptr,
						two_body_layout.num_2N_unco_states,
						two_body_layout.num_2N_coup_states,
						E_bound,
						fwp_states_,
						swp_states,
						pw_states_,
						run_parameters_);
		printf("   - Using E_bound = %.5f MeV \n", E_bound);
		printf("   - Done \n");

		/* Start of code segment for locating on-shell indices */
		/* Setup struct containing all indexing relvant for desired on-shell U-matrix solutions */
		channel_os_indexing chn_os_indexing;

		/* Locate and index on-shell bins from input energies */
		printf(" - Identifying on-shell channels for given input ... \n");
		find_on_shell_bins(solve_config_,
						   chn_os_indexing,
						   pw_substates,
						   fwp_states_,
						   swp_states,
						   run_parameters_);

		/* Locate and index deuteron channels */
		find_deuteron_channels(solve_config_,
							   pw_states_);
		printf("   - Done \n");

		/* The use of both solve_config and chn_os_indexing is redundant, hence this copying.
		*  Should be merged in future. */
		chn_os_indexing.num_T_lab			= solve_config_.num_T_lab;
		chn_os_indexing.q_com_idx_array		= solve_config_.q_com_idx_array;
		chn_os_indexing.num_deuteron_states = solve_config_.deuteron_num_array[chn_3N];
		chn_os_indexing.deuteron_idx_array  = solve_config_.deuteron_idx_arrays[chn_3N];

		/* Start of code segment for resolvent matrix (diagonal array) construction */
		// [EN] Once the SWP spectrum is known, the channel resolvent is diagonal in this basis; step 8 of the
		// workflow therefore reduces to filling one complex propagator value per (alpha, q, p) cell and energy.
		// / [CN] 一旦 SWP 谱已知，通道分辨算符在该基中就是对角的；因此工作流的第 8 步只需为每个 (alpha,q,p) 单元
		// 和能量填入一个复传播子值。
		/* Resolvent array */
		size_t dense_dim = pw_substates.Nalpha * swp_states.Np_WP * swp_states.Nq_WP;
		std::vector<cdouble> G_array(dense_dim * chn_os_indexing.num_T_lab);

		printf(" - Constructing 3N resolvents ... \n");
		for (int j=0; j<chn_os_indexing.num_T_lab; j++){
			double E = solve_config_.E_com_array[j] + swp_states.E_bound;
			calculate_resolvent_array_in_SWP_basis(&G_array[j*dense_dim],
													E,
													swp_states,
													pw_substates,
													run_parameters_);
		}
		printf("   - Done \n");

		/* Start of code segment for iterations of elastic Faddeev equations */
		std::string file_identification = make_file_identification(swp_states,
																two_J_3N,
																P_3N,
																J_2N_max_);

		std::vector<cdouble> U_array(chn_os_indexing.num_T_lab
								   * chn_os_indexing.num_deuteron_states
								   * chn_os_indexing.num_deuteron_states);

		std::vector<cdouble> U_BU_array;
		cdouble* U_BU_ptr = nullptr;
		if (run_parameters_.include_breakup_channels){
			U_BU_array.resize(chn_os_indexing.num_T_lab * chn_os_indexing.num_BU_chns);
			U_BU_ptr = U_BU_array.data();
		}

		// [EN] Step 9 is the expensive AGS solve: build the kernel action on the required on-shell rows, iterate
		// the Neumann rescattering series, and use Padé acceleration before step 10 writes the on-shell U values.
		// / [CN] 第 9 步是最昂贵的 AGS 求解：先在所需的 on-shell 行上构造核作用，再迭代 Neumann 再散射级数，
		// 然后做 Padé 加速，最后由第 10 步写出 on-shell 的 U 元素。
		printf(" - Solving Faddeev equation ... \n");
		solve_faddeev_equations(U_array.data(),
								U_BU_ptr,
								G_array.data(),
								P123_sparse_val_array,
								P123_sparse_row_array,
								P123_sparse_col_array_csc,
								P123_sparse_dim,
								V_WP_unco_array.data(),
								V_WP_coup_ptr,
								swp_states,
								fwp_states_,
								chn_os_indexing,
								pw_substates,
								file_identification,
								run_parameters_,
								tnf_ptr_.get());
		printf("   - Done \n");

		/* Start of code segment for storing on-shell U-matrix solutions */
		std::string U_mat_filename_t = run_parameters_.output_folder + "/" + "U_PW_elements"
																	+ file_identification
																	+ "_PSI_" + std::to_string(idx_param_set)
																	+ ".txt";

		store_U_matrix_elements_txt(U_array.data(),
									solve_config_,
									chn_os_indexing,
									run_parameters_,
									swp_states,
									pw_substates,
									U_mat_filename_t);

		if (run_parameters_.include_breakup_channels){
			std::string U_mat_BU_filename_t = run_parameters_.output_folder + "/" + "U_PW_breakup_elements"
																		+ file_identification
																		+ "_PSI_" + std::to_string(idx_param_set)
																		+ ".txt";

			store_U_BU_matrix_elements_txt(U_BU_ptr,
										solve_config_,
										chn_os_indexing,
										run_parameters_,
										fwp_states_,
										swp_states,
										pw_substates,
										U_mat_BU_filename_t);
		}
	}

private:
	// ---- program-level state (owned by the initialization layer) ------------------
	run_params                                       run_parameters_;
	int                                              Np_WP_ = 0;
	int                                              Nq_WP_ = 0;
	int                                              J_2N_max_ = 0;
	fwp_statespace                                   fwp_states_;
	pw_3N_statespace                                 pw_states_;
	solution_configuration                           solve_config_;
	std::unique_ptr<potential_model>                 pot_ptr_;
	std::unique_ptr<three_nucleon_force_model>        tnf_ptr_;
	std::vector<double>                              parameter_vector_;
	int                                              num_params_to_loop_ = 0;
	int                                              num_params_in_file_ = 0;
	int                                              num_model_params_ = 0;
	std::chrono::system_clock::time_point             program_start_;
};

} // namespace

int run_solver(int argc, char* argv[]) {
	SolverPipeline pipeline;
	return pipeline.run(argc, argv);
}

} // namespace tictac::app
