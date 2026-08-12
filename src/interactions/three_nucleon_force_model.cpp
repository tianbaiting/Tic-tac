
#include "three_nucleon_force_model.h"

#include <iostream>

#include "three_nucleon_force_none.h"
#include "three_nucleon_force_gaussian_stub.h"
#include "chiral_N2LO_3NF.h"
#include "chiral_N2LO_3NF_full_reference.h"
#include "../utils/error_management.h"
#include "constants.h"

// [EN] Helper: build the c1/c3/cD/cE approximation with c4 EXPLICITLY dropped.
// The c4 term (Epelbaum 2002 eq. 2.2-2.3) is not implemented in this model
// (see docs/three_nf_equation_contract.md §8); the factory therefore passes
// c4=0 to the constructor and prints a one-time info message so the drop is
// documented, not silent. / [CN] 构造 c1/c3/cD/cE 近似模型，显式丢弃 c4。
static chiral_N2LO_3NF* make_chiral_N2LO_approx(run_params rp,
                                                 double c1, double c3, double c4_inherited) {
	if (c4_inherited != 0.0) {
		// c4 is inherited from the 2NF LEC set but not implemented in this model.
		// Document the drop explicitly (one-time) — the model name already says
		// "c1c3cDcE_approx", so this is not a silent drop.
		static bool informed = false;
		if (!informed) {
			std::fprintf(stderr,
				"[chiral_N2LO_c1c3cDcE_approx] NOTE: c4 = %.6f GeV^-1 inherited from the "
				"2NF LEC set is DROPPED because the c_4 2PE cross-product term is not "
				"implemented in this model. The 3NF is evaluated with c1/c3/cD/cE only. "
				"This is an approximation to the full N2LO 3NF; see "
				"docs/three_nf_equation_contract.md §8.\n", c4_inherited);
			informed = true;
		}
	}
	return new chiral_N2LO_3NF(rp.c_D, rp.c_E, rp.Lambda_3NF, c1, c3, /*c4=*/0.0);
}

std::unique_ptr<three_nucleon_force_model> three_nucleon_force_model::create(run_params run_parameters){

	const std::string& model = run_parameters.three_nucleon_force;

	if (model=="none" || model==""){
		return std::make_unique<three_nucleon_force_none>();
	}

	if (model=="gaussian_stub"){
		return std::make_unique<three_nucleon_force_gaussian_stub>();
	}

	if (model=="chiral_N2LO"){
		// FAIL-CLOSED: this unqualified name is reserved for the converged,
		// scalable complete model.  A complete O(N^5) validation projector now
		// exists, but routing realistic WPCD runs to it implicitly would be both
		// prohibitively expensive and numerically unaudited.
		raise_error("three_nucleon_force='chiral_N2LO' requests the FULL N2LO 3NF, "
		            "but the converged scalable implementation is not yet available. "
		            "Use 'chiral_N2LO_full_5d_reference' only for explicitly small "
		            "validation runs, or 'chiral_N2LO_c1c3cDcE_approx' for the documented "
		            "incomplete fast model. See docs/complete_n2lo_3nf_status.md.");
	}

	if (model=="chiral_N2LO_full_5d_reference") {
		// Complete but deliberately slow direct-Jj five-angle projector.  Keep it
		// distinct from the reserved production name `chiral_N2LO`: realistic
		// WPCD cache runs require the factorised implementation and convergence
		// evidence before that name can be enabled.
		double c1 = c1_idaho_n3lo, c3 = c3_idaho_n3lo, c4 = c4_idaho_n3lo;
		const std::string& pot = run_parameters.potential_model;
		if (pot == "N2LOopt") {
			c1 = c1_n2lo_opt; c3 = c3_n2lo_opt; c4 = c4_n2lo_opt;
		} else if (pot == "Idaho_N3LO") {
			c1 = c1_idaho_n3lo; c3 = c3_idaho_n3lo; c4 = c4_idaho_n3lo;
		}
		return std::make_unique<chiral_N2LO_3NF_full_reference>(
			run_parameters.c_D, run_parameters.c_E, run_parameters.Lambda_3NF,
			c1, c3, c4, run_parameters.Nangle_3NF);
	}

	if (model=="chiral_N2LO_c1c3cDcE_approx" || model=="chiral_N2LO_without_c4"){
		// [EN] Legacy alias `chiral_N2LO_without_c4` maps to the honest
		// approximate name with a deprecation notice. / [CN] 旧别名映射到近似模型。
		if (model=="chiral_N2LO_without_c4") {
			static bool depwarned = false;
			if (!depwarned) {
				std::fprintf(stderr,
					"[deprecation] three_nucleon_force='chiral_N2LO_without_c4' is a "
					"legacy alias for 'chiral_N2LO_c1c3cDcE_approx'. Please switch to "
					"the honest name. See docs/three_nf_equation_contract.md §8.\n");
				depwarned = true;
			}
		}
		// Select c_i LECs. For chiral 2NFs, inherit the matching set; for
		// phenomenological 2NFs we fall back to the Idaho_N3LO set
		// (c1=-0.81, c3=-3.2, c4=+5.4 GeV^-1) which is the standard Witała/Nogga
		// choice when pairing chiral N2LO 3NF with a phenomenological 2NF.
		double c1 = c1_idaho_n3lo, c3 = c3_idaho_n3lo, c4 = c4_idaho_n3lo;
		const std::string& pot = run_parameters.potential_model;
		if (pot == "N2LOopt") {
			c1 = c1_n2lo_opt; c3 = c3_n2lo_opt; c4 = c4_n2lo_opt;
		} else if (pot == "Idaho_N3LO") {
			c1 = c1_idaho_n3lo; c3 = c3_idaho_n3lo; c4 = c4_idaho_n3lo;
		}
		return std::unique_ptr<three_nucleon_force_model>(
			make_chiral_N2LO_approx(run_parameters, c1, c3, c4));
	}

	std::cout << "Unknown three_nucleon_force=\"" << model << "\". "
		  << "Supported values: \"none\", \"gaussian_stub\", "
		  << "\"chiral_N2LO_c1c3cDcE_approx\", "
		  << "\"chiral_N2LO_full_5d_reference\" (complete but slow; "
		  << "\"chiral_N2LO\" remains reserved for the converged scalable model). Exiting ..." << std::endl;
	exit(-1);
	return nullptr;
}

three_nucleon_force_model* three_nucleon_force_model::fetch(run_params run_parameters){
	// [EN] Compatibility wrapper preserving the legacy raw-pointer API. New
	// production call sites should use create() and hold a unique_ptr. / [CN]
	// 兼容包装，保留旧式裸指针接口；新生产调用点应使用 create() 持有 unique_ptr。
	return create(run_parameters).release();
}
