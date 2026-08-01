
#include "three_nucleon_force_model.h"

#include <iostream>

#include "three_nucleon_force_none.h"
#include "three_nucleon_force_gaussian_stub.h"
#include "chiral_N2LO_3NF.h"
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

three_nucleon_force_model* three_nucleon_force_model::fetch(run_params run_parameters){

	const std::string& model = run_parameters.three_nucleon_force;

	if (model=="none" || model==""){
		return new three_nucleon_force_none();
	}

	if (model=="gaussian_stub"){
		return new three_nucleon_force_gaussian_stub();
	}

	if (model=="chiral_N2LO"){
		// [EN] FAIL-CLOSED: the string `chiral_N2LO` claims the FULL N²LO 3NF
		// (c1, c3, c4, cD, cE). Since c4 is not implemented, requesting the
		// full model is a hard error — never silently degrade to the
		// approximation. Use `chiral_N2LO_c1c3cDcE_approx` for the c1/c3/cD/cE
		// subset. / [CN] 硬阻断：`chiral_N2LO` 声称完整 N²LO 3NF，但 c4 未实现，
		// 故请求该模型为硬错误。请改用 chiral_N2LO_c1c3cDcE_approx。
		raise_error("three_nucleon_force='chiral_N2LO' requests the FULL N2LO 3NF, "
		            "but the c_4 term (Epelbaum 2002 eq. 2.2-2.3) is NOT implemented "
		            "in this build. Use three_nucleon_force='chiral_N2LO_c1c3cDcE_approx' "
		            "for the c1/c3/cD/cE subset (c4 will be dropped with a documented "
		            "notice). See docs/three_nf_equation_contract.md §8.");
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
		return make_chiral_N2LO_approx(run_parameters, c1, c3, c4);
	}

	std::cout << "Unknown three_nucleon_force=\"" << model << "\". "
			  << "Supported values: \"none\", \"gaussian_stub\", "
			  << "\"chiral_N2LO_c1c3cDcE_approx\" (c4 NOT implemented; "
			  << "\"chiral_N2LO\" is rejected as incomplete). Exiting ..." << std::endl;
	exit(-1);
}
