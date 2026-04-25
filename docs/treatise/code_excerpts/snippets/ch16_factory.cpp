// ===============================================================
// 抽取自仓库 [current]: src/interactions/three_nucleon_force_model.cpp
// 行号区段：1..41
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================

#include "three_nucleon_force_model.h"

#include <iostream>

#include "three_nucleon_force_none.h"
#include "three_nucleon_force_gaussian_stub.h"
#include "chiral_N2LO_3NF.h"
#include "constants.h"

three_nucleon_force_model* three_nucleon_force_model::fetch(run_params run_parameters){

	const std::string& model = run_parameters.three_nucleon_force;

	if (model=="none" || model==""){
		return new three_nucleon_force_none();
	}

	if (model=="gaussian_stub"){
		return new three_nucleon_force_gaussian_stub();
	}

	if (model=="chiral_N2LO"){
		// Select c₁,c₃,c₄ from the 2NF potential model
		double c1 = 0.0, c3 = 0.0, c4 = 0.0;
		const std::string& pot = run_parameters.potential_model;
		if (pot == "N2LOopt") {
			c1 = c1_n2lo_opt; c3 = c3_n2lo_opt; c4 = c4_n2lo_opt;
		} else if (pot == "Idaho_N3LO") {
			c1 = c1_idaho_n3lo; c3 = c3_idaho_n3lo; c4 = c4_idaho_n3lo;
		}
		return new chiral_N2LO_3NF(run_parameters.c_D,
								   run_parameters.c_E,
								   run_parameters.Lambda_3NF,
								   c1, c3, c4);
	}

	std::cout << "Unknown three_nucleon_force=\"" << model << "\". "
			  << "Supported values: \"none\", \"gaussian_stub\", \"chiral_N2LO\". Exiting ..." << std::endl;
	exit(-1);
}
