
#include "potential_model.h"

#include "constants.h"

/* Include present potentials */
#include "chiral_LO_internal.h"
#include "chiral_twobody.h"
#include "chiral_N2LOopt.h"
#include "chiral_Idaho_N3LO.h"
#include "malfliet_tjon.h"
#include "nijmegen.h"


potential_model::potential_model(){
}

std::unique_ptr<potential_model> potential_model::create(run_params run_parameters){

	std::string model = run_parameters.potential_model;

	if (model=="LO_internal"){
		auto pot_ptr = std::make_unique<chiral_LO_internal>(MN, 0, 100);

		std::cout << "Beware that the chiral_LO_internal potential precalculates the Legendre polynomials for a given range [Jmin, Jmax] of J-values. \n"
				  << "This range is by default set to Jmin=0 and Jmax=100. If you really require calculations for Jmax>100, you can change this \n"
				  << "BEFORE compilation in potential_model.cpp " << std::endl;

		return pot_ptr;
	}
	else if (model=="IS_LO" 	 ||
			 model=="IS_NLO" 	 ||
			 model=="IS_N2LO" 	 ||
			 model=="IS_N3LO"){
		/* This potential class is more complicated to call since it allows for efficient potential
		 * construction for varying parameters by storing past calculations */
		auto pot_ptr = std::make_unique<chiral_twobody>();
		pot_ptr->call_preset(model);
		pot_ptr->set_run_parameters(run_parameters);
		return pot_ptr;
	}
	else if (model=="N2LOopt"){
		auto pot_ptr = std::make_unique<chiral_N2LOopt>();
		return pot_ptr;
	}
	else if (model=="Idaho_N3LO"){
		auto pot_ptr = std::make_unique<chiral_Idaho_N3LO>();
		return pot_ptr;
	}
	else if (model=="malfliet_tjon"){
		auto pot_ptr = std::make_unique<malfliet_tjon>();
		return pot_ptr;
	}
	else if (model=="nijmegen"){
		auto pot_ptr = std::make_unique<nijmegen>();
		return pot_ptr;
	}
	else{
		std::cout << "Invalid potential model entered. Exiting ..." << std::endl;
		exit(-1);
	}
	return nullptr;
}

potential_model *potential_model::fetch_potential_ptr(run_params run_parameters){
	// [EN] Compatibility wrapper preserving the legacy raw-pointer API. New
	// production call sites should use create() and hold a unique_ptr. / [CN]
	// 兼容包装，保留旧式裸指针接口；新生产调用点应使用 create() 持有 unique_ptr。
	return create(run_parameters).release();
}
