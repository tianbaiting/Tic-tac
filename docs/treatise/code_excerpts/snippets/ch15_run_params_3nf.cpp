// ===============================================================
// 抽取自仓库 [current]: src/config/set_run_parameters.cpp
// 行号区段：217..230
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	else if (option == "three_nucleon_force"){
		run_parameters.three_nucleon_force = input;
	}
	else if (option == "c_D"){
		run_parameters.c_D = std::stod(input);
	}
	else if (option == "c_E"){
		run_parameters.c_E = std::stod(input);
	}
	else if (option == "Lambda_3NF"){
		run_parameters.Lambda_3NF = std::stod(input);
	}
	else if (option == "w1_scale"){
		run_parameters.w1_scale = std::stod(input);
