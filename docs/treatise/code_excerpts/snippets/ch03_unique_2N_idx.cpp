// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_pw_symm_states.cpp
// 行号区段：4..60
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
int unique_2N_idx(int L_2N, int S_2N, int J_2N, int T_2N, bool coupled, run_params run_parameters){
	
	if (coupled==true && run_parameters.tensor_force==false){
		raise_error("Cannot have coupled states without a tensor force.");
	}

	/* If isospin-breaking in 1S0 is enables (isospin_breaking_1S0=true)
	 * we change the unique index to move 1S0 from an uncoupled to a coupled state.
	 * This is because with isospin-breaking 1S0 essentially becomes coupled via T_3N: 1/2 <-> 3/2 */
	
	bool state_1S0 = (S_2N==0 && J_2N==0 && L_2N==0 && T_2N==1);

	int unique_idx = -1;
	if (run_parameters.tensor_force==true){
		if (coupled){
			/* Unique index for all coupled states if tensor force is on */
			unique_idx = J_2N-1;

			/* We give room to 1S0 as a coupled state if isospin-breaking is enabled.
			 * All other coupled states are moved up and 1S0 gets index 0. */
			if (state_1S0==false && run_parameters.isospin_breaking_1S0==true){
				unique_idx += 1;
			}
			if (state_1S0==true  && run_parameters.isospin_breaking_1S0==true){
				unique_idx  = 0;
			}
		}
		else{
			/* Unique index for all uncoupled 2N states if tensor force is on */
			unique_idx = 2*J_2N + S_2N;

			/* We remove the uncoupled slot given to 1S0 if it is coupled */
			if (run_parameters.isospin_breaking_1S0==true){
				unique_idx -= 1;
			}
		}
	}
	else{
		if (coupled){
			/* 1S0 will be the only coupled state if the tensor force is off */
			if (run_parameters.isospin_breaking_1S0==true){
				unique_idx = 0;
			}
		}
		else{
			/* Unique index for all uncoupled 2N states if tensor force is off */
			unique_idx = J_2N*(4*(J_2N>1) + 2) + (S_2N==1) + (S_2N==1 && J_2N!=0)*(J_2N-L_2N + 1);

			/* We remove the uncoupled slot given to 1S0 if it is coupled */
			if (run_parameters.isospin_breaking_1S0==true){
				unique_idx -= 1;
			}
		}
	}

	return unique_idx;
}
