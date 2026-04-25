// ===============================================================
// 抽取自仓库 [current]: src/core/potential/make_potential_matrix.cpp
// 行号区段：110..213
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	/* Row state */
	for (int idx_alpha_r=0; idx_alpha_r<Nalpha; idx_alpha_r++){
		int L_r = L_2N_array[idx_alpha_r];
		int S_r = S_2N_array[idx_alpha_r];
		int J_r = J_2N_array[idx_alpha_r];
		int T_r = T_2N_array[idx_alpha_r];

		int two_T_3N_r = two_T_3N_array[idx_alpha_r];
		
		/* Column state */
		for (int idx_alpha_c=0; idx_alpha_c<Nalpha; idx_alpha_c++){
			int L_c = L_2N_array[idx_alpha_c];
			int S_c = S_2N_array[idx_alpha_c];
			int J_c = J_2N_array[idx_alpha_c];
			int T_c = T_2N_array[idx_alpha_c];

			int two_T_3N_c = two_T_3N_array[idx_alpha_c];

			bool check_T = T_r==T_c;
			bool check_J = J_r==J_c;
			bool check_S = S_r==S_c;
			bool check_L = ( (tensor_force_true && abs(L_r-L_c)<=2) || L_r==L_c);

			/* Check if possible channel through interaction */
			if (check_T && check_J && check_S && check_L){

				/* Type of isoscalar interaction. The numbering is ordered in the same
				 * order as eq. (160) in Glockle et al., Phys. Rep. 274, 107 (1996) */
				int isoscalar_type = -1;

				/* Determine kind of isoscalar.
				 * For now the isospin-breaking terms have special treatment.
				 * Hopefully future code-iterations will make this completely general */
				bool state_1S0 = (S_r==0 && J_r==0 && L_r==0);
				if (state_1S0 && two_T_3N_r==3 && two_T_3N_c==3){
					isoscalar_type = 2; // Interaction can be either np or nn with IS symmetry-breaking
				}
				else if (state_1S0 && (two_T_3N_r==3 && two_T_3N_c==1) || (two_T_3N_r==1 && two_T_3N_c==3)){
					isoscalar_type = 3; // Interaction can be either np or nn with IS symmetry-conservation
				}
				else if (T_r==0 && two_T_3N_r==two_T_3N_c && two_T_3N_r==1){
					isoscalar_type = 0; // Interaction must be np
				}
				else if (T_r==1 && two_T_3N_r==two_T_3N_c && two_T_3N_r==1){
					isoscalar_type = 1; // Interaction can be either np or nn with IS symmetry-conservation
				}
				else{
					raise_error("Unknown isoscalar-potential encountered in potential-matrix construction!");
				}

				/* Detemine if this is a coupled channel.
				 * !!! With isospin symmetry-breaking we count 1S0 as a coupled matrix via T_3N-coupling !!! */
				bool coupled_matrix = false;
				bool coupled_model  = false;
				bool coupled_via_L_2N = (tensor_force_true && (L_r!=L_c || (L_r==L_c && L_r!=J_r && J_r!=0)));
				bool coupled_via_T_3N = (state_1S0==true && run_parameters.isospin_breaking_1S0==true);
				if (coupled_via_L_2N && coupled_via_T_3N){
					raise_error("Warning! Code has not been written to handle isospin-breaking in coupled channels!");
				}
				if (coupled_via_L_2N || coupled_via_T_3N){ // This counts 3P0 as uncoupled; used in matrix structure
					coupled_matrix  = true;
				}
				if (L_r!=L_c || (L_r==L_c && L_r!=J_r)){   // This counts 3P0 as coupled; used in potential models
					coupled_model   = true;
				}

				/* Unique 2N-channel indices */
				int chn_idx_V_coup = 0;
				int chn_idx_V_unco = 0;

				if (coupled_matrix){
					chn_idx_V_coup = unique_2N_idx(L_r, S_r, J_r, T_r, coupled_matrix, run_parameters);
				}
				else{
					chn_idx_V_unco = unique_2N_idx(L_r, S_r, J_r, T_r, coupled_matrix, run_parameters);
				}

				/* Check if calculation has already been performed in some other alpha'-alpha interaction */
				if (coupled_matrix){
					if (matrix_calculated_coup_array[chn_idx_V_coup]==true){
						continue;
					}
					else{
						matrix_calculated_coup_array[chn_idx_V_coup] = true;
					}
				}
				else{
					if (matrix_calculated_unco_array[chn_idx_V_unco]==true){
						continue;
					}
					else{
						matrix_calculated_unco_array[chn_idx_V_unco] = true;
					}
				}

				printf("     - Working on matrix <L'=%d|V(S=%d, J=%d, T=%d)|L=%d> ", L_r, S_r, J_r, T_r, L_c);
				if (coupled_via_L_2N){
					printf("(coupled matrix via L_2N -> calculation includes all couplings of L' and L.)");
				}
				if (coupled_via_T_3N){
					printf("(coupled matrix via T_3N -> isospin-breaking.)");
				}
				printf("\n");

