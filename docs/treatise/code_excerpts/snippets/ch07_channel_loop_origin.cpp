// ===============================================================
// 抽取自仓库 [origin]: CPP/make_pw_symm_states.cpp
// 行号区段：121..230
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	/* two_J_3N loop */
	for (int two_J_3N=1; two_J_3N<two_J_3N_max+1; two_J_3N+=2){
		/* Parity loop: P_3N_remainder is the remainder of (l+L)/2 */
		for (int P_3N_remainder=0; P_3N_remainder<2; P_3N_remainder++){
			if (print_content){
				if (P_3N_remainder==0){
					printf ("Constructing channel %d with JP=%d/2-: \n", N_chn_3N_temp+1, two_J_3N);
				}
				else{
					printf ("Constructing channel %d with JP=%d/2+: \n", N_chn_3N_temp+1, two_J_3N);
				}
			}
			/* Used to cound how many PW states are in current channel (J_3N, P_3N) */
			int Nalpha_in_current_chn = 0;
			/* two_T_3N loop */
			for (int two_T_3N=1; two_T_3N<two_T_3N_max+1; two_T_3N+=2){
				/* J_2N loop */
				for (int J_2N=0; J_2N<J_2N_max+1; J_2N++){
					/* S_2N loop */
					for (int S_2N=0; S_2N<2; S_2N++){
						/* Triangle inequality for L_2N */
						L_2N_min = (int) abs(J_2N - S_2N);
						L_2N_max = J_2N + S_2N;
						/* L_2N loop */
						for (int L_2N=L_2N_min; L_2N<L_2N_max+1; L_2N++){
							/* Triangle inequality for T_2N */
							T_2N_min = (int) abs(two_T_3N - 1)/2;
							T_2N_max = (two_T_3N + 1)/2;
							/* T can't be greater than 1 */
							if (T_2N_max>1){
								T_2N_max=1;
							}
							/* T_2N loop */
							for (int T_2N=T_2N_min; T_2N<T_2N_max+1; T_2N++){
								/* Generalised Pauli exclusion principle: L_2N + S_2N + T_2N must be odd */
								if ( (L_2N + S_2N + T_2N)%2 ){
									/* Triangle inequality for J_1N */
									two_J_1N_min = (int) abs(two_J_3N - 2*J_2N);
									two_J_1N_max = two_J_3N + 2*J_2N;
									/* two_J_1N loop */
									for (int two_J_1N=two_J_1N_min; two_J_1N<two_J_1N_max+1; two_J_1N+=2){
										/* Triangle inequality for L_1N */
										L_1N_min = (int) (abs(two_J_1N - 1))/2;
										L_1N_max = (int) (two_J_1N + 1)/2;
										/* L_1N loop */
										for (int L_1N=L_1N_min; L_1N<L_1N_max+1; L_1N++){
											/* Check 3N-system total parity given by P_3N */
											if ( ((L_2N+L_1N)%2)==P_3N_remainder){

												if (two_T_3N==3){
												    if ( (S_2N==0 && L_2N==0 && J_2N==0)==false ){
												        continue;
												    }
												}

												/* Only S-waves */
												//if (L_2N!=0 || L_1N!=0){
												//	continue;
												//}

												/* Quartet channel for Malfliet-Tjon debugging */
												//if (L_2N!=0 ||
												//	S_2N!=0 ||
												//	J_2N!=0 ||
												//	T_2N!=1){
												//	continue;
												//}

												/* We've found a physical state.
												 * Append to temporary vectors */
												L_2N_temp.push_back(L_2N);
												S_2N_temp.push_back(S_2N);
												J_2N_temp.push_back(J_2N);
												T_2N_temp.push_back(T_2N);
												L_1N_temp.push_back(L_1N);
												two_J_1N_temp.push_back(two_J_1N);
												two_J_3N_temp.push_back(two_J_3N);
												two_T_3N_temp.push_back(two_T_3N);

												/* +1 if P_3N_remainder=0, -1 if P_3N_remainder=1 */
												P_3N = 1 - 2*P_3N_remainder;
												P_3N_temp.push_back(P_3N);

												/* Prints in the same order as table 1 of Glockle et al., Phys. Rep. 274 (1996) 107-285 */
												if (print_content){
													printf (print_table_format_ints, Nalpha_temp, L_2N, S_2N, J_2N, L_1N, two_J_1N,2, T_2N, two_T_3N,2, two_J_3N,2, P_3N);
												}

												/* Increment state counters */
												Nalpha_temp           += 1;
												Nalpha_in_current_chn += 1;
											}
										}
									}
								}
							}
						}
					}
				}
			}
			/* Append 3N channel if there exists states in the given channel (J_3N, T_3N, P_3N)
			 * Using Nalpha_temp_prev rather than Nalpha_temp means we get the starting index
			 * for the current channel. It makes for simpler indexing throughout the code. */
			if (Nalpha_in_current_chn!=0){
				N_chn_3N_temp += 1;
				chn_idx_temp.push_back(Nalpha_temp_prev);
			}
			Nalpha_temp_prev = Nalpha_temp;
		}
	}
