// ===============================================================
// 抽取自仓库 [current]: src/core/state_space/make_pw_symm_states.cpp
// 行号区段：154..175
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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

