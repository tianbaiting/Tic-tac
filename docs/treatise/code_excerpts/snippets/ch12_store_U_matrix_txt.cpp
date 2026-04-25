// ===============================================================
// 抽取自仓库 [current]: src/io/disk_io_routines.cpp
// 行号区段：759..913
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void store_U_matrix_elements_txt(std::complex<double>*    U_array,
								 solution_configuration solve_config,
								 channel_os_indexing chn_os_indexing,
								 run_params run_parameters,
								 swp_statespace swp_states,
							     pw_3N_statespace pw_states,
							     std::string filename){
	
	/* Make local variable for run-parameters */
	std::string potential_model = run_parameters.potential_model;
	/* Make local pointers for solution-configuration */
	double*	T_lab_array = solve_config.T_lab_array;
	double*	E_com_array = solve_config.E_com_array;
	/* Make local pointers & variables for pw-statespace */
	int	 two_J			= pw_states.two_J_3N_array[0];
	int  P_3N 			= pw_states.P_3N_array[0];
	int* L_1N_array		= pw_states.L_1N_array;
	int* two_J_1N_array = pw_states.two_J_1N_array;
	/* Make local pointers & variables for SWP-statespace */
	int	   Np_WP   = swp_states.Np_WP;
	int	   Nq_WP   = swp_states.Nq_WP;
	double E_bound = swp_states.E_bound;
	/* Make local pointers & variables for on-shell channel-indexing */
	int*   q_com_idx_array		= chn_os_indexing.q_com_idx_array;
	int*   deuteron_idx_array	= chn_os_indexing.deuteron_idx_array;
	size_t num_q_com			= (size_t) chn_os_indexing.num_T_lab;
	size_t num_deuteron_states	= (size_t) chn_os_indexing.num_deuteron_states;

	/* Open file*/
	std::ofstream result_file;
	result_file.open(filename);

	std::string parity_sgn = "";
	if      (P_3N==+1){parity_sgn = "+";}
	else if	(P_3N==-1){parity_sgn = "-";}
	else			  {raise_error("Recieved illegal parity in store_U_matrix_elements_txt().");}

	result_file << std::setprecision(8);

	result_file << "# Elastic Nd-scattering U-matrix elements in a wave-packet, Jj-scheme representation. \n";
	result_file << "# \n";
	result_file << "# The calculations below were done for: \n";
	result_file << "# JP:          " << two_J << "/2" << parity_sgn << "\n";
	result_file << "# Potential:   " << potential_model << "\n";
	result_file << "# Np:          " << Np_WP << " \n";
	result_file << "# Nq:          " << Nq_WP << " \n";
	result_file << "# Particles:   " << "nd-scattering" << "\n";
	result_file << "# Deuteron BE: " << E_bound << " MeV \n";
	result_file << "# \n";
	result_file << "# Symbol definitions: \n";
	result_file << "# Uij:   U-matrix element in MeV with row-idx i and col-idx j (NOTE in wave-packet representation!)  \n";
	result_file << "# l:     Spectator nucleon orbital angular momentum  \n";
	result_file << "# j:     Spectator nucleon total angular momentum  \n";
		result_file << "# Tlab:  Solver laboratory scattering/kinetic energy [MeV]  \n";
		result_file << "# Ecm:   Solver centre-of-mass scattering/kinetic energy [MeV] (not total relativistic energy)  \n";
	result_file << "# q_idx: Index of on-shell q-momentum bin correponding to Tlab/Ecm  \n";
	result_file << "# \n";
	result_file << "# ################################################################################################### \n";
	//result_file << "#     Name    row-idx    col-idx     l'   2*j'      l    2*j  \n";

	result_file << std::right
				<< std::scientific
				<< std::setprecision(16);

	/* Write headers for table 1 */
	std::vector<std::string> headers1 = {"Name", "row-idx", "col-idx", "l'", "2*j'", "l", "2*j"};
	result_file << "#";
	for (size_t idx_header=0; idx_header<headers1.size(); idx_header++){
		result_file << std::right << std::setw(10) << headers1[idx_header];
	}
	result_file << "\n";

	/* Loop over deuteron row-indices ("dp"=deuteron prime) */
	for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
		size_t idx_alpha_row = deuteron_idx_array[idx_d_row];
		int l_row 	  = L_1N_array[idx_alpha_row];
		int two_j_row = two_J_1N_array[idx_alpha_row];

		/* Loop over deuteron column-indices ("d"=deuteron) */
		for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
			size_t idx_alpha_col = deuteron_idx_array[idx_d_col];
			int l_col 	  = L_1N_array[idx_alpha_col];
			int two_j_col = two_J_1N_array[idx_alpha_col];

			std::string label = "U" + std::to_string(idx_d_row) + std::to_string(idx_d_col);

			result_file << " ";
			result_file << std::right << std::setw(10) << label;
			result_file << std::right << std::setw(10) << idx_d_row;
			result_file << std::right << std::setw(10) << idx_d_col;
			result_file << std::right << std::setw(10) << l_row;
			result_file << std::right << std::setw(10) << two_j_row;
			result_file << std::right << std::setw(10) << l_col;
			result_file << std::right << std::setw(10) << two_j_col;
			result_file << "\n";
		}
	}
	result_file << "# ################################################################################################### \n";

	/* Write headers for table 2 */
	std::vector<std::string> headers2 = {"Tlab [MeV]"  , "Ecm [MeV]"  , "q_idx"      ,
										 "U00 [MeV]"   , "U01 [MeV]"  , "U02 [MeV]"  ,
										 "U10 [MeV]"   , "U11 [MeV]"  , "U12 [MeV]"  ,
										 "U20 [MeV]"   , "U21 [MeV]"  , "U22 [MeV]"   };
	result_file << "#";
	for (size_t idx_header=0; idx_header<headers2.size(); idx_header++){
		if (idx_header<2){
			result_file << std::right << std::setw(24) << headers2[idx_header];
		}
		else if (idx_header==2){
			result_file << std::right << std::setw(10) << headers2[idx_header];
		}
		else{
			result_file << std::right << std::setw(50) << headers2[idx_header];
		}
	}
	result_file << "\n";

	/* Loop over on-shell q-bins */
	for (size_t q_idx=0; q_idx<num_q_com; q_idx++){
		size_t q_WP_idx = q_com_idx_array[q_idx];

		double T_lab = T_lab_array[q_idx];
		double E_com = E_com_array[q_idx];

		result_file << " ";
		result_file << std::right << std::setw(24) << T_lab;
		result_file << std::right << std::setw(24) << E_com;
		result_file << std::right << std::setw(10) << q_WP_idx;

		/* Loop over deuteron row-indices ("dp"=deuteron prime) */
		for (size_t idx_d_row=0; idx_d_row<num_deuteron_states; idx_d_row++){
			/* Loop over deuteron column-indices ("d"=deuteron) */
			for (size_t idx_d_col=0; idx_d_col<num_deuteron_states; idx_d_col++){
				size_t U_idx = idx_d_row*num_deuteron_states*num_q_com + idx_d_col*num_q_com + q_idx;
				
				std::string U_string_real = to_string_with_precision_and_sign(U_array[U_idx].real(), 16);
				std::string U_string_imag = to_string_with_precision_and_sign(U_array[U_idx].imag(), 16);
				std::string U_string = U_string_real + U_string_imag + "j";
				
				/* Append U-array element in Python numpy-style complex format */
    			result_file << std::right << std::setw(50) <<U_string;
			}
		}
		result_file << "\n";
	}

	result_file << "# ################################################################################################### \n";

	/* Close writing session */
	result_file << std::endl;
	
	/* Close files */
	result_file.close();
}
