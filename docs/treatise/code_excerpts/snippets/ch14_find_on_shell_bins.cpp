// ===============================================================
// 抽取自仓库 [current]: src/config/run_organizer.cpp
// 行号区段：27..100
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	    /* Use q-bin midpoints to derive on-shell solver Tlab values [MeV] for bin selection */
	std::vector<size_t> q_WP_idx_vec;

	/* Pointer to either p_SWP_unco_array or p_SWP_coup_array,
	 * which is determined by whether the channel is coupled or not */
	double* e_SWP_array_ptr = NULL;

	/* Tensor forces affect indexing. This test is necessary */
	bool tensor_force_true = (run_parameters.tensor_force==true);

	/* Special condition to reduce number of on-shell calculations */
	double Eq_lower = 0;
	double Eq_upper = 0;
	double E_com    = 0;
	double q_m		= 0;
	std::vector<bool>   midpoint_idx_vector   (Nq_WP-1, false);
	std::vector<double> T_lab_midpoint_vector (Nq_WP, false);
	for (size_t q_WP_idx=0; q_WP_idx<Nq_WP; q_WP_idx++){
		Eq_lower = 0.5*(q_WP_array[q_WP_idx]   * q_WP_array[q_WP_idx])  /mu1(E_bound);
		Eq_upper = 0.5*(q_WP_array[q_WP_idx+1] * q_WP_array[q_WP_idx+1])/mu1(E_bound);
		E_com	 = 0.5*(Eq_upper + Eq_lower);
		q_m      = com_energy_to_com_q_momentum(E_com);
		T_lab_midpoint_vector[q_WP_idx] = com_momentum_to_lab_energy(q_m, E_bound);
	}

	double* energy_input_array = NULL;
	int		num_energy_input   = 0;
	read_input_energies(energy_input_array, num_energy_input, run_parameters.energy_input_file);
		
	/* Find all bins with an on-shell Tlab */
	for (size_t q_WP_idx=0; q_WP_idx<Nq_WP-1; q_WP_idx++){
		double Eq_lab_lower = T_lab_midpoint_vector[q_WP_idx];
		double Eq_lab_upper = T_lab_midpoint_vector[q_WP_idx+1];

		/* Find all bins with an on-shell Tlab */
		for (int i=0; i<num_energy_input; i++){
			double T_lab_input = energy_input_array[i];
			/* See if input energy lies between two bin mid-points */
			if (Eq_lab_lower<=T_lab_input && T_lab_input<=Eq_lab_upper){
				midpoint_idx_vector[q_WP_idx] = true;
			}
		}
	}

	/* Use on-shell midpoints to set on-shell bins */
	std::vector<bool>   bin_idx_vector   (Nq_WP, false);
	for (size_t q_WP_idx=0; q_WP_idx<Nq_WP-1; q_WP_idx++){
		if (midpoint_idx_vector[q_WP_idx]==true){
			bin_idx_vector[q_WP_idx]   = true;
			bin_idx_vector[q_WP_idx+1] = true;
		}
	}
	/* Append on-shell bin indices to q_WP_idx_vec */
	for (size_t q_WP_idx=0; q_WP_idx<Nq_WP; q_WP_idx++){
		if (bin_idx_vector[q_WP_idx]==true){
			q_WP_idx_vec.push_back(q_WP_idx);
		}
	}
	num_T_lab = q_WP_idx_vec.size();

	solve_config.T_lab_array = new double [num_T_lab];

	for (size_t Tlab_idx=0; Tlab_idx<num_T_lab; Tlab_idx++){
		size_t q_WP_idx = q_WP_idx_vec[Tlab_idx];
		double Eq_lower = 0.5*(q_WP_array[q_WP_idx]   * q_WP_array[q_WP_idx])  /mu1(E_bound);
		double Eq_upper = 0.5*(q_WP_array[q_WP_idx+1] * q_WP_array[q_WP_idx+1])/mu1(E_bound);
		double E_com 	= 0.5*(Eq_upper + Eq_lower);
		double q 	 	= com_energy_to_com_q_momentum(E_com);
		double T_lab 	= com_momentum_to_lab_energy(q, E_bound);
		solve_config.T_lab_array[Tlab_idx] = T_lab;
	}

    solve_config.num_T_lab = num_T_lab;

