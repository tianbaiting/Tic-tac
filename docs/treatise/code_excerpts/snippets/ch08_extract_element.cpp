// ===============================================================
// 抽取自仓库 [current]: src/core/potential/make_potential_matrix.cpp
// 行号区段：4..48
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
double extract_potential_element_from_array(int L, int Lp, int J, int S, bool coupled, double* V_array){
	if (J<abs(S-L) || J>S+L || J<abs(S-Lp) || J>S+Lp){
		std::cout << L <<" "<< Lp <<" "<< J <<" "<< S << std::endl;
		raise_error("Encountered unphysical state in LS-solver");
	}

	double potential_element = NAN;

	/* Sign-convention on off-diagonal coupled elements */
	double sgn = -1; // -1 USED IN BOUND STATE BENCHMARK

	if (coupled){
		if (L==Lp && L<J){         // --
			potential_element =  V_array[2];
		}
		else if (L==Lp && L>J){    // ++
			potential_element =  V_array[5];
		}
		else if (L<Lp){             // -+
			potential_element = sgn*V_array[4]; // 3 USED IN BOUND STATE BENCHMARK
		}
		else{                       // +-
			potential_element = sgn*V_array[3]; //4 USED IN BOUND STATE BENCHMARK
		}
	}
	else{
		if (L==J+1){         // 3P0 or (++)
			potential_element = V_array[5];
		}
		else if (L==J-1){			// (--)
			potential_element = V_array[2];
		}
		else if (S==0){             // S=0
			potential_element = V_array[0];
		}
		else if (S==1 && L==J){    // S=1
			potential_element = V_array[1];
		}
		else{
			potential_element = 0;
		}
	}

	return potential_element;
}
