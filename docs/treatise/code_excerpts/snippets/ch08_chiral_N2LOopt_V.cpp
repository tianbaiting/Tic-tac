// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_N2LOopt.cpp
// 行号区段：19..43
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void chiral_N2LOopt::V(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray){
	
	double tempVarray [6];
	
	/* Empty out last set of elements */
	tempVarray[0] = 0;
	tempVarray[1] = 0;
	tempVarray[2] = 0;
	tempVarray[3] = 0;
	tempVarray[4] = 0;
	tempVarray[5] = 0;
	
	int coupledState = coupled;
    int Tz_reverse = -Tz;
	
	__idaho_chiral_potential_MOD_chp(&qi, &qo, &coupledState, &S, &J, &T, &Tz_reverse, tempVarray);

	/* Change order of array to fit rest of code */
	Varray[0] = tempVarray[0];	//S=0
	Varray[1] = tempVarray[1];	//S=1
	Varray[2] = tempVarray[3];	// Li==Lo<J
	Varray[3] = tempVarray[5];	// Li<Lo
	Varray[4] = tempVarray[4];	// Li>Lo
	Varray[5] = tempVarray[2];	// Li==Lo>J
}
