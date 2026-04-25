// ===============================================================
// 抽取自仓库 [current]: src/interactions/malfliet_tjon.cpp
// 行号区段：27..66
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void malfliet_tjon::V(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray){
	
	Varray[0] = 0;
	Varray[1] = 0;
	Varray[2] = 0;
	Varray[3] = 0;
	Varray[4] = 0;
	Varray[5] = 0;
	
    /* Potential parameters.
     * "R" is for repulsive.
     * "A" is for attractive. */
    double VA_I   =  513.968; // [MeV fm]
    double VA_III =  626.885; // [MeV fm]
    double VR     = 1438.720; // [MeV fm]
    double muA    =    1.550; // [no units]
    double muR    =    3.110; // [no units]

    muA *= hbarc;             // [MeV fm]
    muR *= hbarc;             // [MeV fm]

    /* "num" - numerator
     * "den" - denominator */
    if (J==0){  // MT-I
        double num_A = (qo+qi)*(qo+qi) + muA*muA;
        double num_R = (qo+qi)*(qo+qi) + muR*muR;
        double den_A = (qo-qi)*(qo-qi) + muA*muA;
        double den_R = (qo-qi)*(qo-qi) + muR*muR;

        Varray[0] = ( (VR/hbarc)*std::log(num_R/den_R) - (VA_I/hbarc)*std::log(num_A/den_A) ) / (2*M_PI*qi*qo);
    }
    else if (J==1){  // MT-III
        double num_A = (qo+qi)*(qo+qi) + muA*muA;
        double num_R = (qo+qi)*(qo+qi) + muR*muR;
        double den_A = (qo-qi)*(qo-qi) + muA*muA;
        double den_R = (qo-qi)*(qo-qi) + muR*muR;

        Varray[2] = ( (VR/hbarc)*std::log(num_R/den_R) - (VA_III/hbarc)*std::log(num_A/den_A) ) / (2*M_PI*qi*qo);
    }
}
