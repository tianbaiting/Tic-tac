// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_LO_internal.cpp
// 行号区段：26..67
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
void chiral_LO_internal::V(int i, int j, double &qi, double &qo, bool coupled, int &S, int &J, int &T, int &Tz, double *Varray){
	
	Varray[0] = 0;
	Varray[1] = 0;
	Varray[2] = 0;
	Varray[3] = 0;
	Varray[4] = 0;
	Varray[5] = 0;
	
	/* Pion-exhange */
	pionExchange->potential(qi, qo, J, Varray);
	
	/* Contact terms */
	if (J==0){
		Varray[0] += parameters_LO[0];
	}
	else if (J==1){
		Varray[2] += parameters_LO[1];
	}
	
	/* Minimal relativity factors */
	double Epi = sqrt(M*M + qi*qi); 	// relativistic energy of in-going particle
	double Epo = sqrt(M*M + qo*qo); 	// relativistic energy of out-going particle
	double relFactor_i = sqrt(M/Epi);	// relativistic factor of in-going particle
	double relFactor_o = sqrt(M/Epo);	// relativistic factor of out-going particle
	
	/* Regulator-functions and Fourier transform constants */
	double temp1 = qi/Lambda;
	double temp2 = temp1*temp1*temp1*temp1*temp1*temp1;	// i.e. we have (qi/Lambda)^6
	double f1 	 = exp(-temp2);
	temp1 = qo/Lambda;
	temp2 = temp1*temp1*temp1*temp1*temp1*temp1;		// i.e. we have (qo/Lambda)^6
	double f2 	 = exp(-temp2);
	double coeff = f1*f2*relFactor_i*relFactor_o/(8*pi*pi*pi);
	
	Varray[0] *= coeff;
	Varray[1] *= coeff;
	Varray[2] *= coeff;
	Varray[3] *= coeff;
	Varray[4] *= coeff;
	Varray[5] *= coeff;
}
