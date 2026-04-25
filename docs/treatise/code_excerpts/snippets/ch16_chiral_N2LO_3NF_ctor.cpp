// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_N2LO_3NF.h
// 行号区段：38..58
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
	chiral_N2LO_3NF(double c_D, double c_E, double Lambda_3NF_MeV,
				    double c1 = 0.0, double c3 = 0.0, double c4 = 0.0)
		: m_c_D(c_D)
		, m_c_E(c_E)
		, m_Lambda(Lambda_3NF_MeV / hbarc)           // Convert MeV → fm⁻¹
		, m_Lambda_chi(700.0 / hbarc)                 // Λ_χ = 700 MeV → fm⁻¹
		, m_fpi4((fpi/hbarc)*(fpi/hbarc)*(fpi/hbarc)*(fpi/hbarc))  // (f_π/ħc)⁴ in fm⁻⁴
		, m_gA(gA)                                    // Axial coupling constant (dimensionless)
		, m_mpi_fm(mpi / hbarc)                       // Pion mass in fm⁻¹
		, m_c1(c1 * hbarc / 1000.0)                   // c₁ LEC: GeV⁻¹ → fm
		, m_c3(c3 * hbarc / 1000.0)                   // c₃ LEC: GeV⁻¹ → fm
		, m_c4(c4 * hbarc / 1000.0)                   // c₄ LEC: GeV⁻¹ → fm — reserved for tensor part
		, m_gl_x(N_GL)
		, m_gl_w(N_GL)
	{
		// [EN] Pre-compute Gauss-Legendre nodes and weights on [-1, +1].
		// These are reused for every W1_1pe_contact and W1_2pe call.
		// The gauss() function from src/utils/gauss_legendre.h fills arrays
		// of length N_GL with the standard GL abscissae and weights.
		gauss(m_gl_x.data(), m_gl_w.data(), N_GL);
	}
