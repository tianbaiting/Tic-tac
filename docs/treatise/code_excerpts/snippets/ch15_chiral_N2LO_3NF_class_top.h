// ===============================================================
// 抽取自仓库 [current]: src/interactions/chiral_N2LO_3NF.h
// 行号区段：31..98
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
class chiral_N2LO_3NF : public three_nucleon_force_model
{
public:
	// [EN] Number of Gauss-Legendre quadrature points for x = cos(q̂·q̂') integration.
	// 24 points provide < 0.01% accuracy for the smooth 1/(Q²+mπ²) kernel.
	static constexpr int N_GL = 24;

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

	bool enabled() const override { return m_c_D != 0.0 || m_c_E != 0.0 || m_c1 != 0.0 || m_c3 != 0.0 || m_c4 != 0.0; }
	std::string name() const override { return "chiral_N2LO"; }

	void update_parameters(const double* parameters) override
	{
		// parameters[0] = c_D, parameters[1] = c_E (for LEC fitting walks)
		if (parameters){
			m_c_D = parameters[0];
			m_c_E = parameters[1];
		}
	}

	// [EN] Evaluate the W^(1) matrix element in the partial-wave Jacobi-momentum basis.
	// Implements all three N2LO 3NF contributions with true partial-wave projection:
	//   - 3N contact term (c_E): diagonal in alpha, proportional to τ₂·τ₃
	//   - 1PE-CT (c_D): rank-0 + rank-2 decomposition of (σ₁·q̂)(σ₃·q̂), x-quadrature
	//   - 2PE (c₁, c₃): rank-0 + rank-2 decomposition of (σ₂·q₂)(σ₃·q₃), x-quadrature;
	//     off-diagonal α_r ≠ α_c (e.g. 3S1↔3D1) are now included via the rank-2 tensor.
	// The c₄ cross-product term (isospin τ₁·(τ₂×τ₃)) is deferred.
	//
	// / [CN] 计算 partial-wave Jacobi 动量基下的 W^(1) 矩阵元。包含所有三项 N2LO 3NF：
	// c_E 接触项、c_D 1PE-CT (rank-0+rank-2)、c₁/c₃ 2PE (rank-0+rank-2+off-diagonal)。
	// c₄ 叉积项推迟实现。
	double W1_element(int alpha_r, int alpha_c,
					  double p_r, double q_r,
					  double p_c, double q_c,
					  const pw_3N_statespace& pw_states) const override
	{
		// 3N conserved quantum numbers
		if (pw_states.two_J_3N_array[alpha_r] != pw_states.two_J_3N_array[alpha_c]) return 0.0;
		if (pw_states.two_T_3N_array[alpha_r] != pw_states.two_T_3N_array[alpha_c]) return 0.0;
		if (pw_states.P_3N_array[alpha_r]     != pw_states.P_3N_array[alpha_c])     return 0.0;

		double result = 0.0;
		result += W1_contact(alpha_r, alpha_c, p_r, q_r, p_c, q_c, pw_states);
		result += W1_1pe_contact(alpha_r, alpha_c, p_r, q_r, p_c, q_c, pw_states);
		result += W1_2pe(alpha_r, alpha_c, p_r, q_r, p_c, q_c, pw_states);
		return result;
	}
