#ifndef CHIRAL_N2LO_3NF_H
#define CHIRAL_N2LO_3NF_H

#include "three_nucleon_force_model.h"
#include "constants.h"
#include "spin_isospin_algebra.h"
#include "chiral_3nf_pw_kernels.h"
#include "../utils/chiral_3nf_recoupling.h"
#include "gauss_legendre.h"
#include <cmath>
#include <vector>

// [EN] Chiral N2LO three-nucleon force: the leading 3NF in chiral effective field theory.
// Three contributions at this order:
//   1. Two-pion exchange (2PE):   proportional to c₁, c₃, c₄ (inherited from 2NF)
//   2. One-pion exchange contact (1PE-CT): proportional to c_D
//   3. Three-nucleon contact (CT): proportional to c_E
//
// W^(1) is the spectator-1 decomposition: particle 1 is the spectator, particles 2 and 3
// form the interacting pair. The full 3NF is recovered via W = (1+P) W^(1) where P = P₁₂₃ + P₁₃₂.
//
// References:
//   Epelbaum, Glöckle, Meissner, PRC 66 (2002) 064001
//   Witała et al., PRC 77 (2008) 034004
//   Hebeler et al., PRC 91 (2015) 024003
//
// / [CN] 手征 N2LO 三体核力：手征有效场论中领头阶的 3NF。此阶有三项贡献：
//   1. 双 π 交换 (2PE)：正比于 c₁, c₃, c₄（从 2NF 继承）
//   2. 单 π 交换接触 (1PE-CT)：正比于 c_D
//   3. 三核子接触 (CT)：正比于 c_E
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
	// Currently implements:
	//   - 3N contact term (c_E): diagonal in alpha, proportional to τ₂·τ₃
	//   - 1PE-CT scalar term (c_D): σ₁·σ₃ channel with monopole pion propagator
	//   - 2PE scalar term (c₁, c₃): two-pion exchange with monopole propagators
	//
	// The tensor (rank-2) parts of 1PE-CT and 2PE, and the c₄ cross-product term, are deferred.
	//
	// / [CN] 计算 partial-wave Jacobi 动量基下的 W^(1) 矩阵元。
	// 当前实现：3N 接触项 (c_E)、1PE-CT 标量项 (c_D)、2PE 标量项 (c₁, c₃)。
	// 1PE-CT 和 2PE 的张量(rank-2)部分以及 c₄ 叉积项推迟实现。
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

private:
	double m_c_D;
	double m_c_E;
	double m_Lambda;       // 3NF cutoff in fm⁻¹
	double m_Lambda_chi;   // Chiral symmetry breaking scale in fm⁻¹
	double m_fpi4;         // f_π⁴ in fm⁻⁴ = (f_π[MeV]/ħc)⁴
	double m_gA;           // Axial coupling constant (dimensionless)
	double m_mpi_fm;       // Pion mass in fm⁻¹
	double m_c1;           // c₁ LEC in fm (converted from GeV⁻¹)
	double m_c3;           // c₃ LEC in fm (converted from GeV⁻¹)
	double m_c4;           // c₄ LEC in fm (converted from GeV⁻¹) — reserved for tensor part
	std::vector<double> m_gl_x; // Gauss-Legendre nodes on [-1, +1]
	std::vector<double> m_gl_w; // Gauss-Legendre weights on [-1, +1]

	// [EN] 3N contact term (c_E): partial-wave matrix element assembled from
	// the factored pieces
	//   W^(1)_CT(α',p',q'; α,p,q)
	//     = A_rank0(α', α)                      [recoupling_3nf_scalar]
	//     × [ -½ c_E / (f_π⁴ Λ_χ) ]              [kernel_contact]
	//     × f_R(p',q') f_R(p,q)                 [regulator_gauss, squared Gaussian]
	//
	// The rank-0 recoupling enforces diagonality in all pair & spectator quantum
	// numbers and returns (σ_2·σ_3)(τ_2·τ_3) as the scalar eigenvalue product
	// (Task 1 helper convention; see chiral_3nf_recoupling.h).  For 3S1
	// (S_2N=1, T_2N=0) the triplet σ·σ = +1 is inert, so the c_E piece reduces
	// to the (τ_2·τ_3) eigenvalue on top of the LEC factor and regulators.
	//
	// References:
	//   Epelbaum et al. PRC 66 (2002) 064001, eqs. (2.10), (3.19), (A-4);
	//   tools/check_3nf_normalization/formula_reference.md §1.
	//
	// / [CN] 3N 接触项 (c_E)：由三部分组合的分波矩阵元——重耦合系数
	// (Task 1, 对角选择定则 + σ·σ × τ·τ)、LEC 核 (Task 2, -½ c_E / (fπ⁴ Λ_χ))、
	// 以及 E2002 eq. 3.19 的平方高斯正规化子。
	double W1_contact(int alpha_r, int alpha_c,
					  double p_r, double q_r,
					  double p_c, double q_c,
					  const pw_3N_statespace& pw_states) const
	{
		// Short-circuit when the LEC is off (pure 2NF runs still build this object).
		if (m_c_E == 0.0) return 0.0;

		// Rank-0 recoupling: pair scalar (sigma2.sigma3)(tau2.tau3) with full
		// Kronecker selection on pair and spectator quantum numbers.
		const double recoup = recoupling_3nf_scalar(
			pw_states.L_2N_array[alpha_r], pw_states.S_2N_array[alpha_r],
			pw_states.J_2N_array[alpha_r], pw_states.T_2N_array[alpha_r],
			pw_states.L_1N_array[alpha_r], pw_states.two_J_1N_array[alpha_r],
			pw_states.two_J_3N_array[alpha_r],
			pw_states.L_2N_array[alpha_c], pw_states.S_2N_array[alpha_c],
			pw_states.J_2N_array[alpha_c], pw_states.T_2N_array[alpha_c],
			pw_states.L_1N_array[alpha_c], pw_states.two_J_1N_array[alpha_c],
			pw_states.two_T_3N_array[alpha_r]);
		if (recoup == 0.0) return 0.0;

		// Squared-Gaussian regulator per E2002 eq. (3.19), applied to bra and ket.
		const double f_bra = chiral_3nf::regulator_gauss(p_r, q_r, m_Lambda);
		const double f_ket = chiral_3nf::regulator_gauss(p_c, q_c, m_Lambda);

		// LEC kernel: -½ c_E / (fπ⁴ Λ_χ) (E2002 eq. 2.10, spectator-1 component).
		const double lec = chiral_3nf::kernel_contact(m_c_E, m_fpi4, m_Lambda_chi);

		return recoup * lec * f_bra * f_ket;
	}

	// [EN] 1PE-CT term (c_D): one-pion exchange between spectator (particle 1) and pair
	// particle 3, with a contact interaction in the pair (2,3).
	//
	// Implements the true partial-wave matrix element via rank-0 + rank-2 decomposition
	// of (σ₁·q̂₃)(σ₃·q̂₃):
	//   (σ₁·q̂)(σ₃·q̂) = ⅓(σ₁·σ₃)            [rank-0]
	//                   + [σ₁⊗σ₃]₂·[q̂⊗q̂]₂   [rank-2]
	//
	// The momentum variable Q² = q² + q'² − 2qq'x is the spectator momentum transfer
	// (|Δq|² in E2002 eq. A-2), with x = cos(q̂·q̂').
	//
	// Operator structure (E2002 eq. 2.10):
	//   V^(1)_D = −(g_A D / 8) × Σ_j (τ₁·τ_j)(σ_j·q̂_j)(σ₁·q̂_j)/(|q_j|²+m_π²)
	//   D = c_D / (f_π² Λ_χ)  →  overall coeff = −g_A c_D / (8 f_π⁴ Λ_χ) × 2 (j=2,3 sum)
	//
	// Selection rules from E2002 eq. A-1 (contact pair vertex):
	//   L_2N = L_2N' = 0   (contact pair vertex → S-wave pair only)
	//   L_1N = L_1N'        (σ operators preserve spectator orbital, enforced by recoupling)
	//
	// Angular integration: 24-point Gauss-Legendre quadrature for x ∈ [−1, +1].
	//   rank-0: ∫dx (1/(8π³)) × 1/(Q²+m_π²)           [P_0(x)=1 weight]
	//   rank-2: ∫dx (1/(8π³)) × P_2(x)/(Q²+m_π²)      [P_2(x)=(3x²−1)/2 weight]
	//           Note: rank-2 vanishes for l_1N=l_1N'=0 by CG(0,0;2,0|0,0)=0.
	//
	// Note on sign: the coeff carries a minus sign from E2002 eq. (2.10). For the
	// sign of the final matrix element vs reference: same convention issue as c_E
	// (sign may be flipped; see Task 3 discussion). Do not fix here.
	//
	// / [CN] 1PE-CT 项 (c_D)：真实分波矩阵元，包含 rank-0 和 rank-2 分解
	// 以及24点 Gauss-Legendre x 积分。选择定则：L_2N=0（接触配对顶点）。
	double W1_1pe_contact(int alpha_r, int alpha_c,
						  double p_r, double q_r,
						  double p_c, double q_c,
						  const pw_3N_statespace& pw_states) const
	{
		if (m_c_D == 0.0) return 0.0;

		// Rank-0 recoupling: (1/3)(σ₁·σ₃)(τ₁·τ₃) in 3N Jj basis.
		// Selection rules: L_2N=L_2N'=0, L_1N=L_1N' (enforced inside helper).
		const double recoup0 = recoupling_3nf_1pe_ct_scalar(
			pw_states.L_2N_array[alpha_r], pw_states.S_2N_array[alpha_r],
			pw_states.J_2N_array[alpha_r], pw_states.T_2N_array[alpha_r],
			pw_states.L_1N_array[alpha_r], pw_states.two_J_1N_array[alpha_r],
			pw_states.two_J_3N_array[alpha_r],
			pw_states.L_2N_array[alpha_c], pw_states.S_2N_array[alpha_c],
			pw_states.J_2N_array[alpha_c], pw_states.T_2N_array[alpha_c],
			pw_states.L_1N_array[alpha_c], pw_states.two_J_1N_array[alpha_c],
			pw_states.two_T_3N_array[alpha_r]);

		// Rank-2 recoupling: [σ₁⊗σ₃]₂·Y₂(q̂) coefficient.
		// For l_1N=0 channels (dominant triton configs): returns 0 by CG selection rule.
		const double recoup2 = recoupling_3nf_1pe_ct_rank2(
			pw_states.L_2N_array[alpha_r], pw_states.S_2N_array[alpha_r],
			pw_states.J_2N_array[alpha_r], pw_states.T_2N_array[alpha_r],
			pw_states.L_1N_array[alpha_r], pw_states.two_J_1N_array[alpha_r],
			pw_states.two_J_3N_array[alpha_r],
			pw_states.L_2N_array[alpha_c], pw_states.S_2N_array[alpha_c],
			pw_states.J_2N_array[alpha_c], pw_states.T_2N_array[alpha_c],
			pw_states.L_1N_array[alpha_c], pw_states.two_J_1N_array[alpha_c],
			pw_states.two_T_3N_array[alpha_r]);

		if (recoup0 == 0.0 && recoup2 == 0.0) return 0.0;

		// Squared-Gaussian regulator per E2002 eq. (3.19).
		const double f_bra = chiral_3nf::regulator_gauss(p_r, q_r, m_Lambda);
		const double f_ket = chiral_3nf::regulator_gauss(p_c, q_c, m_Lambda);

		// Gauss-Legendre x-integration for the 1PE pion propagator.
		// Q²(x) = q² + q'² − 2qq'x  (|Δq|²  per E2002 eq. A-2)
		// rank-0 weight: P_0(x) = 1
		// rank-2 weight: P_2(x) = (3x²−1)/2
		double integ0 = 0.0, integ2 = 0.0;
		for (int ix = 0; ix < N_GL; ++ix) {
			const double x  = m_gl_x[ix];
			const double wx = m_gl_w[ix];
			const double k  = chiral_3nf::kernel_1pe_contact(p_r, q_r, p_c, q_c,
			                                                  x, m_mpi_fm);
			if (recoup0 != 0.0) integ0 += wx * k;               // P_0 = 1
			if (recoup2 != 0.0) integ2 += wx * k * (1.5*x*x - 0.5); // P_2(x)
		}

		// Overall coefficient: −g_A c_D / (8 f_π⁴ Λ_χ) × 2
		// (factor of 2 from summing j=2 and j=3 in the operator).
		const double coeff = -m_gA * m_c_D / (8.0 * m_fpi4 * m_Lambda_chi) * 2.0;

		return coeff * f_bra * f_ket * (recoup0 * integ0 + recoup2 * integ2);
	}

	// [EN] 2PE term (c₁, c₃): two-pion exchange between particles 2,3 (pair endpoints)
	// through particle 1 (spectator/center). Implements the scalar (rank-0) channel
	// with monopole (angle-averaged) pion propagators.
	// The tensor (rank-2) part and c₄ cross-product term are deferred.
	//
	// Scalar channel of (σ₂·q₂)(σ₃·q₃) = (1/3)(σ₂·σ₃)(q₂·q₃)
	// Since σ₂·σ₃ and τ₂·τ₃ are pair eigenvalues, this contribution is diagonal in alpha.
	//
	// / [CN] 2PE 项 (c₁, c₃)：粒子2和3（配对端点）通过粒子1（旁观者/中心）的双π交换。
	// 此处实现标量（rank-0）通道，使用单极（角度平均）π传播子。
	// 张量（rank-2）部分和 c₄ 叉积项推迟实现。
	double W1_2pe(int alpha_r, int alpha_c,
				  double p_r, double q_r,
				  double p_c, double q_c,
				  const pw_3N_statespace& pw_states) const
	{
		if (m_c1 == 0.0 && m_c3 == 0.0) return 0.0;

		// Scalar channel: diagonal in all quantum numbers
		if (alpha_r != alpha_c) return 0.0;

		int S_2N = pw_states.S_2N_array[alpha_r];
		int T_2N = pw_states.T_2N_array[alpha_r];

		// σ₂·σ₃ eigenvalue: 2S(S+1) - 3
		double sig23 = 2.0 * S_2N * (S_2N + 1.0) - 3.0;

		// τ₂·τ₃ eigenvalue: 2T(T+1) - 3
		double tau23 = 2.0 * T_2N * (T_2N + 1.0) - 3.0;

		// Angle-averaged momentum transfer quantities (fm⁻²)
		double A_sq = p_r*p_r + p_c*p_c;           // <|Δp|²>
		double B_sq = 0.25*(q_r*q_r + q_c*q_c);    // <|Δq/2|²>
		double Q_sq = A_sq + B_sq;                   // <q₂²> = <q₃²>
		double q2_dot_q3 = B_sq - A_sq;             // <q₂·q₃>

		// Propagator product: 1/((Q² + m_π²)²)
		double denom = Q_sq + m_mpi_fm * m_mpi_fm;
		double prop = 1.0 / (denom * denom);

		// c₁ contribution: (-4c₁m²) × <q₂·q₃>
		double c1_part = -4.0 * m_c1 * m_mpi_fm * m_mpi_fm * q2_dot_q3;

		// c₃ contribution: (2c₃) × <q₂·q₃>²
		double c3_part = 2.0 * m_c3 * q2_dot_q3 * q2_dot_q3;

		// Scalar factor: (1/3) × σ₂·σ₃ × τ₂·τ₃ × (c₁_part + c₃_part) × propagator
		double scalar = (1.0/3.0) * sig23 * tau23 * (c1_part + c3_part) * prop;

		// Overall coefficient: g_A² / (4 f_π⁴)
		double coeff = m_gA * m_gA / (4.0 * m_fpi4);

		// Gaussian regulator per E2002 eq. (3.19):
		//   f_R(p,q) = exp(-((4p² + 3q²)/(4Λ²))²) = exp(-((p² + ¾q²)/Λ²)²)
		double inv_L2 = 1.0 / (m_Lambda * m_Lambda);
		double a_bra  = (p_r*p_r + 0.75*q_r*q_r) * inv_L2;
		double a_ket  = (p_c*p_c + 0.75*q_c*q_c) * inv_L2;
		double f_bra = std::exp(-a_bra * a_bra);
		double f_ket = std::exp(-a_ket * a_ket);

		return coeff * scalar * f_bra * f_ket;
	}
};

#endif // CHIRAL_N2LO_3NF_H
