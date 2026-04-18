#ifndef CHIRAL_N2LO_3NF_H
#define CHIRAL_N2LO_3NF_H

#include "three_nucleon_force_model.h"
#include "constants.h"
#include "spin_isospin_algebra.h"
#include <cmath>

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
	chiral_N2LO_3NF(double c_D, double c_E, double Lambda_3NF_MeV,
				    double c1 = 0.0, double c3 = 0.0, double c4 = 0.0)
		: m_c_D(c_D)
		, m_c_E(c_E)
		, m_Lambda(Lambda_3NF_MeV / hbarc)           // Convert MeV → fm⁻¹
		, m_Lambda_chi(700.0 / hbarc)                 // Λ_χ = 700 MeV → fm⁻¹
		, m_fpi4_inv((fpi/hbarc)*(fpi/hbarc)*(fpi/hbarc)*(fpi/hbarc))  // (f_π/ħc)⁴ in fm⁻⁴
		, m_gA(gA)                                    // Axial coupling constant (dimensionless)
		, m_mpi_fm(mpi / hbarc)                       // Pion mass in fm⁻¹
		, m_c1(c1)                                    // c₁ LEC (GeV⁻¹, from 2NF)
		, m_c3(c3)                                    // c₃ LEC (GeV⁻¹, from 2NF)
		, m_c4(c4)                                    // c₄ LEC (GeV⁻¹, from 2NF) — reserved for tensor part
	{
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
	double m_fpi4_inv;     // (f_π/ħc)⁴ in fm⁻⁴
	double m_gA;           // Axial coupling constant (dimensionless)
	double m_mpi_fm;       // Pion mass in fm⁻¹
	double m_c1;           // c₁ LEC (GeV⁻¹, from 2NF)
	double m_c3;           // c₃ LEC (GeV⁻¹, from 2NF)
	double m_c4;           // c₄ LEC (GeV⁻¹, from 2NF) — reserved for tensor part

	// [EN] 3N contact term (c_E): diagonal in all quantum numbers α.
	//   W^(1)_CT(α',p',q'; α,p,q) = c_E/(f_π⁴ Λ_χ) × [2T_2N(T_2N+1) - 3]
	//                                × f_Λ(p',q') × f_Λ(p,q) × δ_{α',α}
	//
	// where:
	//   τ₂·τ₃ = 2T_2N(T_2N+1) - 3  (eigenvalue of pair isospin operator)
	//   f_Λ(p,q) = exp(-(p² + ¾q²)/Λ²)  (Gaussian regulator)
	//
	// / [CN] 3N 接触项 (c_E)：在所有量子数 α 上对角。
	double W1_contact(int alpha_r, int alpha_c,
					  double p_r, double q_r,
					  double p_c, double q_c,
					  const pw_3N_statespace& pw_states) const
	{
		// Contact term is diagonal in all quantum numbers
		if (alpha_r != alpha_c) return 0.0;

		// Skip if c_E is zero (pure 2NF mode with 3NF object still constructed)
		if (m_c_E == 0.0) return 0.0;

		// Pair isospin for this partial-wave state
		int T_2N = pw_states.T_2N_array[alpha_r];

		// τ₂·τ₃ eigenvalue: 2T(T+1) - 3
		//   T_2N=0 (singlet): -3
		//   T_2N=1 (triplet): +1
		double tau_dot_tau = 2.0 * T_2N * (T_2N + 1.0) - 3.0;

		// Gaussian regulator: f_Λ(p, q) = exp(-(p² + ¾q²) / Λ²)
		double inv_L2 = 1.0 / (m_Lambda * m_Lambda);
		double f_bra = std::exp(-(p_r*p_r + 0.75*q_r*q_r) * inv_L2);
		double f_ket = std::exp(-(p_c*p_c + 0.75*q_c*q_c) * inv_L2);

		// Overall coefficient: c_E / (f_π⁴ × Λ_χ) in fm⁵ (all in fm⁻¹ units)
		double coeff = m_c_E / (m_fpi4_inv * m_Lambda_chi);

		return coeff * tau_dot_tau * f_bra * f_ket;
	}

	// [EN] 1PE-CT term (c_D): one-pion exchange between spectator (particle 1) and a pair
	// particle (particle 3), with a contact interaction in the pair (2,3). This implements
	// only the scalar (rank-0) part proportional to σ₁·σ₃ with a monopole (angle-averaged)
	// pion propagator. The tensor (rank-2) part will be added in a follow-up task.
	//
	// W^(1)_D,scalar = -c_D gA / (8 f_π⁴ Λ_χ) × 2
	//   × f_Λ(p',q') f_Λ(p,q) × τ₁·τ₃ × σ₁·σ₃ × (1/3) × <q₃²/(q₃² + m_π²)>
	//
	// The factor of 2 accounts for pair symmetry (summing j=2,3 contributions).
	// The monopole pion propagator uses <q₃²> = p² + p'² + (q² + q'²)/4 which arises
	// from angle-averaging q₃ = (p⃗ - p⃗') + (q⃗' - q⃗)/2 with uncorrelated directions.
	//
	// / [CN] 1PE-CT 项 (c_D)：旁观者(粒子1)与配对粒子(粒子3)之间的单 π 交换，配对(2,3)间
	// 为接触相互作用。此处仅实现正比于 σ₁·σ₃ 的标量(rank-0)部分，使用单极(角度平均)
	// π 传播子。张量(rank-2)部分将在后续任务中添加。
	double W1_1pe_contact(int alpha_r, int alpha_c,
						  double p_r, double q_r,
						  double p_c, double q_c,
						  const pw_3N_statespace& pw_states) const
	{
		if (m_c_D == 0.0) return 0.0;

		// Extract quantum numbers
		int L_r = pw_states.L_2N_array[alpha_r];  int L_c = pw_states.L_2N_array[alpha_c];
		int S_r = pw_states.S_2N_array[alpha_r];  int S_c = pw_states.S_2N_array[alpha_c];
		int J_r = pw_states.J_2N_array[alpha_r];  int J_c = pw_states.J_2N_array[alpha_c];
		int T_r = pw_states.T_2N_array[alpha_r];  int T_c = pw_states.T_2N_array[alpha_c];
		int l_r = pw_states.L_1N_array[alpha_r];  int l_c = pw_states.L_1N_array[alpha_c];
		int tj_r = pw_states.two_J_1N_array[alpha_r]; int tj_c = pw_states.two_J_1N_array[alpha_c];
		int tJ3 = pw_states.two_J_3N_array[alpha_r]; // same for both (conserved)

		// Selection rules for σ₁·σ₃ scalar part:
		// Orbital angular momenta conserved: L'=L, l'=l
		if (L_r != L_c || l_r != l_c) return 0.0;

		// Compute spin recoupling: σ₁·σ₃
		double sig13 = reduced_me_sigma1_dot_sigma3(
			L_r, S_r, J_r, l_r, tj_r, tJ3,
			L_c, S_c, J_c, l_c, tj_c);
		if (std::abs(sig13) < 1e-15) return 0.0;

		// Compute isospin recoupling: τ₁·τ₃
		int tT3 = pw_states.two_T_3N_array[alpha_r];
		double tau13 = tau1_dot_tau3(T_r, T_c, tT3);
		if (std::abs(tau13) < 1e-15) return 0.0;

		// Gaussian regulator: f_Λ(p,q) = exp(-(p² + ¾q²)/Λ²)
		double inv_L2 = 1.0 / (m_Lambda * m_Lambda);
		double f_bra = std::exp(-(p_r*p_r + 0.75*q_r*q_r) * inv_L2);
		double f_ket = std::exp(-(p_c*p_c + 0.75*q_c*q_c) * inv_L2);

		// Angle-averaged pion propagator (scalar part):
		// q₃ = (p⃗ - p⃗') + (q⃗' - q⃗)/2, so q₃² = |Δp|² + |Δq|²/4 + (Δp)·(Δq)
		// For the monopole approximation with uncorrelated angular directions:
		//   <|Δp|²> = p² + p'²,  <|Δq|²> = q² + q'²,  <Δp·Δq> = 0
		// Hence <q₃²> = p² + p'² + (q² + q'²)/4
		// The scalar channel factor is (1/3) × <q₃²/(q₃² + m_π²)>
		// ≈ (1/3) × <q₃²> / (<q₃²> + m_π²)
		double q3_sq_avg = p_r*p_r + p_c*p_c + 0.25*(q_r*q_r + q_c*q_c);
		double pion_factor = q3_sq_avg / (q3_sq_avg + m_mpi_fm * m_mpi_fm) / 3.0;

		// Overall coefficient: -c_D * gA / (8 * f_π⁴ * Λ_χ) × 2 (pair symmetry)
		double coeff = -m_c_D * m_gA / (8.0 * m_fpi4_inv * m_Lambda_chi) * 2.0;

		return coeff * f_bra * f_ket * tau13 * sig13 * pion_factor;
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
		double coeff = m_gA * m_gA / (4.0 * m_fpi4_inv);

		// Gaussian regulator
		double inv_L2 = 1.0 / (m_Lambda * m_Lambda);
		double f_bra = std::exp(-(p_r*p_r + 0.75*q_r*q_r) * inv_L2);
		double f_ket = std::exp(-(p_c*p_c + 0.75*q_c*q_c) * inv_L2);

		return coeff * scalar * f_bra * f_ket;
	}
};

#endif // CHIRAL_N2LO_3NF_H
