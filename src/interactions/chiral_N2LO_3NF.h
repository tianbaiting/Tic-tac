#ifndef CHIRAL_N2LO_3NF_H
#define CHIRAL_N2LO_3NF_H

#include "three_nucleon_force_model.h"
#include "constants.h"
#include "spin_isospin_algebra.h"
#include "chiral_3nf_pw_kernels.h"
#include "../utils/chiral_3nf_recoupling.h"
#include "../utils/error_management.h"
#include "gauss_legendre.h"
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

// [EN] Chiral N2LO three-nucleon force: the leading 3NF in chiral effective field theory.
// Three contributions at this order:
//   1. Two-pion exchange (2PE):   proportional to c₁, c₃ (inherited from 2NF)
//   2. One-pion exchange contact (1PE-CT): proportional to c_D
//   3. Three-nucleon contact (CT): proportional to c_E
//
// **IMPLEMENTED TERMS**: c_E, c_D, c₁, c₃ ONLY. The c₄ term
// (Epelbaum 2002 eq. 2.2-2.3, isospin τ₁·(τ₂×τ₃), momentum cross-product
// σ·(qᵢ×qⱼ)) is **NOT implemented**. The model is therefore an APPROXIMATION
// to the full N²LO 3NF and is honestly named `chiral_N2LO_c1c3cDcE_approx`
// (see docs/three_nf_equation_contract.md §8).
//
// **FAIL-CLOSED CONTRACT** (fix/3nf-physics-contract Phase 1):
//   * The constructor REJECTS c₄ ≠ 0 by raising a fatal error. The c₄ term
//     must NEVER be silently dropped — passing a non-zero c₄ means the user
//     expects a full N²LO 3NF, which this model is not.
//   * The factory `three_nucleon_force_model::fetch` rejects the string
//     `chiral_N2LO` (which would claim a complete N²LO 3NF) and only accepts
//     `chiral_N2LO_c1c3cDcE_approx` (plus the legacy alias `chiral_N2LO_without_c4`).
//   * `c4_implemented()` always returns false.
//
// W^(1) is the spectator-1 decomposition: particle 1 is the spectator, particles 2 and 3
// form the interacting pair. The full 3NF is recovered via W = (1+P+P²) W^(1) where
// P = P₁₂₃ + P₁₃₂. The AGS KERNEL uses W^(1)·(1+P) (left W^(1), right (1+P)) —
// see three_nucleon_force_model.h and tests/cpp/test_faddeev_operator_order.cpp.
//
// References:
//   Epelbaum, Glöckle, Meissner, PRC 66 (2002) 064001
//   Witała et al., PRC 77 (2008) 034004
//   Hebeler et al., PRC 91 (2015) 044001
//
// / [CN] 手征 N2LO 三体核力。已实现：c_E, c_D, c₁, c₃（近似模型，名为
// chiral_N2LO_c1c3cDcE_approx）。c₄ 项未实现——构造时若 c₄≠0 将抛出致命错误，
// 绝不静默丢弃。
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
		, m_c4(c4 * hbarc / 1000.0)                   // c₄ LEC: GeV⁻¹ → fm — NOT IMPLEMENTED
		, m_c4_input(c4)                              // raw c₄ value (GeV⁻¹) for the error message
		, m_gl_x(N_GL)
		, m_gl_w(N_GL)
	{
		// [EN] Pre-compute Gauss-Legendre nodes and weights on [-1, +1].
		gauss(m_gl_x.data(), m_gl_w.data(), N_GL);

		// [EN] FAIL-CLOSED on c₄ (Phase 1, fix/3nf-physics-contract). The c₄
		// 2PE cross-product term (Epelbaum 2002 eq. 2.2-2.3,
		// τ₁·(τ₂×τ₃) × σ·(qᵢ×qⱼ)) is NOT implemented in this model. A non-zero
		// c₄ is REJECTED here — it must never be silently dropped, because the
		// user passing c₄ ≠ 0 expects a complete N²LO 3NF that this model is
		// not. The full c₄ PWD requires the 5D angular integral (Golak 2010)
		// with the τ·(τ×τ) isospin recoupling and the (σ×σ)·(q×q) momentum
		// cross-product; see docs/three_nf_equation_contract.md §8 and
		// docs/3nf_audit_2026-06-21.md §B2 for the implementation plan.
		// / [CN] c₄ 硬阻断：c₄ 项未实现，c₄≠0 时抛出致命错误，绝不静默丢弃。
		if (m_c4 != 0.0) {
			std::fprintf(stderr,
				"*** chiral_N2LO_3NF FATAL: c4 = %.6f GeV^-1 was supplied but the c_4 term "
				"(Epelbaum 2002 eq. 2.2-2.3, tau_1.(tau_2 x tau_3) x sigma.(q_i x q_j)) "
				"is NOT IMPLEMENTED in this model (chiral_N2LO_c1c3cDcE_approx). "
				"Refusing to silently drop c_4 — this would produce a wrong physics result. "
				"Either set c4=0 (and use only c_E, c_D, c_1, c_3) or implement the c_4 PWD. "
				"See docs/three_nf_equation_contract.md §8.\n", c4);
			raise_error("chiral_N2LO_3NF constructed with non-zero c_4 (not implemented). "
			            "Set c4=0 for the c1/c3/cD/cE approximation; see "
			            "docs/three_nf_equation_contract.md §8.");
		}
	}

	bool enabled() const override {
		// c_4 is rejected at construction, so it cannot enable the model.
		return m_c_D != 0.0 || m_c_E != 0.0 || m_c1 != 0.0 || m_c3 != 0.0;
	}

	std::string name() const override {
		// Honest name: this is the c1/c3/cD/cE approximation, NOT the full N²LO 3NF.
		return "chiral_N2LO_c1c3cDcE_approx";
	}

	// True iff every N²LO term is implemented. Currently always false because
	// c_4 is not implemented (and rejected at construction).
	virtual bool c4_implemented() const { return false; }

	// LEC accessors (GeV⁻¹ units, matching run_params convention) for cache key.
	// See three_nucleon_force_model.h for the rationale (3NF audit B5).
	double lec_c1_gev() const override { return m_c1 * 1000.0 / hbarc; }
	double lec_c3_gev() const override { return m_c3 * 1000.0 / hbarc; }
	double lec_c4_gev() const override { return m_c4_input; }  // raw input value (GeV⁻¹)

	// Status string for run-metadata output.
	virtual std::string capabilities() const {
		return std::string("c_E=implemented, c_D=implemented, ")
		     + "c_1=implemented, c_3=implemented, "
		     + "c_4=" + (c4_implemented() ? "implemented" : "NOT implemented (blocked)");
	}

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
	double m_c4;           // c₄ LEC in fm (converted from GeV⁻¹) — NOT IMPLEMENTED, hard-blocked
	double m_c4_input;     // Raw c₄ value (GeV⁻¹) as supplied by the user, for metadata
	bool   m_c4_blocked;   // True iff c₄ ≠ 0 was supplied (triggers name change + warning)
	std::vector<double> m_gl_x; // Gauss-Legendre nodes on [-1, +1]
	std::vector<double> m_gl_w; // Gauss-Legendre weights on [-1, +1]

	// [EN] 3N contact term (c_E): partial-wave matrix element assembled from
	// the factored pieces
	//   W^(1)_CT(α',p',q'; α,p,q)
	//     = A_cE(α', α)                          [recoupling_3nf_contact_cE]
	//     × [ -½ c_E / (f_π⁴ Λ_χ) ]              [kernel_contact]
	//     × f_R(p',q') f_R(p,q)                  [regulator_gauss, squared Gaussian]
	//
	// Per Epelbaum 2002 eq. (2.10)/(A-4), the c_E contact is a PURE SPIN SCALAR.
	// Its spin-isospin recoupling A_cE depends only on the pair isospin T_2N and
	// the spectator 1/2 ↔ T_3N coupling — NO σ_2·σ_3 factor. The previous
	// implementation reused recoupling_3nf_scalar (now recoupling_3nf_2pe_scalar)
	// which erroneously multiplies by σ_2·σ_3 = -3 for S_2N=0 or +1 for S_2N=1.
	// See docs/3nf_audit_2026-06-21.md §B1 for the bug analysis.
	//
	// For 3S1 (S_2N=1, T_2N=0): A_cE = (+1)·(-1)·(-3)·(-1/2) = -3/2 (NEW, correct)
	//                            old: (+1)·(-3) = -3            (BUG: extra σ·σ)
	// For 1S0 (S_2N=0, T_2N=1): A_cE = (-1)·(+1)·(+1)·(+1/4) = -1/4 (NEW, correct)
	//                            old: (-3)·(+1) = -3              (BUG: extra σ·σ × wrong 6j)
	//
	// References:
	//   Epelbaum et al. PRC 66 (2002) 064001, eqs. (2.10), (3.19), (A-4);
	//   tools/check_3nf_normalization/formula_reference.md §1.
	//
	// / [CN] 3N 接触项 (c_E)：由重耦合系数（Epelbaum A-4，纯 τ·τ）、LEC 核
	// (-½ c_E/(f_π⁴ Λ_χ))、E2002 eq. 3.19 平方高斯正规化子组合而成。
	// c_E 是纯自旋标量——不依赖 σ_2·σ_3（修复 B1）。
	double W1_contact(int alpha_r, int alpha_c,
					  double p_r, double q_r,
					  double p_c, double q_c,
					  const pw_3N_statespace& pw_states) const
	{
		// Short-circuit when the LEC is off (pure 2NF runs still build this object).
		if (m_c_E == 0.0) return 0.0;

		// c_E contact recoupling (Epelbaum A-4): pure τ_2·τ_3 spin-scalar.
		// NO σ_2·σ_3 factor — see docs/3nf_audit_2026-06-21.md §B1.
		const double recoup = recoupling_3nf_contact_cE(
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

	// [EN] 2PE term (c₁, c₃): two-pion exchange between pair particles 2,3 via
	// spectator 1. Implements the TRUE partial-wave matrix element via rank-0 + rank-2
	// decomposition of (σ₂·q₂)(σ₃·q₃):
	//
	//   (σ₂·q₂)(σ₃·q₃) = ⅓(σ₂·σ₃)(q₂·q₃)         [rank-0, scalar]
	//                   + [σ₂⊗σ₃]₂·[q₂⊗q₃]₂         [rank-2, tensor]
	//
	// The (τ₂·τ₃) isospin factor is pair-diagonal (eigenvalue) for both ranks.
	//
	// Operator (from [G2010] eq. 18, [E2002] eq. 2.2-2.3):
	//   V^(1)_2π = (gA/2fπ)² × (τ₂·τ₃) × (σ₂·q₂)(σ₃·q₃) × F₁(q₂,q₃)
	//   F₁ = [-4c₁mπ²/fπ² + 2c₃(q₂·q₃)/fπ²] / [(q₂²+mπ²)(q₃²+mπ²)]
	//
	// Rank-0 recoupling: recoupling_3nf_2pe_scalar returns (σ₂·σ₃)(τ₂·τ₃) pair
	// eigenvalues with full channel-diagonal selection rules. The ⅓ rank-0
	// coefficient is applied explicitly in the outer coefficient.
	//
	// Rank-2 recoupling: recoupling_3nf_rank2 returns τ₂·τ₃ × √30 × hat × 9j × CG
	// (pair-only operator, opens L_2N=0 ↔ L_2N=2 transitions). Combined with
	// kernel_2pe_c1c3_rank2 (which supplies the spatial factor 2/√30 × pp²/prop²),
	// the total rank-2 contribution is τ₂·τ₃ × 2 × hat × 9j × CG × lec × pp²/prop².
	//
	// Key change from legacy: the α_r != α_c early-exit is REMOVED because the rank-2
	// tensor couples L_2N=0 ↔ L_2N=2 (e.g. 3S1 ↔ 3D1 in the pair), introducing
	// off-diagonal elements that were previously dropped. These off-diagonal matrix
	// elements are the primary cause of the sign flip in the 3NF normalization check.
	//
	// Angular integration: 24-point Gauss-Legendre quadrature over x=cos(q̂·q̂')∈[-1,+1].
	//   rank-0: ∫dx kernel_2pe_c1c3(p_r,q_r,p_c,q_c,x,...)  [1/(8π³) included in kernel]
	//   rank-2: ∫dx kernel_2pe_c1c3_rank2(p_r,q_r,p_c,q_c,x,...) [spatial p_c² factor]
	//
	// Overall prefactor: (gA/2fπ)² = gA²/(4fπ²). Note: fπ in fm⁻¹, so fπ² = fm⁻²
	// and (gA/2fπ)² has units fm². The fπ² is NOT fπ⁴ — see kernel_2pe_c1c3 which
	// already includes one fπ² factor inside lec_bracket.
	//
	// Note on c_E / c_D sign conventions: same deferred issue applies here — do NOT
	// flip signs to make X positive; the sign conventions for c_1 and c_3 follow
	// [G2010] eq. (18) with no additional adjustments.
	//
	// / [CN] 2PE 项 (c₁, c₃)：真实分波矩阵元，包含 rank-0 和 rank-2 分解以及
	// 24点 Gauss-Legendre x 积分。关键：删除 alpha_r != alpha_c 早退，因为
	// rank-2 张量耦合 L_2N=0 ↔ L_2N=2（如 3S1↔3D1），引入非对角矩阵元。
	double W1_2pe(int alpha_r, int alpha_c,
				  double p_r, double q_r,
				  double p_c, double q_c,
				  const pw_3N_statespace& pw_states) const
	{
		if (m_c1 == 0.0 && m_c3 == 0.0) return 0.0;

		// [EN] Rank-0 recoupling: (σ₂·σ₃)(τ₂·τ₃) pair eigenvalues, diagonal in all
		// pair and spectator quantum numbers. Returns 0 for off-diagonal α pairs.
		// Note: recoupling_3nf_2pe_scalar returns sigma*sigma × tau*tau WITHOUT the ⅓
		// rank-0 coefficient — we apply 1/3 explicitly in the overall coefficient below.
		const double recoup0 = recoupling_3nf_2pe_scalar(
			pw_states.L_2N_array[alpha_r], pw_states.S_2N_array[alpha_r],
			pw_states.J_2N_array[alpha_r], pw_states.T_2N_array[alpha_r],
			pw_states.L_1N_array[alpha_r], pw_states.two_J_1N_array[alpha_r],
			pw_states.two_J_3N_array[alpha_r],
			pw_states.L_2N_array[alpha_c], pw_states.S_2N_array[alpha_c],
			pw_states.J_2N_array[alpha_c], pw_states.T_2N_array[alpha_c],
			pw_states.L_1N_array[alpha_c], pw_states.two_J_1N_array[alpha_c],
			pw_states.two_T_3N_array[alpha_r]);

		// [EN] Rank-2 recoupling: (τ₂·τ₃) × √30 × hat × 9j × CG(L_c,0;2,0|L_r,0).
		// Non-zero for L_r - L_c = ±2 with S_r=S_c=1 (triplet pairs only).
		// Enables 3S1 ↔ 3D1 coupling and the associated sign flip.
		const double recoup2 = recoupling_3nf_rank2(
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

		// fπ in fm⁻¹ (needed for (gA/2fπ)² = gA²/(4fπ²) prefactor).
		// Computed from m_fpi4 = fπ⁴ stored in the class.
		const double fpi_fm = std::sqrt(std::sqrt(m_fpi4));

		// [EN] Rank-2 Hermitian symmetrization: compute the transposed recoupling
		// recoup2_T = recoupling with bra and ket swapped (alpha_c→alpha_r convention).
		// The spatial kernel kernel_2pe_c1c3_rank2 has a pp² = p_c² dependence that
		// makes W1(alpha_r,alpha_c,p_r,q_r,p_c,q_c) ≠ W1(alpha_c,alpha_r,p_c,q_c,p_r,q_r)
		// when only one side is evaluated. Explicit symmetrization:
		//   rank2_contribution = ½ [recoup2(r,c) × integ2(p_r,q_r,p_c,q_c)
		//                          + recoup2(c,r) × integ2(p_c,q_c,p_r,q_r)]
		// ensures the full W1 satisfies W1(r,c,p_r,q_r,p_c,q_c) = W1(c,r,p_c,q_c,p_r,q_r).
		const double recoup2_T = recoupling_3nf_rank2(
			pw_states.L_2N_array[alpha_c], pw_states.S_2N_array[alpha_c],
			pw_states.J_2N_array[alpha_c], pw_states.T_2N_array[alpha_c],
			pw_states.L_1N_array[alpha_c], pw_states.two_J_1N_array[alpha_c],
			pw_states.two_J_3N_array[alpha_c],
			pw_states.L_2N_array[alpha_r], pw_states.S_2N_array[alpha_r],
			pw_states.J_2N_array[alpha_r], pw_states.T_2N_array[alpha_r],
			pw_states.L_1N_array[alpha_r], pw_states.two_J_1N_array[alpha_r],
			pw_states.two_T_3N_array[alpha_c]);

		const bool need_rank2 = (recoup2 != 0.0 || recoup2_T != 0.0);

		// Gauss-Legendre x-quadrature for pion propagator kernels.
		// x = cos(q̂·q̂') ∈ [-1, +1]; Q²(x) = q² + q'² - 2qq'x = |Δq|²(x).
		double integ0 = 0.0, integ2 = 0.0, integ2_T = 0.0;
		for (int ix = 0; ix < N_GL; ++ix) {
			const double x  = m_gl_x[ix];
			const double wx = m_gl_w[ix];
			if (recoup0 != 0.0) {
				integ0 += wx * chiral_3nf::kernel_2pe_c1c3(
					p_r, q_r, p_c, q_c, x, m_mpi_fm, m_c1, m_c3, fpi_fm);
			}
			if (need_rank2) {
				// Direct: kernel with pp = p_c (ket pair momentum)
				const double k2 = chiral_3nf::kernel_2pe_c1c3_rank2(
					p_r, q_r, p_c, q_c, x, m_mpi_fm, m_c1, m_c3, fpi_fm);
				// Transposed: swap (p_r,q_r) ↔ (p_c,q_c) so pp = p_r (bra pair momentum)
				const double k2_T = chiral_3nf::kernel_2pe_c1c3_rank2(
					p_c, q_c, p_r, q_r, x, m_mpi_fm, m_c1, m_c3, fpi_fm);
				if (recoup2   != 0.0) integ2   += wx * k2;
				if (recoup2_T != 0.0) integ2_T += wx * k2_T;
			}
		}

		// Overall coefficient: (gA/2fπ)² = gA²/(4fπ²).
		// The fπ² here matches [G2010] eq. (18): (gA/2fπ)² prefactor in F₁.
		// Note: this is fπ² NOT fπ⁴ (the kernels already carry one fπ² inside lec_bracket).
		// Units: gA² dimensionless, fπ² in fm⁻², so coeff has units fm².
		const double coeff = m_gA * m_gA / (4.0 * fpi_fm * fpi_fm);

		// Rank-0 carries a ⅓ from the (σ₂·q₂)(σ₃·q₃) = ⅓(σ₂·σ₃)(q₂·q₃) decomposition.
		// Rank-2 does NOT carry ⅓ — it is the full [σ₂⊗σ₃]₂·[q₂⊗q₃]₂ piece.
		// The ½ symmetrization factor ensures Hermitian symmetry: each side contributes
		// half of the combined recoup2×integ2 + recoup2_T×integ2_T.
		return coeff * f_bra * f_ket
		       * ((1.0/3.0) * recoup0 * integ0
		          + 0.5 * (recoup2 * integ2 + recoup2_T * integ2_T));
	}
};

#endif // CHIRAL_N2LO_3NF_H
