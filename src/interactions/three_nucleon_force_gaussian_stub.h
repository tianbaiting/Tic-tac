#ifndef THREE_NUCLEON_FORCE_GAUSSIAN_STUB_H
#define THREE_NUCLEON_FORCE_GAUSSIAN_STUB_H

#include "three_nucleon_force_model.h"
#include <cmath>

// [EN] Gaussian stub 3NF for data-flow testing. Returns a separable Gaussian form factor
// that is diagonal in alpha (no partial-wave mixing). The strength is deliberately tiny
// so that the Padé/Neumann iteration still converges normally. This model exists solely
// to verify that the 3NF kernel injection plumbing is wired correctly — it has no physics
// content. / [CN] 用于数据流测试的 Gaussian 桩 3NF。返回一个在 alpha 上对角（无分波混合）的
// 可分离 Gaussian 形状因子。强度故意设得极小以确保 Padé/Neumann 迭代仍正常收敛。
// 此模型仅用于验证 3NF 核注入管道正确接线，没有物理内容。
class three_nucleon_force_gaussian_stub : public three_nucleon_force_model
{
public:
	three_nucleon_force_gaussian_stub() = default;

	bool enabled() const override { return true; }
	std::string name() const override { return "gaussian_stub"; }

	// W^(1)(alpha_r, alpha_c, p_r, q_r, p_c, q_c) = V0 * delta(alpha_r, alpha_c)
	//     * exp(-(p_r^2 + p_c^2) / Lambda^2) * exp(-(q_r^2 + q_c^2) / Lambda^2)
	// V0 = 0.001 MeV fm^3  (tiny), Lambda = 2.0 fm^{-1}
	double W1_element(int alpha_r, int alpha_c,
					  double p_r, double q_r,
					  double p_c, double q_c,
					  const pw_3N_statespace& pw_states) const override
	{
		if (alpha_r != alpha_c) return 0.0;

		constexpr double V0     = 0.001;   // MeV fm^3
		constexpr double Lambda = 2.0;     // fm^{-1}
		constexpr double inv_L2 = 1.0 / (Lambda * Lambda);

		return V0 * std::exp(-(p_r*p_r + p_c*p_c) * inv_L2)
				   * std::exp(-(q_r*q_r + q_c*q_c) * inv_L2);
	}
};

#endif // THREE_NUCLEON_FORCE_GAUSSIAN_STUB_H
