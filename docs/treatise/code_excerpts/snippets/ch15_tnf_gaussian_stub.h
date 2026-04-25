// ===============================================================
// 抽取自仓库 [current]: src/interactions/three_nucleon_force_gaussian_stub.h
// 行号区段：14..41
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
