#include "cache_layer.h"
#include "chiral_N2LO_3NF_factorized.h"
#include "constants.h"
#include "w1_pw_cache.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <string>

namespace {

int failures = 0;

void expect_true(const char* label, bool condition)
{
	if (!condition) {
		std::printf("FAIL %s\n", label);
		++failures;
	} else {
		std::printf("PASS %s\n", label);
	}
}

void expect_relative(const char* label, double value, double expected,
	                 double relative_tolerance)
{
	const double scale = std::max(std::abs(expected), 1.0e-300);
	if (std::abs(value - expected) > relative_tolerance * scale) {
		std::printf("FAIL %-34s value=% .17e expected=% .17e\n",
		            label, value, expected);
		++failures;
	} else {
		std::printf("PASS %-34s value=% .17e\n", label, value);
	}
}

class counting_factorized_tnf final : public three_nucleon_force_model
{
public:
	counting_factorized_tnf()
		: m_inner(-0.2, -0.205, 500.0, -0.81, -3.2, 5.4, 6)
	{}

	bool enabled() const override { return m_inner.enabled(); }
	std::string name() const override { return m_inner.name(); }

	double W1_element(int alpha_r, int alpha_c,
	                  double p_r, double q_r, double p_c, double q_c,
	                  const pw_3N_statespace& pw) const override
	{
		m_calls.fetch_add(1, std::memory_order_relaxed);
		return m_inner.W1_element(alpha_r, alpha_c, p_r, q_r, p_c, q_c, pw);
	}

	double lec_c1_gev() const override { return m_inner.lec_c1_gev(); }
	double lec_c3_gev() const override { return m_inner.lec_c3_gev(); }
	double lec_c4_gev() const override { return m_inner.lec_c4_gev(); }
	int angular_order_3nf() const override { return m_inner.angular_order_3nf(); }
	double axial_coupling_3nf() const override { return m_inner.axial_coupling_3nf(); }
	double pion_decay_constant_mev_3nf() const override
	{
		return m_inner.pion_decay_constant_mev_3nf();
	}
	double pion_mass_mev_3nf() const override { return m_inner.pion_mass_mev_3nf(); }
	double chiral_scale_mev_3nf() const override { return m_inner.chiral_scale_mev_3nf(); }
	double hbarc_mev_fm_3nf() const override { return m_inner.hbarc_mev_fm_3nf(); }

	void reset_calls() const { m_calls.store(0, std::memory_order_relaxed); }
	std::size_t calls() const { return m_calls.load(std::memory_order_relaxed); }

private:
	chiral_N2LO_3NF_factorized m_inner;
	mutable std::atomic<std::size_t> m_calls{0};
};

struct test_space
{
	// Minimal subspace of the production J=1/2+ basis: alpha 4 (3S1) and
	// alpha 0 (1S0) from construct_symmetric_pw_states.  Restricting the cache
	// test to these two channels keeps its all-block direct integration fast.
	int l_pair[2] = {0, 0};
	int s_pair[2] = {1, 0};
	int j_pair[2] = {1, 0};
	int t_pair[2] = {0, 1};
	int lambda[2] = {0, 0};
	int two_i[2] = {1, 1};
	int two_j[2] = {1, 1};
	int two_t[2] = {1, 1};
	int parity[2] = {1, 1};
	pw_3N_statespace pw{};

	test_space()
	{
		pw.Nalpha = 2;
		pw.L_2N_array = l_pair;
		pw.S_2N_array = s_pair;
		pw.J_2N_array = j_pair;
		pw.T_2N_array = t_pair;
		pw.L_1N_array = lambda;
		pw.two_J_1N_array = two_i;
		pw.two_J_3N_array = two_j;
		pw.two_T_3N_array = two_t;
		pw.P_3N_array = parity;
	}
};

run_params make_cache_parameters()
{
	run_params parameters{};
	parameters.Np_WP = 1;
	parameters.Nq_WP = 1;
	parameters.J_2N_max = 1;
	parameters.two_J_3N_max = 1;
	parameters.Np_per_WP_W1 = 2;
	parameters.Nq_per_WP_W1 = 2;
	parameters.Nangle_3NF = 6;
	parameters.potential_model = "N2LOopt";
	parameters.three_nucleon_force = "chiral_N2LO_full_factorized";
	parameters.c_D = -0.2;
	parameters.c_E = -0.205;
	parameters.Lambda_3NF = 500.0;
	parameters.chebyshev_s = 100.0;
	parameters.chebyshev_t = 1.0;
	parameters.tensor_force = true;
	parameters.isospin_breaking_1S0 = false;
	return parameters;
}

} // namespace

int main()
{
	counting_factorized_tnf tnf;
	test_space space;
	const pw_3N_statespace& pw = space.pw;
	run_params parameters = make_cache_parameters();

	// A deliberately wide cell in each radial variable, expressed in the MeV
	// convention consumed by W1_PW_cache.  The independent Python integration
	// uses the equivalent [0.2,0.8] fm^-1 bounds.
	const double p_boundaries[2] = {0.2 * hbarc, 0.8 * hbarc};
	const double q_boundaries[2] = {0.2 * hbarc, 0.8 * hbarc};

	// No cache layer has been initialized yet: this is the direct quadrature path.
	W1_PW_cache cache_off;
	cache_off.build(tnf, p_boundaries, 1, q_boundaries, 1, pw, parameters);
	expect_true("cache-off evaluates W1", tnf.calls() > 0);

	char root_template[] = "/tmp/tictac-w1-cache-test-XXXXXX";
	char* cache_root = ::mkdtemp(root_template);
	if (cache_root == nullptr) {
		std::perror("mkdtemp");
		return 1;
	}
	tictac::cache::initialize(cache_root);

	// First initialized build is a miss and writes every allowed channel block.
	tnf.reset_calls();
	W1_PW_cache cache_store;
	cache_store.build(tnf, p_boundaries, 1, q_boundaries, 1, pw, parameters);
	expect_true("cache miss evaluates W1", tnf.calls() > 0);
	expect_true("cache miss stores W1 blocks", tictac::cache::summary().n_w1 > 0);

	// The identical second build must be served entirely from the HDF5 cache.
	tnf.reset_calls();
	W1_PW_cache cache_hit;
	cache_hit.build(tnf, p_boundaries, 1, q_boundaries, 1, pw, parameters);
	expect_true("cache hit makes zero W1 calls", tnf.calls() == 0);

	for (int alpha_r = 0; alpha_r < pw.Nalpha; ++alpha_r) {
		for (int alpha_c = 0; alpha_c < pw.Nalpha; ++alpha_c) {
			const double uncached = cache_off.get(alpha_r, alpha_c, 0, 0, 0, 0);
			const double stored = cache_store.get(alpha_r, alpha_c, 0, 0, 0, 0);
			const double loaded = cache_hit.get(alpha_r, alpha_c, 0, 0, 0, 0);
			expect_true("cache-off/store bitwise equality", uncached == stored);
			expect_true("cache-off/hit bitwise equality", uncached == loaded);
		}
	}

	// Independent Python production-driver quadrature for the complete c4+cD
	// 3S1->1S0 transition at the same radial and transfer orders.
	const double transition = cache_hit.get(0, 1, 0, 0, 0, 0);
	expect_relative("N=2 wide transition Python oracle", transition,
	                3.196424777580e-2, 3.0e-10);

	tictac::cache::shutdown();
	std::filesystem::remove_all(cache_root);
	return failures == 0 ? 0 : 1;
}
