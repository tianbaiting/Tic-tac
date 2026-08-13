#include "chiral_N2LO_3NF_factorized.h"
#include "chiral_N2LO_3NF_full_reference.h"

#include <algorithm>
#include <cmath>
#include <cstdio>

namespace {

int failures = 0;

void expect_close(const char* label, double value, double expected,
	                double relative_tolerance = 2.0e-10)
{
	const double tolerance = relative_tolerance * std::max(1.0, std::abs(expected));
	if (std::abs(value - expected) > tolerance) {
		std::printf("FAIL %-32s value=% .17e expected=% .17e\n", label, value, expected);
		++failures;
	} else {
		std::printf("PASS %-32s value=% .17e\n", label, value);
	}
}

void expect_true(const char* label, bool condition)
{
	if (!condition) {
		std::printf("FAIL %s\n", label);
		++failures;
	} else {
		std::printf("PASS %s\n", label);
	}
}

struct test_space {
	int l_pair[6] = {1, 0, 0, 2, 1, 0};
	int s_pair[6] = {1, 1, 0, 1, 1, 0};
	int j_pair[6] = {1, 1, 0, 1, 2, 0};
	int t_pair[6] = {1, 0, 1, 0, 1, 1};
	int lambda[6] = {1, 0, 0, 0, 0, 1};
	int two_i[6] = {1, 1, 1, 1, 1, 3};
	int two_j[6] = {1, 1, 1, 1, 3, 3};
	int two_t[6] = {1, 1, 1, 1, 1, 3};
	int parity[6] = {1, 1, 1, 1, -1, -1};
	pw_3N_statespace pw{};

	test_space()
	{
		pw.Nalpha = 6;
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

} // namespace

int main()
{
	test_space space;
	constexpr double p_bra = 0.9;
	constexpr double q_bra = 1.1;
	constexpr double p_ket = 0.4;
	constexpr double q_ket = 0.7;
	constexpr int order = 8;

	// Golden values are from the independently tested Python Hebeler
	// factorization at the same order, cutoff, constants, and momenta.
	chiral_N2LO_3NF_factorized c1(0.0, 0.0, 500.0, -0.81, 0.0, 0.0, order);
	chiral_N2LO_3NF_factorized c3(0.0, 0.0, 500.0, 0.0, -3.2, 0.0, order);
	chiral_N2LO_3NF_factorized c4(0.0, 0.0, 500.0, 0.0, 0.0, 5.4, order);
	chiral_N2LO_3NF_factorized cD(-0.2, 0.0, 500.0, 0.0, 0.0, 0.0, order);
	chiral_N2LO_3NF_factorized cE(0.0, -0.205, 500.0, 0.0, 0.0, 0.0, order);

	expect_close("S diagonal c1", c1.W1_element(
		1, 1, p_bra, q_bra, p_ket, q_ket, space.pw), 1.0573286327355724e-3);
	expect_close("S diagonal c3", c3.W1_element(
		1, 1, p_bra, q_bra, p_ket, q_ket, space.pw), 3.7548841401283677e-3);
	expect_close("S transition c4", c4.W1_element(
		1, 2, p_bra, q_bra, p_ket, q_ket, space.pw), 1.3130109844946473e-2);
	expect_close("S transition cD", cD.W1_element(
		1, 2, p_bra, q_bra, p_ket, q_ket, space.pw), 6.533793487265248e-4);
	expect_close("S diagonal cE", cE.W1_element(
		1, 1, p_bra, q_bra, p_ket, q_ket, space.pw), 8.6321225934281e-3);

	chiral_N2LO_3NF_factorized c1_order6(
		0.0, 0.0, 500.0, -0.81, 0.0, 0.0, 6);
	chiral_N2LO_3NF_factorized c3_order6(
		0.0, 0.0, 500.0, 0.0, -3.2, 0.0, 6);
	expect_close("P diagonal c1", c1_order6.W1_element(
		0, 0, p_bra, q_bra, p_ket, q_ket, space.pw), 2.0387073870411527e-5);
	expect_close("P diagonal c3", c3_order6.W1_element(
		0, 0, p_bra, q_bra, p_ket, q_ket, space.pw), 6.207736537116497e-4);
	expect_close("S-to-pair-D c1", c1_order6.W1_element(
		1, 3, p_bra, q_bra, p_ket, q_ket, space.pw), 1.4441424077447324e-4);
	expect_close("S-to-pair-D c3", c3_order6.W1_element(
		1, 3, p_bra, q_bra, p_ket, q_ket, space.pw), 3.0771522905696383e-3);
	expect_close("pair-D-to-S c1", c1_order6.W1_element(
		3, 1, p_bra, q_bra, p_ket, q_ket, space.pw), 4.503649205515546e-3);
	expect_close("pair-D-to-S c3", c3_order6.W1_element(
		3, 1, p_bra, q_bra, p_ket, q_ket, space.pw), 1.156232962137228e-2);
	expect_close("c1 Hermitian reverse", c1_order6.W1_element(
		1, 3, p_ket, q_ket, p_bra, q_bra, space.pw), 4.503649205515546e-3);
	expect_close("c3 Hermitian reverse", c3_order6.W1_element(
		1, 3, p_ket, q_ket, p_bra, q_bra, space.pw), 1.156232962137228e-2);
	expect_close("J=3/2 negative-parity c1", c1_order6.W1_element(
		4, 4, p_bra, q_bra, p_ket, q_ket, space.pw), -2.979704746377356e-5);
	expect_close("J=3/2 negative-parity c3", c3_order6.W1_element(
		4, 4, p_bra, q_bra, p_ket, q_ket, space.pw), -3.373887409567732e-5);

	chiral_N2LO_3NF_factorized all(
		-0.2, -0.205, 500.0, -0.81, -3.2, 5.4, order);
	const double off_diagonal = all.W1_element(
		1, 2, p_bra, q_bra, p_ket, q_ket, space.pw);
	const double off_diagonal_reverse = all.W1_element(
		2, 1, p_ket, q_ket, p_bra, q_bra, space.pw);
	expect_close("all-term transition", off_diagonal,
		1.3783489193672998e-2);
	expect_close("exact Hermitian reverse", off_diagonal_reverse,
		off_diagonal, 2.0e-12);
	const std::vector<std::pair<int, int>> batch_channels{{1, 2}, {1, 1}, {0, 3}};
	std::vector<double> batch_values;
	all.W1_elements_for_channels(
		batch_channels, p_bra, q_bra, p_ket, q_ket, space.pw, batch_values);
	for (std::size_t index = 0; index < batch_channels.size(); ++index) {
		const auto channel = batch_channels[index];
		const double scalar_value = all.W1_element(
			channel.first, channel.second,
			p_bra, q_bra, p_ket, q_ket, space.pw);
		expect_true("batch/scalar bitwise equality", batch_values[index] == scalar_value);
	}
	chiral_N2LO_3NF_factorized all_order6(
		-0.2, -0.205, 500.0, -0.81, -3.2, 5.4, 6);
	expect_close("T=3/2 negative-parity all", all_order6.W1_element(
		5, 5, p_bra, q_bra, p_ket, q_ket, space.pw),
		1.7269084962373656e-3);

	// The independent direct-Jj five-angle reference converges toward the
	// factorized result without any Hermitian averaging.
	chiral_N2LO_3NF_full_reference reference4(
		-0.2, -0.205, 500.0, -0.81, -3.2, 5.4, 4);
	chiral_N2LO_3NF_full_reference reference6(
		-0.2, -0.205, 500.0, -0.81, -3.2, 5.4, 6);
	const double reference_value4 = reference4.W1_element(
		1, 2, p_bra, q_bra, p_ket, q_ket, space.pw);
	const double reference_value6 = reference6.W1_element(
		1, 2, p_bra, q_bra, p_ket, q_ket, space.pw);
	expect_true("five-angle N=6 closer than N=4",
		std::abs(reference_value6 - off_diagonal)
		< std::abs(reference_value4 - off_diagonal));
	expect_true("five-angle N=6 agrees below 1e-5",
		std::abs(reference_value6 - off_diagonal) < 1.0e-5);

	if (failures != 0) return 1;
	std::printf("All complete factorized N2LO C++ tests passed.\n");
	return 0;
}
