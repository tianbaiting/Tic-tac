#include "chiral_N2LO_3NF_factorized.h"

#include "constants.h"
#include "coupling_coefficients.h"
#include "gauss_legendre.h"

#include <gsl/gsl_sf_legendre.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <exception>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using complex = std::complex<double>;
using vector3 = std::array<double, 3>;
using state8 = std::array<complex, 8>;

constexpr double pi_value = 3.141592653589793238462643383279502884;
constexpr complex imaginary_unit{0.0, 1.0};

struct ls_channel {
	int l_pair;
	int s_pair;
	int lambda_spectator;
	int total_L;
	int two_total_S;
	int two_total_J;
	int t_pair;
	int two_total_T;
};

struct jj_channel {
	int l_pair;
	int s_pair;
	int j_pair;
	int lambda_spectator;
	int two_j_spectator;
	int two_total_J;
	int t_pair;
	int two_total_T;
};

std::array<int, 8> channel_key(const ls_channel& channel)
{
	return {channel.l_pair, channel.s_pair, channel.lambda_spectator,
	        channel.total_L, channel.two_total_S, channel.two_total_J,
	        channel.t_pair, channel.two_total_T};
}

std::array<int, 8> channel_key(const jj_channel& channel)
{
	return {channel.l_pair, channel.s_pair, channel.j_pair,
	        channel.lambda_spectator, channel.two_j_spectator,
	        channel.two_total_J, channel.t_pair, channel.two_total_T};
}

jj_channel make_jj_channel(int alpha, const pw_3N_statespace& pw)
{
	if (alpha < 0 || alpha >= pw.Nalpha) {
		throw std::out_of_range("factorized 3NF channel index is out of range");
	}
	return {
		pw.L_2N_array[alpha], pw.S_2N_array[alpha], pw.J_2N_array[alpha],
		pw.L_1N_array[alpha], pw.two_J_1N_array[alpha],
		pw.two_J_3N_array[alpha], pw.T_2N_array[alpha],
		pw.two_T_3N_array[alpha]
	};
}

using ls_expansion_table = std::vector<std::pair<ls_channel, double>>;

ls_expansion_table build_ls_expansion(const jj_channel& channel)
{
	ls_expansion_table result;
	double norm = 0.0;
	for (int total_L = std::abs(channel.l_pair - channel.lambda_spectator);
	     total_L <= channel.l_pair + channel.lambda_spectator; ++total_L) {
		for (int two_total_S = std::abs(2 * channel.s_pair - 1);
		     two_total_S <= 2 * channel.s_pair + 1; two_total_S += 2) {
			if (std::abs(2 * total_L - two_total_S) > channel.two_total_J) continue;
			if (2 * total_L + two_total_S < channel.two_total_J) continue;
			const double coefficient = std::sqrt(
				(2.0 * channel.j_pair + 1.0)
				* (channel.two_j_spectator + 1.0)
				* (2.0 * total_L + 1.0)
				* (two_total_S + 1.0))
				* wigner_9j(
					2 * channel.l_pair, 2 * channel.s_pair, 2 * channel.j_pair,
					2 * channel.lambda_spectator, 1, channel.two_j_spectator,
					2 * total_L, two_total_S, channel.two_total_J);
			if (std::abs(coefficient) < 2.0e-15) continue;
			result.push_back({{
				channel.l_pair, channel.s_pair, channel.lambda_spectator,
				total_L, two_total_S, channel.two_total_J,
				channel.t_pair, channel.two_total_T
			}, coefficient});
			norm += coefficient * coefficient;
		}
	}
	if (std::abs(norm - 1.0) > 2.0e-12) {
		throw std::runtime_error("factorized 3NF Jj-to-LS recoupling is not unitary");
	}
	return result;
}

std::shared_ptr<const ls_expansion_table> get_ls_expansion(const jj_channel& channel)
{
	using key_type = std::array<int, 8>;
	static std::mutex cache_mutex;
	static std::map<key_type, std::shared_ptr<const ls_expansion_table>> cache;
	const key_type key = channel_key(channel);
	{
		std::lock_guard<std::mutex> lock(cache_mutex);
		const auto found = cache.find(key);
		if (found != cache.end()) return found->second;
	}
	auto candidate = std::make_shared<ls_expansion_table>(build_ls_expansion(channel));
	std::lock_guard<std::mutex> lock(cache_mutex);
	return cache.emplace(key, candidate).first->second;
}

int product_index(int two_m1, int two_m2, int two_m3)
{
	auto bit = [](int two_m) {
		if (two_m == 1) return 0;
		if (two_m == -1) return 1;
		throw std::invalid_argument("spin/isospin projection is not +/-1/2");
	};
	return 4 * bit(two_m1) + 2 * bit(two_m2) + bit(two_m3);
}

state8 zero_state()
{
	state8 result{};
	result.fill(complex{0.0, 0.0});
	return result;
}

complex inner_product(const state8& bra, const state8& ket)
{
	complex result{0.0, 0.0};
	for (int index = 0; index < 8; ++index) result += std::conj(bra[index]) * ket[index];
	return result;
}

state8 coupled_three_half_state(int pair_angular_momentum,
	                              int two_total, int two_m_total)
{
	state8 result = zero_state();
	const int two_pair = 2 * pair_angular_momentum;
	for (int two_m_pair = -two_pair; two_m_pair <= two_pair; two_m_pair += 2) {
		for (int two_m1 : {-1, 1}) {
			const double outer = clebsch_gordan(
				two_pair, 1, two_total, two_m_pair, two_m1, two_m_total);
			if (outer == 0.0) continue;
			for (int two_m2 : {-1, 1}) {
				for (int two_m3 : {-1, 1}) {
					const double inner = clebsch_gordan(
						1, 1, two_pair, two_m2, two_m3, two_m_pair);
					if (inner != 0.0) {
						result[product_index(two_m1, two_m2, two_m3)] += outer * inner;
					}
				}
			}
		}
	}
	return result;
}

struct spin_axis {
	int particle;
	int axis;
};

state8 apply_pauli_axis(const state8& input, int particle, int axis)
{
	state8 result = zero_state();
	const int stride = particle == 1 ? 4 : (particle == 2 ? 2 : 1);
	for (int index = 0; index < 8; ++index) {
		const bool down = (index & stride) != 0;
		if (axis == 0) {
			result[index ^ stride] += input[index];
		} else if (axis == 1) {
			result[index ^ stride] += (down ? -imaginary_unit : imaginary_unit) * input[index];
		} else {
			result[index] += (down ? -1.0 : 1.0) * input[index];
		}
	}
	return result;
}

state8 apply_spin_axes(state8 state, const std::vector<spin_axis>& axes)
{
	for (const spin_axis& operation : axes) {
		state = apply_pauli_axis(state, operation.particle, operation.axis);
	}
	return state;
}

complex isospin_matrix_element(const ls_channel& bra, const ls_channel& ket, int kind)
{
	const int two_m_t = bra.two_total_T;
	const state8 bra_state = coupled_three_half_state(bra.t_pair, bra.two_total_T, two_m_t);
	const state8 ket_state = coupled_three_half_state(ket.t_pair, ket.two_total_T, two_m_t);
	state8 operated = zero_state();
	if (kind >= 0 && kind <= 2) {
		const int first = kind == 0 ? 2 : 1;
		const int second = kind == 0 ? 3 : (kind == 1 ? 2 : 3);
		for (int axis = 0; axis < 3; ++axis) {
			state8 term = apply_pauli_axis(ket_state, second, axis);
			term = apply_pauli_axis(term, first, axis);
			for (int index = 0; index < 8; ++index) operated[index] += term[index];
		}
	} else {
		const int permutations[6][4] = {
			{0, 1, 2, +1}, {1, 2, 0, +1}, {2, 0, 1, +1},
			{1, 0, 2, -1}, {2, 1, 0, -1}, {0, 2, 1, -1}
		};
		for (const auto& permutation : permutations) {
			state8 term = apply_pauli_axis(ket_state, 3, permutation[2]);
			term = apply_pauli_axis(term, 2, permutation[1]);
			term = apply_pauli_axis(term, 1, permutation[0]);
			for (int index = 0; index < 8; ++index) {
				operated[index] += static_cast<double>(permutation[3]) * term[index];
			}
		}
	}
	return inner_product(bra_state, operated);
}

complex spherical_harmonic(int l_value, int m_value, const vector3& direction)
{
	const double radius = std::sqrt(
		direction[0] * direction[0] + direction[1] * direction[1]
		+ direction[2] * direction[2]);
	const double cosine = std::max(-1.0, std::min(1.0, direction[2] / radius));
	const double phi = std::atan2(direction[1], direction[0]);
	if (m_value >= 0) {
		return gsl_sf_legendre_sphPlm(l_value, m_value, cosine)
		     * std::exp(imaginary_unit * static_cast<double>(m_value) * phi);
	}
	const int positive_m = -m_value;
	const complex positive = gsl_sf_legendre_sphPlm(l_value, positive_m, cosine)
	                       * std::exp(imaginary_unit * static_cast<double>(positive_m) * phi);
	return (positive_m % 2 == 0 ? 1.0 : -1.0) * std::conj(positive);
}

complex bipolar_harmonic(int l_first, int l_second, int total_l,
	                       const vector3& first, const vector3& second)
{
	complex result{0.0, 0.0};
	for (int m_first = -l_first; m_first <= l_first; ++m_first) {
		const int m_second = -m_first;
		if (std::abs(m_second) > l_second) continue;
		const double coefficient = clebsch_gordan(
			2 * l_first, 2 * l_second, 2 * total_l,
			2 * m_first, 2 * m_second, 0);
		result += coefficient
		        * spherical_harmonic(l_first, m_first, first)
		        * spherical_harmonic(l_second, m_second, second);
	}
	return result;
}

struct harmonic_term {
	int l_value;
	int m_value;
	complex coefficient;
};

complex cartesian_y1_coefficient(int axis, int mu)
{
	if (axis == 0 && mu == -1) return std::sqrt(2.0 * pi_value / 3.0);
	if (axis == 0 && mu == +1) return -std::sqrt(2.0 * pi_value / 3.0);
	if (axis == 1 && (mu == -1 || mu == +1)) {
		return imaginary_unit * std::sqrt(2.0 * pi_value / 3.0);
	}
	if (axis == 2 && mu == 0) return std::sqrt(4.0 * pi_value / 3.0);
	return 0.0;
}

std::vector<harmonic_term> multiply_once(int l_value, int m_value, int axis)
{
	std::map<std::pair<int, int>, complex> result;
	for (int mu = -1; mu <= 1; ++mu) {
		const complex cartesian = cartesian_y1_coefficient(axis, mu);
		if (cartesian == complex{0.0, 0.0}) continue;
		const int new_m = m_value + mu;
		for (int new_l = std::abs(l_value - 1); new_l <= l_value + 1; ++new_l) {
			if (std::abs(new_m) > new_l) continue;
			const double product = std::sqrt(
				3.0 / (4.0 * pi_value) * (2.0 * l_value + 1.0) / (2.0 * new_l + 1.0))
				* clebsch_gordan(2 * l_value, 2, 2 * new_l, 0, 0, 0)
				* clebsch_gordan(2 * l_value, 2, 2 * new_l,
				                  2 * m_value, 2 * mu, 2 * new_m);
			const complex value = cartesian * product;
			if (std::abs(value) > 1.0e-15) result[{new_l, new_m}] += value;
		}
	}
	std::vector<harmonic_term> terms;
	for (const auto& item : result) terms.push_back({item.first.first, item.first.second, item.second});
	return terms;
}

std::vector<harmonic_term> multiply_axes(
	int l_value, int m_value, const std::vector<int>& axes, bool conjugate)
{
	std::vector<harmonic_term> states{{l_value, m_value, 1.0}};
	for (int axis : axes) {
		std::map<std::pair<int, int>, complex> updated;
		for (const harmonic_term& state : states) {
			for (const harmonic_term& next : multiply_once(state.l_value, state.m_value, axis)) {
				updated[{next.l_value, next.m_value}] += state.coefficient * next.coefficient;
			}
		}
		states.clear();
		for (const auto& item : updated) {
			states.push_back({item.first.first, item.first.second, item.second});
		}
	}
	if (conjugate) {
		for (harmonic_term& state : states) state.coefficient = std::conj(state.coefficient);
	}
	return states;
}

struct harmonic_cache_key {
	int l_value;
	int m_value;
	std::vector<int> axes;
	bool conjugate;

	bool operator<(const harmonic_cache_key& other) const
	{
		return std::tie(l_value, m_value, axes, conjugate)
		     < std::tie(other.l_value, other.m_value, other.axes, other.conjugate);
	}
};

const std::vector<harmonic_term>& get_harmonic_terms(
	int l_value, int m_value, const std::vector<int>& axes, bool conjugate)
{
	static thread_local std::map<harmonic_cache_key, std::vector<harmonic_term>> cache;
	const harmonic_cache_key key{l_value, m_value, axes, conjugate};
	const auto found = cache.find(key);
	if (found != cache.end()) return found->second;
	return cache.emplace(key, multiply_axes(l_value, m_value, axes, conjugate))
	            .first->second;
}

struct channel_term {
	int m_pair;
	int m_spectator;
	int spin_index;
	complex coefficient;
};

std::vector<channel_term> channel_terms(const ls_channel& channel, int two_m_j)
{
	std::map<std::array<int, 3>, complex> accumulated;
	for (int m_pair = -channel.l_pair; m_pair <= channel.l_pair; ++m_pair) {
		for (int m_spectator = -channel.lambda_spectator;
		     m_spectator <= channel.lambda_spectator; ++m_spectator) {
			const int m_total_l = m_pair + m_spectator;
			const double orbital = clebsch_gordan(
				2 * channel.l_pair, 2 * channel.lambda_spectator, 2 * channel.total_L,
				2 * m_pair, 2 * m_spectator, 2 * m_total_l);
			if (orbital == 0.0) continue;
			for (int two_m_pair_spin = -2 * channel.s_pair;
			     two_m_pair_spin <= 2 * channel.s_pair; two_m_pair_spin += 2) {
				for (int two_m_one : {-1, 1}) {
					const int two_m_s = two_m_pair_spin + two_m_one;
					const double spin = clebsch_gordan(
						2 * channel.s_pair, 1, channel.two_total_S,
						two_m_pair_spin, two_m_one, two_m_s);
					const double total = clebsch_gordan(
						2 * channel.total_L, channel.two_total_S, channel.two_total_J,
						2 * m_total_l, two_m_s, two_m_j);
					if (spin == 0.0 || total == 0.0) continue;
					for (int two_m_two : {-1, 1}) {
						for (int two_m_three : {-1, 1}) {
							const double pair = clebsch_gordan(
								1, 1, 2 * channel.s_pair,
								two_m_two, two_m_three, two_m_pair_spin);
							if (pair == 0.0) continue;
							const int index = product_index(two_m_one, two_m_two, two_m_three);
							accumulated[{m_pair, m_spectator, index}] +=
								orbital * spin * total * pair;
						}
					}
				}
			}
		}
	}
	std::vector<channel_term> result;
	for (const auto& item : accumulated) {
		if (std::abs(item.second) > 1.0e-15) {
			result.push_back({item.first[0], item.first[1], item.first[2], item.second});
		}
	}
	return result;
}

const std::vector<channel_term>& get_channel_terms(
	const ls_channel& channel, int two_m_j)
{
	using key_type = std::array<int, 9>;
	static thread_local std::map<key_type, std::vector<channel_term>> cache;
	key_type key{};
	const auto channel_values = channel_key(channel);
	std::copy(channel_values.begin(), channel_values.end(), key.begin());
	key[8] = two_m_j;
	const auto found = cache.find(key);
	if (found != cache.end()) return found->second;
	return cache.emplace(key, channel_terms(channel, two_m_j)).first->second;
}

enum coordinate_id { p_bra = 0, q_bra = 1, p_ket = 2, q_ket = 3 };

struct coordinate_axis {
	coordinate_id coordinate;
	int axis;
};

struct cartesian_term {
	complex coefficient;
	std::vector<coordinate_axis> coordinate_axes;
	std::vector<spin_axis> spin_axes;
};

std::vector<cartesian_term> q_transfer_dot_terms(int second_particle)
{
	const std::array<std::pair<coordinate_id, double>, 2> sources{{
		{q_bra, +1.0}, {q_ket, -1.0}
	}};
	std::vector<cartesian_term> result;
	for (const auto& first : sources) {
		for (const auto& second : sources) {
			for (int first_axis = 0; first_axis < 3; ++first_axis) {
				for (int second_axis = 0; second_axis < 3; ++second_axis) {
					result.push_back({
						first.second * second.second,
						{{first.first, first_axis}, {second.first, second_axis}},
						{{1, first_axis}, {second_particle, second_axis}}
					});
				}
			}
		}
	}
	return result;
}

std::vector<cartesian_term> q2_q3_dot_terms()
{
	const std::array<std::pair<coordinate_id, double>, 4> q2_sources{{
		{p_bra, +1.0}, {p_ket, -1.0}, {q_bra, -0.5}, {q_ket, +0.5}
	}};
	const std::array<std::pair<coordinate_id, double>, 4> q3_sources{{
		{p_bra, -1.0}, {p_ket, +1.0}, {q_bra, -0.5}, {q_ket, +0.5}
	}};
	std::vector<cartesian_term> result;
	for (const auto& q2 : q2_sources) {
		for (const auto& q3 : q3_sources) {
			for (int q2_axis = 0; q2_axis < 3; ++q2_axis) {
				for (int q3_axis = 0; q3_axis < 3; ++q3_axis) {
					result.push_back({
						q2.second * q3.second,
						{{q2.first, q2_axis}, {q3.first, q3_axis}},
						{{2, q2_axis}, {3, q3_axis}}
					});
				}
			}
		}
	}
	return result;
}

std::vector<cartesian_term> c4_cartesian_terms()
{
	const std::array<std::pair<coordinate_id, double>, 4> q2_sources{{
		{p_bra, +1.0}, {p_ket, -1.0}, {q_bra, -0.5}, {q_ket, +0.5}
	}};
	const std::array<std::pair<coordinate_id, double>, 4> q3_sources{{
		{p_bra, -1.0}, {p_ket, +1.0}, {q_bra, -0.5}, {q_ket, +0.5}
	}};
	const std::array<std::pair<coordinate_id, double>, 2> dp_sources{{
		{p_bra, +1.0}, {p_ket, -1.0}
	}};
	const std::array<std::pair<coordinate_id, double>, 2> dq_sources{{
		{q_bra, +1.0}, {q_ket, -1.0}
	}};
	const int epsilon[6][4] = {
		{0, 1, 2, +1}, {1, 2, 0, +1}, {2, 0, 1, +1},
		{1, 0, 2, -1}, {2, 1, 0, -1}, {0, 2, 1, -1}
	};
	std::vector<cartesian_term> result;
	result.reserve(3456);
	for (const auto& q2 : q2_sources) for (const auto& q3 : q3_sources)
	for (const auto& dp : dp_sources) for (const auto& dq : dq_sources)
	for (int q2_axis = 0; q2_axis < 3; ++q2_axis)
	for (int q3_axis = 0; q3_axis < 3; ++q3_axis)
	for (const auto& entry : epsilon) {
		result.push_back({
			-q2.second * q3.second * dp.second * dq.second * entry[3],
			{{q2.first, q2_axis}, {q3.first, q3_axis},
			 {dp.first, entry[0]}, {dq.first, entry[1]}},
			{{2, q2_axis}, {3, q3_axis}, {1, entry[2]}}
		});
	}
	return result;
}

enum class algebra_kind { identity = 0, q23 = 1, c4 = 2, d12 = 3, d13 = 4 };

const std::vector<cartesian_term>& cartesian_terms(algebra_kind kind)
{
	static const std::vector<cartesian_term> identity_terms{{1.0, {}, {}}};
	static const std::vector<cartesian_term> q23_terms = q2_q3_dot_terms();
	static const std::vector<cartesian_term> c4_terms = c4_cartesian_terms();
	static const std::vector<cartesian_term> d12_terms = q_transfer_dot_terms(2);
	static const std::vector<cartesian_term> d13_terms = q_transfer_dot_terms(3);
	switch (kind) {
		case algebra_kind::identity: return identity_terms;
		case algebra_kind::q23: return q23_terms;
		case algebra_kind::c4: return c4_terms;
		case algebra_kind::d12: return d12_terms;
		case algebra_kind::d13: return d13_terms;
	}
	return identity_terms;
}

state8 spin_basis_state(int index)
{
	state8 result = zero_state();
	result[index] = 1.0;
	return result;
}

complex spin_matrix_element(int bra_index, int ket_index,
	                          const std::vector<spin_axis>& axes)
{
	const state8 operated = apply_spin_axes(spin_basis_state(ket_index), axes);
	return operated[bra_index];
}

struct spin_cache_key {
	int bra_index;
	int ket_index;
	std::vector<int> axes;

	bool operator<(const spin_cache_key& other) const
	{
		return std::tie(bra_index, ket_index, axes)
		     < std::tie(other.bra_index, other.ket_index, other.axes);
	}
};

complex get_spin_matrix_element(int bra_index, int ket_index,
	                              const std::vector<spin_axis>& axes)
{
	static thread_local std::map<spin_cache_key, complex> cache;
	std::vector<int> encoded_axes;
	encoded_axes.reserve(axes.size());
	for (const spin_axis& axis : axes) {
		encoded_axes.push_back(3 * (axis.particle - 1) + axis.axis);
	}
	const spin_cache_key key{bra_index, ket_index, std::move(encoded_axes)};
	const auto found = cache.find(key);
	if (found != cache.end()) return found->second;
	return cache.emplace(key, spin_matrix_element(bra_index, ket_index, axes))
	            .first->second;
}

struct weight_key {
	std::array<int, 8> angular;
	std::array<int, 4> radial_powers;

	bool operator<(const weight_key& other) const
	{
		return std::tie(angular, radial_powers)
		     < std::tie(other.angular, other.radial_powers);
	}
};

using weight_table = std::vector<std::pair<weight_key, complex>>;

weight_table build_weight_table(
	const ls_channel& bra, const ls_channel& ket, algebra_kind kind)
{
	if (bra.two_total_J != ket.two_total_J) return {};
	std::map<weight_key, complex> accumulated;
	const double inverse_m_count = 1.0 / (bra.two_total_J + 1.0);
	for (int two_m_j = -bra.two_total_J;
	     two_m_j <= bra.two_total_J; two_m_j += 2) {
		const auto& bra_channel_terms = get_channel_terms(bra, two_m_j);
		const auto& ket_channel_terms = get_channel_terms(ket, two_m_j);
		for (const cartesian_term& cartesian : cartesian_terms(kind)) {
			std::array<std::vector<int>, 4> axes_by_slot;
			std::array<int, 4> powers{{0, 0, 0, 0}};
			for (const coordinate_axis& factor : cartesian.coordinate_axes) {
				const int slot = static_cast<int>(factor.coordinate);
				axes_by_slot[slot].push_back(factor.axis);
				++powers[slot];
			}
			for (const channel_term& bra_term : bra_channel_terms) {
				const auto& bra_pair = get_harmonic_terms(
					bra.l_pair, bra_term.m_pair, axes_by_slot[p_bra], true);
				const auto& bra_spectator = get_harmonic_terms(
					bra.lambda_spectator, bra_term.m_spectator,
					axes_by_slot[q_bra], true);
				for (const channel_term& ket_term : ket_channel_terms) {
					const complex spin = get_spin_matrix_element(
						bra_term.spin_index, ket_term.spin_index, cartesian.spin_axes);
					if (std::abs(spin) < 1.0e-15) continue;
					const auto& ket_pair = get_harmonic_terms(
						ket.l_pair, ket_term.m_pair, axes_by_slot[p_ket], false);
					const auto& ket_spectator = get_harmonic_terms(
						ket.lambda_spectator, ket_term.m_spectator,
						axes_by_slot[q_ket], false);
					const complex state_coefficient = cartesian.coefficient
						* std::conj(bra_term.coefficient) * ket_term.coefficient
						* spin * inverse_m_count;
					for (const harmonic_term& bp : bra_pair)
					for (const harmonic_term& bq : bra_spectator)
					for (const harmonic_term& kp : ket_pair)
					for (const harmonic_term& kq : ket_spectator) {
						weight_key key{{
							bp.l_value, bp.m_value, bq.l_value, bq.m_value,
							kp.l_value, kp.m_value, kq.l_value, kq.m_value
						}, powers};
						accumulated[key] += state_coefficient
							* bp.coefficient * bq.coefficient
							* kp.coefficient * kq.coefficient;
					}
				}
			}
		}
	}
	weight_table result;
	result.reserve(accumulated.size());
	for (const auto& item : accumulated) {
		if (std::abs(item.second) > 2.0e-14) result.push_back(item);
	}
	return result;
}

std::array<int, 17> weight_cache_key(
	const ls_channel& bra, const ls_channel& ket, algebra_kind kind)
{
	const auto bra_key = channel_key(bra);
	const auto ket_key = channel_key(ket);
	std::array<int, 17> result{};
	std::copy(bra_key.begin(), bra_key.end(), result.begin());
	std::copy(ket_key.begin(), ket_key.end(), result.begin() + 8);
	result[16] = static_cast<int>(kind);
	return result;
}

std::shared_ptr<const weight_table> get_weight_table(
	const ls_channel& bra, const ls_channel& ket, algebra_kind kind)
{
	using key_type = std::array<int, 17>;
	using table_ptr = std::shared_ptr<const weight_table>;
	using table_future = std::shared_future<table_ptr>;
	static std::mutex cache_mutex;
	static std::map<key_type, table_future> cache;
	const key_type key = weight_cache_key(bra, ket, kind);
	table_future result;
	std::shared_ptr<std::promise<table_ptr>> builder;
	{
		std::lock_guard<std::mutex> lock(cache_mutex);
		const auto found = cache.find(key);
		if (found != cache.end()) {
			result = found->second;
		} else {
			builder = std::make_shared<std::promise<table_ptr>>();
			result = builder->get_future().share();
			cache.emplace(key, result);
		}
	}
	if (!builder) return result.get();
	try {
		table_ptr candidate = std::make_shared<weight_table>(
			build_weight_table(bra, ket, kind));
		builder->set_value(candidate);
		return candidate;
	} catch (...) {
		builder->set_exception(std::current_exception());
		throw;
	}
}

struct quadrature_grid {
	std::vector<double> nodes;
	std::vector<double> weights;
};

std::shared_ptr<const quadrature_grid> get_quadrature_grid(int order)
{
	static std::mutex cache_mutex;
	static std::map<int, std::shared_ptr<const quadrature_grid>> cache;
	{
		std::lock_guard<std::mutex> lock(cache_mutex);
		const auto found = cache.find(order);
		if (found != cache.end()) return found->second;
	}
	auto candidate = std::make_shared<quadrature_grid>();
	candidate->nodes.resize(order);
	candidate->weights.resize(order);
	gauss(candidate->nodes.data(), candidate->weights.data(), order);
	std::lock_guard<std::mutex> lock(cache_mutex);
	return cache.emplace(order, candidate).first->second;
}

enum class scalar_kernel_kind { contact = 0, two_pion = 1, c3 = 2, one_pion = 3 };

double legendre_polynomial(int l_value, double x)
{
	return gsl_sf_legendre_Pl(l_value, x);
}

complex scalar_kernel_value(
	scalar_kernel_kind kind, double delta_p, double delta_q,
	double relative_cosine, double pion_mass)
{
	if (kind == scalar_kernel_kind::contact) return 1.0;
	if (kind == scalar_kernel_kind::one_pion) {
		return 1.0 / (delta_q * delta_q + pion_mass * pion_mass);
	}
	const double common = delta_p * delta_p + 0.25 * delta_q * delta_q
	                    + pion_mass * pion_mass;
	const double d2 = common - delta_p * delta_q * relative_cosine;
	const double d3 = common + delta_p * delta_q * relative_cosine;
	const double propagators = 1.0 / (d2 * d3);
	if (kind == scalar_kernel_kind::c3) {
		return (-delta_p * delta_p + 0.25 * delta_q * delta_q) * propagators;
	}
	return propagators;
}

complex uncoupled_orbital_kernel(
	const std::array<int, 8>& angular,
	double p, double q, double pp, double qp,
	scalar_kernel_kind kind, int order, double pion_mass)
{
	const int bra_l_pair = angular[0];
	const int bra_m_pair = angular[1];
	const int bra_lambda = angular[2];
	const int bra_m_lambda = angular[3];
	const int ket_l_pair = angular[4];
	const int ket_m_pair = angular[5];
	const int ket_lambda = angular[6];
	const int ket_m_lambda = angular[7];
	if (ket_m_pair - bra_m_pair != bra_m_lambda - ket_m_lambda) return 0.0;
	const int lbar_min = std::max(std::abs(bra_l_pair - ket_l_pair),
	                              std::abs(bra_lambda - ket_lambda));
	const int lbar_max = std::min(bra_l_pair + ket_l_pair,
	                              bra_lambda + ket_lambda);
	if (lbar_min > lbar_max) return 0.0;
	const auto grid = get_quadrature_grid(order);
	const double dp_lo = std::abs(pp - p);
	const double dp_hi = pp + p;
	const double dq_lo = std::abs(qp - q);
	const double dq_hi = qp + q;
	const double phase = ((ket_m_pair + bra_m_lambda) % 2 == 0) ? 1.0 : -1.0;
	const double prefactor = phase * 2.0 * std::pow(2.0 * pi_value, 4)
	                       / (p * pp * q * qp);
	complex total{0.0, 0.0};
	for (int lbar = lbar_min; lbar <= lbar_max; ++lbar) {
		const double pair_cg = clebsch_gordan(
			2 * bra_l_pair, 2 * ket_l_pair, 2 * lbar,
			-2 * bra_m_pair, 2 * ket_m_pair,
			2 * (-bra_m_pair + ket_m_pair));
		const double spectator_cg = clebsch_gordan(
			2 * bra_lambda, 2 * ket_lambda, 2 * lbar,
			-2 * bra_m_lambda, 2 * ket_m_lambda,
			2 * (-bra_m_lambda + ket_m_lambda));
		const double angular_coefficient = pair_cg * spectator_cg / (2.0 * lbar + 1.0);
		if (angular_coefficient == 0.0) continue;
		complex integral{0.0, 0.0};
		for (int idp = 0; idp < order; ++idp) {
			const double delta_p = 0.5 * ((dp_hi - dp_lo) * grid->nodes[idp] + dp_hi + dp_lo);
			const double weight_p = 0.5 * (dp_hi - dp_lo) * grid->weights[idp];
			const double cosine_p = std::max(-1.0, std::min(1.0,
				(pp * pp - p * p - delta_p * delta_p) / (2.0 * delta_p * p)));
			const double sine_p = std::sqrt(std::max(0.0, 1.0 - cosine_p * cosine_p));
			const vector3 p_hat{sine_p, 0.0, cosine_p};
			const vector3 pp_vector{p * sine_p, 0.0, p * cosine_p + delta_p};
			const double pp_norm = std::sqrt(
				pp_vector[0] * pp_vector[0] + pp_vector[2] * pp_vector[2]);
			const vector3 pp_hat{pp_vector[0] / pp_norm, 0.0, pp_vector[2] / pp_norm};
			const complex pair_bipolar = bipolar_harmonic(
				bra_l_pair, ket_l_pair, lbar, pp_hat, p_hat);
			for (int idq = 0; idq < order; ++idq) {
				const double delta_q = 0.5 * ((dq_hi - dq_lo) * grid->nodes[idq] + dq_hi + dq_lo);
				const double weight_q = 0.5 * (dq_hi - dq_lo) * grid->weights[idq];
				const double cosine_q = std::max(-1.0, std::min(1.0,
					(qp * qp - q * q - delta_q * delta_q) / (2.0 * delta_q * q)));
				const double sine_q = std::sqrt(std::max(0.0, 1.0 - cosine_q * cosine_q));
				const vector3 q_hat{sine_q, 0.0, cosine_q};
				const vector3 qp_vector{q * sine_q, 0.0, q * cosine_q + delta_q};
				const double qp_norm = std::sqrt(
					qp_vector[0] * qp_vector[0] + qp_vector[2] * qp_vector[2]);
				const vector3 qp_hat{qp_vector[0] / qp_norm, 0.0, qp_vector[2] / qp_norm};
				const complex spectator_bipolar = bipolar_harmonic(
					bra_lambda, ket_lambda, lbar, qp_hat, q_hat);
				complex relative_integral{0.0, 0.0};
				for (int ix = 0; ix < order; ++ix) {
					const double x = grid->nodes[ix];
					relative_integral += grid->weights[ix]
						* legendre_polynomial(lbar, x)
						* scalar_kernel_value(kind, delta_p, delta_q, x, pion_mass);
				}
				integral += weight_p * delta_p * weight_q * delta_q
					* pair_bipolar * spectator_bipolar * relative_integral;
			}
		}
		total += angular_coefficient * integral;
	}
	return prefactor * total;
}

struct orbital_cache_key {
	std::array<int, 8> angular;
	scalar_kernel_kind kind;

	bool operator<(const orbital_cache_key& other) const
	{
		return std::tie(angular, kind) < std::tie(other.angular, other.kind);
	}
};

using orbital_cache = std::map<orbital_cache_key, complex>;

complex project_cartesian_algebra(
	const ls_channel& bra, const ls_channel& ket,
	double p, double q, double pp, double qp,
	algebra_kind algebra, scalar_kernel_kind kernel,
	int order, double pion_mass, orbital_cache& cache)
{
	const auto weights = get_weight_table(bra, ket, algebra);
	const std::array<double, 4> radial{{pp, qp, p, q}};
	complex result{0.0, 0.0};
	for (const auto& item : *weights) {
		const weight_key& key = item.first;
		double radial_factor = 1.0;
		for (int slot = 0; slot < 4; ++slot) {
			radial_factor *= std::pow(radial[slot], key.radial_powers[slot]);
		}
		const orbital_cache_key cache_key{key.angular, kernel};
		auto found = cache.find(cache_key);
		if (found == cache.end()) {
			found = cache.emplace(
				cache_key,
				uncoupled_orbital_kernel(
					key.angular, p, q, pp, qp, kernel, order, pion_mass)
			).first;
		}
		result += item.second * radial_factor * found->second;
	}
	return result;
}

double regulator(double p, double q, double cutoff)
{
	const double invariant = p * p + 0.75 * q * q;
	const double ratio = invariant / (cutoff * cutoff);
	return std::exp(-(ratio * ratio));
}

std::array<complex, 5> project_ls_components(
	const ls_channel& bra, const ls_channel& ket,
	double p, double q, double pp, double qp,
	double c_D, double c_E, double c1_fm, double c3_fm, double c4_fm,
	int order, orbital_cache& cache)
{
	std::array<complex, 5> result{};
	result.fill(0.0);
	const double f_pi = fpi / hbarc;
	const double pion_mass = mpi / hbarc;
	const double lambda_chi = 700.0 / hbarc;
	const complex iso23 = isospin_matrix_element(bra, ket, 0);
	const complex iso12 = isospin_matrix_element(bra, ket, 1);
	const complex iso13 = isospin_matrix_element(bra, ket, 2);
	const complex iso_cross = isospin_matrix_element(bra, ket, 3);
	const double common = gA * gA / (4.0 * std::pow(f_pi, 4));
	if ((c1_fm != 0.0 || c3_fm != 0.0) && std::abs(iso23) > 1.0e-15) {
		if (c1_fm != 0.0) {
			const complex spin = project_cartesian_algebra(
				bra, ket, p, q, pp, qp, algebra_kind::q23,
				scalar_kernel_kind::two_pion, order, pion_mass, cache);
			result[0] = common * (-4.0 * c1_fm * pion_mass * pion_mass) * iso23 * spin;
		}
		if (c3_fm != 0.0) {
			const complex spin = project_cartesian_algebra(
				bra, ket, p, q, pp, qp, algebra_kind::q23,
				scalar_kernel_kind::c3, order, pion_mass, cache);
			result[1] = common * (2.0 * c3_fm) * iso23 * spin;
		}
	}
	if (c4_fm != 0.0 && std::abs(iso_cross) > 1.0e-15) {
		const complex spin = project_cartesian_algebra(
			bra, ket, p, q, pp, qp, algebra_kind::c4,
			scalar_kernel_kind::two_pion, order, pion_mass, cache);
		result[2] = common * c4_fm * iso_cross * spin;
	}
	if (c_D != 0.0) {
		complex spin12{0.0, 0.0};
		complex spin13{0.0, 0.0};
		if (std::abs(iso12) > 1.0e-15) {
			spin12 = project_cartesian_algebra(
				bra, ket, p, q, pp, qp, algebra_kind::d12,
				scalar_kernel_kind::one_pion, order, pion_mass, cache);
		}
		if (std::abs(iso13) > 1.0e-15) {
			spin13 = project_cartesian_algebra(
				bra, ket, p, q, pp, qp, algebra_kind::d13,
				scalar_kernel_kind::one_pion, order, pion_mass, cache);
		}
		const double d_lec = c_D / (f_pi * f_pi * lambda_chi);
		result[3] = -gA * d_lec / (8.0 * f_pi * f_pi)
		          * (iso12 * spin12 + iso13 * spin13);
	}
	if (c_E != 0.0 && std::abs(iso23) > 1.0e-15) {
		const complex scalar = project_cartesian_algebra(
			bra, ket, p, q, pp, qp, algebra_kind::identity,
			scalar_kernel_kind::contact, order, pion_mass, cache);
		const double e_lec = c_E / (std::pow(f_pi, 4) * lambda_chi);
		result[4] = e_lec * iso23 * scalar;
	}
	return result;
}

} // namespace

chiral_N2LO_3NF_factorized::chiral_N2LO_3NF_factorized(
	double c_D, double c_E, double Lambda_3NF_MeV,
	double c1, double c3, double c4, int transfer_order)
	: m_c_D(c_D), m_c_E(c_E), m_lambda(Lambda_3NF_MeV / hbarc),
	  m_c1_gev(c1), m_c3_gev(c3), m_c4_gev(c4),
	  m_transfer_order(transfer_order)
{
	if (m_lambda <= 0.0) throw std::invalid_argument("factorized 3NF cutoff must be positive");
	if (m_transfer_order < 1) {
		throw std::invalid_argument("factorized 3NF transfer quadrature order must be positive");
	}
}

bool chiral_N2LO_3NF_factorized::enabled() const
{
	return m_c_D != 0.0 || m_c_E != 0.0 || m_c1_gev != 0.0
	    || m_c3_gev != 0.0 || m_c4_gev != 0.0;
}

void chiral_N2LO_3NF_factorized::update_parameters(const double* parameters)
{
	if (parameters != nullptr) {
		m_c_D = parameters[0];
		m_c_E = parameters[1];
	}
}

double chiral_N2LO_3NF_factorized::W1_element(
	int alpha_r, int alpha_c, double p_r, double q_r,
	double p_c, double q_c, const pw_3N_statespace& pw_states) const
{
	const jj_channel bra = make_jj_channel(alpha_r, pw_states);
	const jj_channel ket = make_jj_channel(alpha_c, pw_states);
	if (bra.two_total_J != ket.two_total_J || bra.two_total_T != ket.two_total_T) return 0.0;
	const auto bra_expansion = get_ls_expansion(bra);
	const auto ket_expansion = get_ls_expansion(ket);
	std::array<complex, 5> totals{};
	totals.fill(0.0);
	orbital_cache cache;
	const double c1_fm = m_c1_gev * hbarc / 1000.0;
	const double c3_fm = m_c3_gev * hbarc / 1000.0;
	const double c4_fm = m_c4_gev * hbarc / 1000.0;
	for (const auto& bra_ls : *bra_expansion) {
		for (const auto& ket_ls : *ket_expansion) {
			const auto components = project_ls_components(
				bra_ls.first, ket_ls.first, p_c, q_c, p_r, q_r,
				m_c_D, m_c_E, c1_fm, c3_fm, c4_fm, m_transfer_order, cache);
			const double coefficient = bra_ls.second * ket_ls.second;
			for (int component = 0; component < 5; ++component) {
				totals[component] += coefficient * components[component];
			}
		}
	}
	complex total{0.0, 0.0};
	for (complex value : totals) total += value;
	const double fourier_normalization = std::pow(2.0 * pi_value, -6);
	total *= fourier_normalization
	       * regulator(p_r, q_r, m_lambda) * regulator(p_c, q_c, m_lambda);
	const double tolerance = 2.0e-9 * std::max(1.0, std::abs(total.real()));
	if (std::abs(total.imag()) > tolerance) {
		throw std::runtime_error("factorized complete 3NF produced a non-real scalar");
	}
	return total.real();
}

double chiral_N2LO_3NF_factorized::axial_coupling_3nf() const { return gA; }
double chiral_N2LO_3NF_factorized::pion_decay_constant_mev_3nf() const { return fpi; }
double chiral_N2LO_3NF_factorized::pion_mass_mev_3nf() const { return mpi; }
double chiral_N2LO_3NF_factorized::chiral_scale_mev_3nf() const { return 700.0; }
double chiral_N2LO_3NF_factorized::hbarc_mev_fm_3nf() const { return hbarc; }
