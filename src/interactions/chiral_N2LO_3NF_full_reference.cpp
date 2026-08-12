#include "chiral_N2LO_3NF_full_reference.h"

#include "constants.h"
#include "coupling_coefficients.h"
#include "gauss_legendre.h"

#include <gsl/gsl_sf_legendre.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <vector>

namespace {

using complex = std::complex<double>;
using vector3 = std::array<double, 3>;
using state8 = std::array<complex, 8>;

constexpr double pi_value = 3.141592653589793238462643383279502884;
constexpr complex imaginary_unit{0.0, 1.0};

struct channel {
	int l_pair;
	int s_pair;
	int j_pair;
	int lambda_spectator;
	int two_j_spectator;
	int two_total_J;
	int t_pair;
	int two_total_T;
};

vector3 subtract(const vector3& a, const vector3& b)
{
	return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

vector3 scale(double factor, const vector3& a)
{
	return {factor * a[0], factor * a[1], factor * a[2]};
}

double dot(const vector3& a, const vector3& b)
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

vector3 cross(const vector3& a, const vector3& b)
{
	return {
		a[1] * b[2] - a[2] * b[1],
		a[2] * b[0] - a[0] * b[2],
		a[0] * b[1] - a[1] * b[0]
	};
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
	for (int i = 0; i < 8; ++i) result += std::conj(bra[i]) * ket[i];
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

state8 pair_spectator_spin_product(int s_pair, int two_m_pair,
	                                 int two_m_spectator)
{
	state8 result = zero_state();
	for (int two_m2 : {-1, 1}) {
		for (int two_m3 : {-1, 1}) {
			const double coefficient = clebsch_gordan(
				1, 1, 2 * s_pair, two_m2, two_m3, two_m_pair);
			if (coefficient != 0.0) {
				result[product_index(two_m_spectator, two_m2, two_m3)] += coefficient;
			}
		}
	}
	return result;
}

complex spherical_harmonic(int l, int m, double cosine, double phi)
{
	const double x = std::max(-1.0, std::min(1.0, cosine));
	if (m >= 0) {
		return gsl_sf_legendre_sphPlm(l, m, x)
		     * std::exp(imaginary_unit * static_cast<double>(m) * phi);
	}
	const int positive_m = -m;
	const complex positive = gsl_sf_legendre_sphPlm(l, positive_m, x)
	                       * std::exp(imaginary_unit * static_cast<double>(positive_m) * phi);
	return (positive_m % 2 == 0 ? 1.0 : -1.0) * std::conj(positive);
}

state8 angular_spin_state_jj(const channel& ch,
	                           double cosine_p, double phi_p,
	                           double cosine_q, double phi_q,
	                           int two_m_j)
{
	state8 result = zero_state();
	for (int m_l = -ch.l_pair; m_l <= ch.l_pair; ++m_l) {
		const complex y_l = spherical_harmonic(ch.l_pair, m_l, cosine_p, phi_p);
		for (int two_m_pair = -2 * ch.s_pair;
		     two_m_pair <= 2 * ch.s_pair; two_m_pair += 2) {
			const int two_m_j_pair = 2 * m_l + two_m_pair;
			if (std::abs(two_m_j_pair) > 2 * ch.j_pair) continue;
			const double pair_cg = clebsch_gordan(
				2 * ch.l_pair, 2 * ch.s_pair, 2 * ch.j_pair,
				2 * m_l, two_m_pair, two_m_j_pair);
			if (pair_cg == 0.0) continue;
			for (int m_lambda = -ch.lambda_spectator;
			     m_lambda <= ch.lambda_spectator; ++m_lambda) {
				const complex y_lambda = spherical_harmonic(
					ch.lambda_spectator, m_lambda, cosine_q, phi_q);
				for (int two_m_spectator : {-1, 1}) {
					const int two_m_j_spectator = 2 * m_lambda + two_m_spectator;
					if (std::abs(two_m_j_spectator) > ch.two_j_spectator) continue;
					const double spectator_cg = clebsch_gordan(
						2 * ch.lambda_spectator, 1, ch.two_j_spectator,
						2 * m_lambda, two_m_spectator, two_m_j_spectator);
					const double total_cg = clebsch_gordan(
						2 * ch.j_pair, ch.two_j_spectator, ch.two_total_J,
						two_m_j_pair, two_m_j_spectator, two_m_j);
					const complex coefficient = pair_cg * spectator_cg * total_cg
					                          * y_l * y_lambda;
					if (coefficient == complex{0.0, 0.0}) continue;
					const state8 spin = pair_spectator_spin_product(
						ch.s_pair, two_m_pair, two_m_spectator);
					for (int i = 0; i < 8; ++i) result[i] += coefficient * spin[i];
				}
			}
		}
	}
	return result;
}

struct angular_grid {
	int order;
	std::vector<double> cosine;
	std::vector<double> cosine_weight;
	std::vector<double> phi;
	std::vector<double> phi_weight;
};

std::shared_ptr<const angular_grid> get_angular_grid(int order)
{
	static std::mutex cache_mutex;
	static std::map<int, std::shared_ptr<const angular_grid>> cache;
	{
		std::lock_guard<std::mutex> lock(cache_mutex);
		const auto found = cache.find(order);
		if (found != cache.end()) return found->second;
	}

	auto candidate = std::make_shared<angular_grid>();
	candidate->order = order;
	candidate->cosine.resize(order);
	candidate->cosine_weight.resize(order);
	candidate->phi.resize(order);
	candidate->phi_weight.resize(order);
	gauss(candidate->cosine.data(), candidate->cosine_weight.data(), order);
	gauss(candidate->phi.data(), candidate->phi_weight.data(), order);
	for (int i = 0; i < order; ++i) {
		candidate->phi[i] = pi_value * (candidate->phi[i] + 1.0);
		candidate->phi_weight[i] *= pi_value;
	}

	std::lock_guard<std::mutex> lock(cache_mutex);
	const auto inserted = cache.emplace(order, candidate);
	return inserted.first->second;
}

struct angular_basis {
	int order;
	int num_m;
	std::vector<state8> ket_states;
	std::vector<state8> bra_states;

	const state8& ket(int iq, int im) const
	{
		return ket_states[static_cast<std::size_t>(iq) * num_m + im];
	}

	const state8& bra(int ipp, int iphip, int iqp, int iphiqp, int im) const
	{
		const std::size_t angular_index =
			((static_cast<std::size_t>(ipp) * order + iphip) * order + iqp) * order + iphiqp;
		return bra_states[angular_index * num_m + im];
	}
};

std::array<int, 9> basis_key(const channel& ch, int order)
{
	return {order, ch.l_pair, ch.s_pair, ch.j_pair, ch.lambda_spectator,
	        ch.two_j_spectator, ch.two_total_J, ch.t_pair, ch.two_total_T};
}

std::shared_ptr<const angular_basis> get_angular_basis(
	const channel& ch, const std::shared_ptr<const angular_grid>& grid)
{
	using key_type = std::array<int, 9>;
	static std::mutex cache_mutex;
	static std::map<key_type, std::shared_ptr<const angular_basis>> cache;
	const key_type key = basis_key(ch, grid->order);
	{
		std::lock_guard<std::mutex> lock(cache_mutex);
		const auto found = cache.find(key);
		if (found != cache.end()) return found->second;
	}

	auto candidate = std::make_shared<angular_basis>();
	candidate->order = grid->order;
	candidate->num_m = ch.two_total_J + 1;
	candidate->ket_states.reserve(
		static_cast<std::size_t>(grid->order) * candidate->num_m);
	for (int iq = 0; iq < grid->order; ++iq) {
		for (int two_m_j = -ch.two_total_J;
		     two_m_j <= ch.two_total_J; two_m_j += 2) {
			candidate->ket_states.push_back(angular_spin_state_jj(
				ch, 1.0, 0.0, grid->cosine[iq], 0.0, two_m_j));
		}
	}

	const std::size_t number_of_bra_angles = static_cast<std::size_t>(grid->order)
	                                             * grid->order * grid->order * grid->order;
	candidate->bra_states.reserve(number_of_bra_angles * candidate->num_m);
	for (int ipp = 0; ipp < grid->order; ++ipp) {
		for (int iphip = 0; iphip < grid->order; ++iphip) {
			for (int iqp = 0; iqp < grid->order; ++iqp) {
				for (int iphiqp = 0; iphiqp < grid->order; ++iphiqp) {
					for (int two_m_j = -ch.two_total_J;
					     two_m_j <= ch.two_total_J; two_m_j += 2) {
						candidate->bra_states.push_back(angular_spin_state_jj(
							ch, grid->cosine[ipp], grid->phi[iphip],
							grid->cosine[iqp], grid->phi[iphiqp], two_m_j));
					}
				}
			}
		}
	}

	std::lock_guard<std::mutex> lock(cache_mutex);
	const auto inserted = cache.emplace(key, candidate);
	return inserted.first->second;
}

state8 apply_pauli_axis(const state8& input, int particle, int axis)
{
	state8 result = zero_state();
	const int mask = 1 << (3 - particle);
	for (int index = 0; index < 8; ++index) {
		const bool down = (index & mask) != 0;
		if (axis == 0) {
			result[index ^ mask] += input[index];
		} else if (axis == 1) {
			result[index ^ mask] += (down ? -imaginary_unit : imaginary_unit) * input[index];
		} else {
			result[index] += (down ? -1.0 : 1.0) * input[index];
		}
	}
	return result;
}

state8 apply_pair_dot(const state8& input, int particle_a, int particle_b)
{
	state8 result = zero_state();
	for (int axis = 0; axis < 3; ++axis) {
		const state8 first = apply_pauli_axis(input, particle_b, axis);
		const state8 second = apply_pauli_axis(first, particle_a, axis);
		for (int i = 0; i < 8; ++i) result[i] += second[i];
	}
	return result;
}

state8 apply_triple_cross_123(const state8& input)
{
	state8 result = zero_state();
	const int permutations[6][4] = {
		{0, 1, 2, +1}, {1, 2, 0, +1}, {2, 0, 1, +1},
		{1, 0, 2, -1}, {2, 1, 0, -1}, {0, 2, 1, -1}
	};
	for (const auto& permutation : permutations) {
		state8 term = apply_pauli_axis(input, 3, permutation[2]);
		term = apply_pauli_axis(term, 2, permutation[1]);
		term = apply_pauli_axis(term, 1, permutation[0]);
		for (int i = 0; i < 8; ++i) result[i] += static_cast<double>(permutation[3]) * term[i];
	}
	return result;
}

struct spin_bilinears {
	std::array<complex, 9> sigma23{};
	std::array<complex, 27> sigma231{};
	std::array<complex, 9> sigma12{};
	std::array<complex, 9> sigma13{};
	complex identity{0.0, 0.0};
};

std::array<int, 18> bilinear_key(const channel& bra, const channel& ket, int order)
{
	const auto bra_key = basis_key(bra, order);
	const auto ket_key = basis_key(ket, order);
	std::array<int, 18> result{};
	std::copy(bra_key.begin(), bra_key.end(), result.begin());
	std::copy(ket_key.begin(), ket_key.end(), result.begin() + 9);
	return result;
}

std::shared_ptr<const std::vector<spin_bilinears>> get_spin_bilinears(
	const channel& bra, const channel& ket,
	const std::shared_ptr<const angular_grid>& grid)
{
	using key_type = std::array<int, 18>;
	static std::mutex cache_mutex;
	static std::map<key_type, std::shared_ptr<const std::vector<spin_bilinears>>> cache;
	const key_type key = bilinear_key(bra, ket, grid->order);
	{
		std::lock_guard<std::mutex> lock(cache_mutex);
		const auto found = cache.find(key);
		if (found != cache.end()) return found->second;
	}

	const auto bra_basis = get_angular_basis(bra, grid);
	const auto ket_basis = get_angular_basis(ket, grid);
	auto candidate = std::make_shared<std::vector<spin_bilinears>>();
	const std::size_t number_of_angles = static_cast<std::size_t>(grid->order)
	                                   * grid->order * grid->order * grid->order * grid->order;
	candidate->resize(number_of_angles);
	std::size_t angular_index = 0;
	for (int iq = 0; iq < grid->order; ++iq) {
		for (int ipp = 0; ipp < grid->order; ++ipp) {
			for (int iphip = 0; iphip < grid->order; ++iphip) {
				for (int iqp = 0; iqp < grid->order; ++iqp) {
					for (int iphiqp = 0; iphiqp < grid->order; ++iphiqp, ++angular_index) {
						spin_bilinears& values = (*candidate)[angular_index];
						for (int im = 0; im < bra_basis->num_m; ++im) {
							const state8& bra_spin = bra_basis->bra(ipp, iphip, iqp, iphiqp, im);
							const state8& ket_spin = ket_basis->ket(iq, im);
							values.identity += inner_product(bra_spin, ket_spin);
							for (int a = 0; a < 3; ++a) {
								for (int b = 0; b < 3; ++b) {
									state8 operated = apply_pauli_axis(ket_spin, 3, b);
									operated = apply_pauli_axis(operated, 2, a);
									values.sigma23[3 * a + b] += inner_product(bra_spin, operated);

									operated = apply_pauli_axis(ket_spin, 2, b);
									operated = apply_pauli_axis(operated, 1, a);
									values.sigma12[3 * a + b] += inner_product(bra_spin, operated);

									operated = apply_pauli_axis(ket_spin, 3, b);
									operated = apply_pauli_axis(operated, 1, a);
									values.sigma13[3 * a + b] += inner_product(bra_spin, operated);

									for (int c = 0; c < 3; ++c) {
										operated = apply_pauli_axis(ket_spin, 1, c);
										operated = apply_pauli_axis(operated, 3, b);
										operated = apply_pauli_axis(operated, 2, a);
										values.sigma231[9 * a + 3 * b + c] +=
											inner_product(bra_spin, operated);
									}
								}
							}
						}
						const double inverse_m_count = 1.0 / bra_basis->num_m;
						values.identity *= inverse_m_count;
						for (complex& value : values.sigma23) value *= inverse_m_count;
						for (complex& value : values.sigma231) value *= inverse_m_count;
						for (complex& value : values.sigma12) value *= inverse_m_count;
						for (complex& value : values.sigma13) value *= inverse_m_count;
					}
				}
			}
		}
	}

	std::lock_guard<std::mutex> lock(cache_mutex);
	const auto inserted = cache.emplace(key, candidate);
	return inserted.first->second;
}

complex contract_rank2(const std::array<complex, 9>& tensor,
	                     const vector3& first, const vector3& second)
{
	complex result{0.0, 0.0};
	for (int a = 0; a < 3; ++a) {
		for (int b = 0; b < 3; ++b) result += first[a] * second[b] * tensor[3 * a + b];
	}
	return result;
}

complex contract_rank3(const std::array<complex, 27>& tensor,
	                     const vector3& first, const vector3& second,
	                     const vector3& third)
{
	complex result{0.0, 0.0};
	for (int a = 0; a < 3; ++a) {
		for (int b = 0; b < 3; ++b) {
			for (int c = 0; c < 3; ++c) {
				result += first[a] * second[b] * third[c] * tensor[9 * a + 3 * b + c];
			}
		}
	}
	return result;
}

channel make_channel(int alpha, const pw_3N_statespace& pw)
{
	return {
		pw.L_2N_array[alpha], pw.S_2N_array[alpha], pw.J_2N_array[alpha],
		pw.L_1N_array[alpha], pw.two_J_1N_array[alpha],
		pw.two_J_3N_array[alpha], pw.T_2N_array[alpha],
		pw.two_T_3N_array[alpha]
	};
}

double regulator(double p, double q, double cutoff)
{
	const double ratio = (p * p + 0.75 * q * q) / (cutoff * cutoff);
	return std::exp(-(ratio * ratio));
}

} // namespace

chiral_N2LO_3NF_full_reference::chiral_N2LO_3NF_full_reference(
	double c_D, double c_E, double Lambda_3NF_MeV,
	double c1_gev, double c3_gev, double c4_gev, int angular_order)
	: m_c_D(c_D), m_c_E(c_E), m_Lambda(Lambda_3NF_MeV / hbarc),
	  m_c1_gev(c1_gev), m_c3_gev(c3_gev), m_c4_gev(c4_gev),
	  m_c1_fm(c1_gev * hbarc / 1000.0),
	  m_c3_fm(c3_gev * hbarc / 1000.0),
	  m_c4_fm(c4_gev * hbarc / 1000.0),
	  m_fpi(fpi / hbarc), m_mpi(mpi / hbarc),
	  m_lambda_chi(700.0 / hbarc), m_angular_order(angular_order)
{
	if (!(m_Lambda > 0.0)) throw std::invalid_argument("Lambda_3NF must be positive");
	if (m_angular_order < 1) throw std::invalid_argument("3NF angular order must be positive");
}

bool chiral_N2LO_3NF_full_reference::enabled() const
{
	return m_c_D != 0.0 || m_c_E != 0.0 || m_c1_fm != 0.0
	    || m_c3_fm != 0.0 || m_c4_fm != 0.0;
}

double chiral_N2LO_3NF_full_reference::axial_coupling_3nf() const { return gA; }
double chiral_N2LO_3NF_full_reference::pion_decay_constant_mev_3nf() const { return fpi; }
double chiral_N2LO_3NF_full_reference::pion_mass_mev_3nf() const { return mpi; }
double chiral_N2LO_3NF_full_reference::hbarc_mev_fm_3nf() const { return hbarc; }

double chiral_N2LO_3NF_full_reference::W1_element(
	int alpha_r, int alpha_c, double p_r, double q_r,
	double p_c, double q_c, const pw_3N_statespace& pw) const
{
	if (pw.two_J_3N_array[alpha_r] != pw.two_J_3N_array[alpha_c]) return 0.0;
	if (pw.two_T_3N_array[alpha_r] != pw.two_T_3N_array[alpha_c]) return 0.0;
	if (pw.P_3N_array[alpha_r] != pw.P_3N_array[alpha_c]) return 0.0;
	if (!enabled()) return 0.0;

	const channel bra = make_channel(alpha_r, pw);
	const channel ket = make_channel(alpha_c, pw);
	const int order = m_angular_order;
	const auto grid = get_angular_grid(order);
	const auto spin_cache = get_spin_bilinears(bra, ket, grid);
	const auto& cosine = grid->cosine;
	const auto& cosine_weight = grid->cosine_weight;
	const auto& phi = grid->phi;
	const auto& phi_weight = grid->phi_weight;

	const int two_m_t = bra.two_total_T;
	const state8 iso_bra = coupled_three_half_state(bra.t_pair, bra.two_total_T, two_m_t);
	const state8 iso_ket = coupled_three_half_state(ket.t_pair, ket.two_total_T, two_m_t);
	const complex iso23 = inner_product(iso_bra, apply_pair_dot(iso_ket, 2, 3));
	const complex iso_cross = inner_product(iso_bra, apply_triple_cross_123(iso_ket));
	const complex iso12 = inner_product(iso_bra, apply_pair_dot(iso_ket, 1, 2));
	const complex iso13 = inner_product(iso_bra, apply_pair_dot(iso_ket, 1, 3));

	std::array<complex, 5> totals{};
	totals.fill(complex{0.0, 0.0});
	const vector3 p_in{0.0, 0.0, p_c};
	const double m_pi_squared = m_mpi * m_mpi;
	const double f_pi_squared = m_fpi * m_fpi;
	const double f_pi_fourth = f_pi_squared * f_pi_squared;
	const double d_lec = m_c_D / (f_pi_squared * m_lambda_chi);
	const double e_lec = m_c_E / (f_pi_fourth * m_lambda_chi);
	std::size_t angular_index = 0;

	for (int iq = 0; iq < order; ++iq) {
		const double sq = std::sqrt(std::max(0.0, 1.0 - cosine[iq] * cosine[iq]));
		const vector3 q_in{q_c * sq, 0.0, q_c * cosine[iq]};
		for (int ipp = 0; ipp < order; ++ipp) {
			const double spp = std::sqrt(std::max(0.0, 1.0 - cosine[ipp] * cosine[ipp]));
			for (int iphip = 0; iphip < order; ++iphip) {
				const vector3 p_out{
					p_r * spp * std::cos(phi[iphip]),
					p_r * spp * std::sin(phi[iphip]),
					p_r * cosine[ipp]
				};
				for (int iqp = 0; iqp < order; ++iqp) {
					const double sqp = std::sqrt(std::max(0.0, 1.0 - cosine[iqp] * cosine[iqp]));
					for (int iphiqp = 0; iphiqp < order; ++iphiqp, ++angular_index) {
						const vector3 q_out{
							q_r * sqp * std::cos(phi[iphiqp]),
							q_r * sqp * std::sin(phi[iphiqp]),
							q_r * cosine[iqp]
						};
						const vector3 k2_out = subtract(p_out, scale(0.5, q_out));
						const vector3 k3_out = subtract(scale(-1.0, p_out), scale(0.5, q_out));
						const vector3 k2_in = subtract(p_in, scale(0.5, q_in));
						const vector3 k3_in = subtract(scale(-1.0, p_in), scale(0.5, q_in));
						const vector3 transfer1 = subtract(q_out, q_in);
						const vector3 transfer2 = subtract(k2_out, k2_in);
						const vector3 transfer3 = subtract(k3_out, k3_in);
						const double d1 = dot(transfer1, transfer1) + m_pi_squared;
						const double d2 = dot(transfer2, transfer2) + m_pi_squared;
						const double d3 = dot(transfer3, transfer3) + m_pi_squared;
						const double common = gA * gA / (4.0 * f_pi_fourth * d2 * d3);
						const double q2q3 = dot(transfer2, transfer3);

						const spin_bilinears& spin = (*spin_cache)[angular_index];
						const complex spin23 = contract_rank2(spin.sigma23, transfer2, transfer3);
						const complex spin4 = contract_rank3(
							spin.sigma231, transfer2, transfer3, cross(transfer2, transfer3));
						const complex spin_d2 = contract_rank2(spin.sigma12, transfer1, transfer1);
						const complex spin_d3 = contract_rank2(spin.sigma13, transfer1, transfer1);

						const double weight = cosine_weight[iq] * cosine_weight[ipp]
						                    * phi_weight[iphip] * cosine_weight[iqp]
						                    * phi_weight[iphiqp];
						totals[0] += weight * common * (-4.0 * m_c1_fm * m_pi_squared)
						           * iso23 * spin23;
						totals[1] += weight * common * (2.0 * m_c3_fm * q2q3)
						           * iso23 * spin23;
						totals[2] += weight * common * m_c4_fm * iso_cross * spin4;
						totals[3] += weight * (-gA * d_lec / (8.0 * f_pi_squared * d1))
						           * (iso12 * spin_d2 + iso13 * spin_d3);
						totals[4] += weight * e_lec * iso23 * spin.identity;
					}
				}
			}
		}
	}

	const double rotational_volume = 8.0 * pi_value * pi_value;
	const double fourier_normalization = std::pow(2.0 * pi_value, -6);
	const double regulator_product = regulator(p_r, q_r, m_Lambda)
	                               * regulator(p_c, q_c, m_Lambda);
	complex result{0.0, 0.0};
	for (const complex value : totals) result += value;
	result *= rotational_volume * fourier_normalization * regulator_product;
	const double imaginary_tolerance = 2.0e-10 * std::max(1.0, std::abs(result.real()));
	if (std::abs(result.imag()) > imaginary_tolerance) {
		std::fprintf(stderr,
			"complete 3NF reference PWD produced non-real scalar: Re=%.17e Im=%.17e\n",
			result.real(), result.imag());
		throw std::runtime_error("complete 3NF reference PWD failed reality check");
	}
	return result.real();
}
