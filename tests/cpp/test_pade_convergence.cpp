#include "core/faddeev_solver/pade_approximant.h"

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

namespace {

int failures = 0;

void check(bool condition, const char* label){
	if (condition){
		std::printf("PASS: %s\n", label);
	} else {
		std::printf("FAIL: %s\n", label);
		++failures;
	}
}

} // namespace

int main(){
	constexpr std::size_t max_order = 14;

	// A nontrivial [24/24] rational series checks the numerically stable
	// coefficient-solve implementation.  The former determinant ratio
	// overflowed or underflowed on similarly scaled production histories.
	constexpr std::size_t rational_order = 24;
	std::vector<cdouble> denominator(1, cdouble(1.0, 0.0));
	for (std::size_t root = 0; root < rational_order; ++root){
		const double factor = 0.02 + 0.006 * static_cast<double>(root);
		std::vector<cdouble> next(denominator.size() + 1, cdouble(0.0, 0.0));
		for (std::size_t k = 0; k < denominator.size(); ++k){
			next[k] += denominator[k];
			next[k + 1] -= factor * denominator[k];
		}
		denominator = next;
	}
	std::vector<cdouble> numerator(rational_order + 1);
	for (std::size_t k = 0; k <= rational_order; ++k){
		numerator[k] = cdouble(
			std::sin(0.7 * static_cast<double>(k)) / static_cast<double>(k + 1),
			std::cos(0.4 * static_cast<double>(k)) / static_cast<double>(2 * k + 1));
	}
	std::vector<cdouble> rational_series(2 * rational_order + 1);
	for (std::size_t n = 0; n < rational_series.size(); ++n){
		rational_series[n] = n <= rational_order ? numerator[n] : cdouble(0.0, 0.0);
		for (std::size_t j = 1; j <= std::min(n, rational_order); ++j){
			rational_series[n] -= denominator[j] * rational_series[n - j];
		}
	}
	cdouble exact_rational = cdouble(0.0, 0.0);
	for (const cdouble value : numerator){
		exact_rational += value;
	}
	cdouble exact_denominator = cdouble(0.0, 0.0);
	for (const cdouble value : denominator){
		exact_denominator += value;
	}
	exact_rational /= exact_denominator;
	const cdouble computed_rational = pade_approximant(
		rational_series.data(), rational_order, rational_order, cdouble(1.0, 0.0));
	check(std::abs(computed_rational - exact_rational) < 1e-10,
	      "[24/24] coefficient solve reproduces a known rational function");

	// Rank-deficient Padé systems occur when the true denominator degree is
	// smaller than the requested table order.  They have a valid generalized
	// approximant and must not become determinant-form 0/0.
	std::vector<cdouble> geometric(13);
	for (std::size_t n = 0; n < geometric.size(); ++n){
		geometric[n] = std::pow(cdouble(0.35, 0.0), static_cast<int>(n));
	}
	const cdouble geometric_pade = pade_approximant(
		geometric.data(), 6, 6, cdouble(1.0, 0.0));
	check(std::isfinite(geometric_pade.real()) && std::isfinite(geometric_pade.imag())
	      && std::abs(geometric_pade - cdouble(1.0 / 0.65, 0.0)) < 1e-11,
	      "rank-deficient Padé table falls back to a finite generalized solution");

	// This is the failure mode of the removed legacy criterion: because every
	// successive difference is > 1, its "best" index stayed at zero and P[0/0]
	// was declared genuinely converged at order five.
	std::vector<cdouble> moving(max_order + 1);
	for (std::size_t n = 0; n <= max_order; ++n){
		moving[n] = cdouble(2.0 * n, -0.5 * n);
	}
	auto early = assess_pade_convergence(moving.data(), 5, max_order);
	check(!early.stop, "moving Padé prefix is not frozen at order five");
	check(!early.genuinely_converged, "moving Padé prefix is not called converged");

	auto moving_final = assess_pade_convergence(moving.data(), max_order, max_order);
	check(moving_final.stop, "maximum order terminates the sequence");
	check(moving_final.max_order_truncated,
	      "moving final tail is labelled max-order truncated");
	check(moving_final.selected_order == 1,
	      "closest finite pair is used only as a labelled optimal truncation");

	std::vector<cdouble> stable(max_order + 1);
	for (std::size_t n = 0; n <= max_order; ++n){
		stable[n] = cdouble(1.25, -0.40);
	}
	stable[11] += cdouble(2e-6, -1e-6);
	stable[12] += cdouble(1e-6, -5e-7);
	stable[13] += cdouble(4e-7, -2e-7);
	stable[14] += cdouble(1e-7, -5e-8);
	auto stable_final = assess_pade_convergence(stable.data(), max_order, max_order);
	check(stable_final.genuinely_converged,
	      "three-step final Padé tail is genuinely converged");
	check(!stable_final.max_order_truncated,
	      "stable final tail is not labelled truncated");

	// A coincidental early plateau is not evidence if the final tail moves.
	std::vector<cdouble> false_plateau(max_order + 1, cdouble(0.5, 0.2));
	false_plateau[12] = cdouble(0.8, 0.1);
	false_plateau[13] = cdouble(0.3, 0.7);
	false_plateau[14] = cdouble(0.9, -0.2);
	auto plateau_final = assess_pade_convergence(false_plateau.data(), max_order, max_order);
	check(plateau_final.max_order_truncated,
	      "early plateau followed by drift is not called converged");

	std::vector<cdouble> singular_tail(max_order + 1);
	for (std::size_t n = 0; n < max_order; ++n){
		const double residual = 1.0 / static_cast<double>(n + 1);
		singular_tail[n] = cdouble(2.0 + residual, -0.3 + 0.5 * residual);
	}
	singular_tail[max_order] = cdouble(std::numeric_limits<double>::quiet_NaN(), 0.0);
	auto singular_final = assess_pade_convergence(singular_tail.data(), max_order, max_order);
	check(singular_final.selected_order == max_order - 1,
	      "non-finite final approximant falls back to latest finite order");
	check(singular_final.max_order_truncated,
	      "non-finite final approximant cannot be certified converged");

	std::printf("\nPadé convergence policy: %d failure(s)\n", failures);
	return failures == 0 ? 0 : 1;
}
