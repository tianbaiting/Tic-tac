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
