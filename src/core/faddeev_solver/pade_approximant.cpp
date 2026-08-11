#include "pade_approximant.h"

#include "utils/matrix_routines.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

cdouble pade_approximant(const cdouble* a_coeff_array,
                         std::size_t N,
                         std::size_t M,
                         cdouble z){
	const std::size_t dim = M + 1;
	std::vector<cdouble> P_array(dim * dim);
	std::vector<cdouble> Q_array(dim * dim);

	for (std::size_t row_idx = 0; row_idx < M; row_idx++){
		for (std::size_t col_idx = 0; col_idx < dim; col_idx++){
			const cdouble value = a_coeff_array[N - M + 1 + row_idx + col_idx];
			P_array[row_idx * dim + col_idx] = value;
			Q_array[row_idx * dim + col_idx] = value;
		}
	}

	for (std::size_t col_idx = 0; col_idx < dim; col_idx++){
		Q_array[M * dim + col_idx] = std::pow(z, M - col_idx);
		P_array[M * dim + col_idx] = 0.0;
		for (std::size_t j = M - col_idx; j < N + 1; j++){
			P_array[M * dim + col_idx] +=
				a_coeff_array[j - (M - col_idx)] * std::pow(z, j);
		}
	}

	const cdouble P_det = determinant(P_array.data(), (int)dim);
	const cdouble Q_det = determinant(Q_array.data(), (int)dim);
	return P_det / Q_det;
}

pade_convergence_decision assess_pade_convergence(
	const cdouble* approximants,
	std::size_t current_order,
	std::size_t max_order,
	std::size_t stable_steps,
	double relative_tolerance,
	double absolute_tolerance){
	pade_convergence_decision decision;
	decision.selected_order = current_order;

	// Do not compact a row or freeze an element from a short, accidentally flat
	// prefix.  The full configured Padé history is cheap compared with building
	// the physical kernel and is needed for an honest convergence statement.
	if (current_order < max_order){
		return decision;
	}
	decision.stop = true;

	// Identify the latest finite value and the closest finite consecutive pair
	// over the COMPLETE history.  The latter is an optimal-truncation heuristic,
	// never a convergence certificate: the final-tail test below alone controls
	// genuinely_converged.
	bool selected_is_final = false;
	bool have_finite_value = false;
	std::size_t latest_finite_order = current_order;
	for (std::size_t offset = 0; offset <= current_order; ++offset){
		const std::size_t order = current_order - offset;
		const cdouble value = approximants[order];
		if (std::isfinite(value.real()) && std::isfinite(value.imag())){
			latest_finite_order = order;
			have_finite_value = true;
			selected_is_final = (order == current_order);
			break;
		}
	}

	bool have_finite_pair = false;
	std::size_t closest_pair_order = latest_finite_order;
	double closest_pair_difference = std::numeric_limits<double>::infinity();
	for (std::size_t upper = 1; upper <= current_order; ++upper){
		const cdouble current = approximants[upper];
		const cdouble previous = approximants[upper - 1];
		if (!std::isfinite(current.real()) || !std::isfinite(current.imag())
		 || !std::isfinite(previous.real()) || !std::isfinite(previous.imag())){
			continue;
		}
		const double difference = std::abs(current - previous);
		if (difference < closest_pair_difference){
			closest_pair_difference = difference;
			closest_pair_order = upper;
			have_finite_pair = true;
		}
	}

	bool stable_tail = selected_is_final && stable_steps > 0
	                && current_order >= stable_steps;
	if (stable_tail){
		for (std::size_t step = 0; step < stable_steps; ++step){
			const std::size_t upper = current_order - step;
			const cdouble current = approximants[upper];
			const cdouble previous = approximants[upper - 1];
			if (!std::isfinite(current.real()) || !std::isfinite(current.imag())
			 || !std::isfinite(previous.real()) || !std::isfinite(previous.imag())){
				stable_tail = false;
				break;
			}
			const double scale = std::max(std::abs(current), std::abs(previous));
			const double tolerance = absolute_tolerance + relative_tolerance * scale;
			if (std::abs(current - previous) > tolerance){
				stable_tail = false;
				break;
			}
		}
	}

	decision.genuinely_converged = stable_tail;
	decision.max_order_truncated = !stable_tail;
	if (stable_tail){
		decision.selected_order = current_order;
	} else if (have_finite_pair){
		decision.selected_order = closest_pair_order;
	} else if (have_finite_value){
		decision.selected_order = latest_finite_order;
	}
	return decision;
}
