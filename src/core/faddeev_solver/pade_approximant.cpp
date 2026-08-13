#include "pade_approximant.h"

#include <lapacke.h>

#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>

cdouble pade_approximant(const cdouble* a_coeff_array,
                         std::size_t N,
                         std::size_t M,
                         cdouble z){
	// Solve directly for Q(z)=1+q_1*z+...+q_M*z^M from
	//   sum_{j=1}^M a_{N+i-j} q_j = -a_{N+i},  i=1,...,M.
	// This is algebraically equivalent to the former det(P)/det(Q) formula,
	// but avoids forming two determinants whose products overflow/underflow at
	// useful diagonal orders (the reduced 3NF diagnostic already fails by
	// [24/24] despite a finite, stable rational approximant).
	std::vector<cdouble> denominator_coefficients(M + 1, cdouble(0.0, 0.0));
	denominator_coefficients[0] = cdouble(1.0, 0.0);

	if (M > 0){
		std::vector<cdouble> system(M * M);
		std::vector<cdouble> rhs(M);
		for (std::size_t row = 0; row < M; ++row){
			const std::size_t i = row + 1;
			rhs[row] = -a_coeff_array[N + i];
			double row_scale = std::abs(rhs[row]);
			for (std::size_t col = 0; col < M; ++col){
				const std::size_t j = col + 1;
				const cdouble value = a_coeff_array[N + i - j];
				system[row * M + col] = value;
				row_scale = std::max(row_scale, std::abs(value));
			}
			if (row_scale > 0.0 && std::isfinite(row_scale)){
				rhs[row] /= row_scale;
				for (std::size_t col = 0; col < M; ++col){
					system[row * M + col] /= row_scale;
				}
			}
		}

		// Keep pristine scaled inputs for the rank-revealing fallback because
		// zgesv overwrites both arrays.
		const std::vector<cdouble> system_copy = system;
		const std::vector<cdouble> rhs_copy = rhs;
		std::vector<lapack_int> pivots(M);
		lapack_int info = LAPACKE_zgesv(
			LAPACK_ROW_MAJOR, static_cast<lapack_int>(M), 1,
			reinterpret_cast<lapack_complex_double*>(system.data()),
			static_cast<lapack_int>(M), pivots.data(),
			reinterpret_cast<lapack_complex_double*>(rhs.data()), 1);

		if (info > 0){
			// A defective Padé table can have a rank-deficient coefficient matrix
			// even when a lower-degree rational representation is well defined.
			// Rank-revealing QR supplies that representation instead of returning
			// the determinant form's indeterminate 0/0.
			system = system_copy;
			rhs = rhs_copy;
			std::vector<lapack_int> column_pivots(M, 0);
			lapack_int rank = 0;
			info = LAPACKE_zgelsy(
				LAPACK_ROW_MAJOR, static_cast<lapack_int>(M),
				static_cast<lapack_int>(M), 1,
				reinterpret_cast<lapack_complex_double*>(system.data()),
				static_cast<lapack_int>(M),
				reinterpret_cast<lapack_complex_double*>(rhs.data()), 1,
				column_pivots.data(), -1.0, &rank);
		}
		if (info != 0){
			const double nan = std::numeric_limits<double>::quiet_NaN();
			return cdouble(nan, nan);
		}
		for (std::size_t j = 1; j <= M; ++j){
			denominator_coefficients[j] = rhs[j - 1];
		}
	}

	std::vector<cdouble> numerator_coefficients(N + 1, cdouble(0.0, 0.0));
	for (std::size_t k = 0; k <= N; ++k){
		for (std::size_t j = 0; j <= std::min(k, M); ++j){
			numerator_coefficients[k] +=
				denominator_coefficients[j] * a_coeff_array[k - j];
		}
	}

	// Horner evaluation avoids the extra dynamic range of explicit powers.
	cdouble numerator = numerator_coefficients[N];
	for (std::size_t k = N; k-- > 0;){
		numerator = numerator * z + numerator_coefficients[k];
	}
	cdouble denominator = denominator_coefficients[M];
	for (std::size_t j = M; j-- > 0;){
		denominator = denominator * z + denominator_coefficients[j];
	}
	return numerator / denominator;
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
