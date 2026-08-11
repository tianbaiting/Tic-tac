#ifndef PADE_APPROXIMANT_H
#define PADE_APPROXIMANT_H

#include <cstddef>

#include "type_defs.h"

// Evaluate the [N/M] Padé approximant at z from the scalar series
// a_0 + a_1 z + ... + a_(N+M) z^(N+M).
cdouble pade_approximant(const cdouble* a_coeff_array,
                         std::size_t N,
                         std::size_t M,
                         cdouble z);

// Honest convergence decision for a scalar diagonal Padé sequence P[n/n].
// The solver evaluates all configured orders before accepting convergence:
// an old "best order has not changed recently" heuristic could select P[0/0]
// even while the later approximants were still moving by O(1).  A result is
// now called converged only when the final consecutive tail is numerically
// stable.  Otherwise the upper member of the closest finite consecutive pair
// in the complete history is returned as a best-effort optimal truncation,
// explicitly labelled max-order truncated by the caller.
struct pade_convergence_decision {
	std::size_t selected_order = 0;
	bool stop = false;
	bool genuinely_converged = false;
	bool max_order_truncated = false;
};

pade_convergence_decision assess_pade_convergence(
	const cdouble* approximants,
	std::size_t current_order,
	std::size_t max_order,
	std::size_t stable_steps = 3,
	double relative_tolerance = 1e-5,
	double absolute_tolerance = 1e-7);

#endif // PADE_APPROXIMANT_H
