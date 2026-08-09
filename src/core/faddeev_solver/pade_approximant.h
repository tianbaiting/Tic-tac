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

#endif // PADE_APPROXIMANT_H
