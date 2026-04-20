#ifndef BUILD_PW_3N_STATESPACE_H
#define BUILD_PW_3N_STATESPACE_H

#include "type_defs.h"
#include "read_triton_psi.h"

// [EN] Build a pw_3N_statespace whose α-table is a direct copy of ψ's α-table.
// Fills two_J_3N_array, two_T_3N_array, P_3N_array from the ψ quantum numbers
// (2J_3N=1 for triton, 2T_3N=1, parity = (-1)^(L_2N+L_1N)). Caller owns the arrays
// allocated inside; free with free_pw_3N_statespace_triton.
pw_3N_statespace build_pw_3n_statespace_from_triton(const triton_wavefunction& w);
void free_pw_3n_statespace_triton(pw_3N_statespace& s);

#endif
