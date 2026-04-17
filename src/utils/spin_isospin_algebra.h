#ifndef SPIN_ISOSPIN_ALGEBRA_H
#define SPIN_ISOSPIN_ALGEBRA_H

#include "coupling_coefficients.h"

// [EN] Partial-wave reduced matrix elements for spin-momentum and isospin operators
// appearing in the N2LO chiral 3NF. All angular momenta use integer conventions
// (NOT twice-integers) except where noted. The coupling_coefficients functions
// (CG, 6j, 9j) internally use twice-integers — conversions happen inside.
//
// Reference: Epelbaum, Glöckle, Meissner, PRC 66 (2002) 064001, Appendix A.
//
// / [CN] N2LO 手征 3NF 中自旋-动量和同位旋算符的分波约化矩阵元。

// tau_2 . tau_3 eigenvalue in the pair isospin basis.
// Returns 2*T*(T+1) - 3:  T=0 → -3,  T=1 → +1.
inline double tau23_eigenvalue(int T_2N) {
    return 2.0 * T_2N * (T_2N + 1.0) - 3.0;
}

// Reduced matrix element of sigma.q_hat between pair partial-wave states.
// <L' S' J_2N || sigma . q_hat || L S J_2N>
// This is nonzero only when S=S'=1, |L'-L|=1, and J_2N is valid for both (L,S) and (L',S').
// q_hat is the unit vector along q; the |q| dependence is factored out.
//
// Formula (for S=S'=1, |L'-L|=1):
//   (-1)^{L'+1+J} * sqrt(3*(2L'+1)*(2L+1)) * {L' 1 J; 1 L 1} * CG(L,0;1,0|L',0)
double reduced_me_sigma_dot_qhat(int L_prime, int S_prime, int J_2N,
                                  int L, int S);

// sigma_1.sigma_3 reduced matrix element in the recoupled spectator-pair spin basis.
// Placeholder — full implementation in Task 4.
double reduced_me_sigma1_dot_sigma3(int L_prime, int S_prime, int J_2N_prime,
                                     int l_prime, int two_j_1N_prime, int two_J_3N,
                                     int L, int S, int J_2N,
                                     int l, int two_j_1N);

// Rank-2 tensor operator decomposition.
// Placeholder — full implementation in Task 6.
double reduced_me_sigma_q_tensor(int L_prime, int S_prime, int J_2N,
                                  int L, int S, int rank);

#endif // SPIN_ISOSPIN_ALGEBRA_H
