#ifndef CHIRAL_3NF_RECOUPLING_H
#define CHIRAL_3NF_RECOUPLING_H

// [EN] Spin-isospin-spatial recoupling coefficients for the chiral N2LO
// three-nucleon force partial-wave matrix elements. These coefficients
// are analytical, channel-dependent, and momentum-independent, so they
// factor out of the radial hot loop in W^(1) evaluation.
//
// Reference formulas: tools/check_3nf_normalization/formula_reference.md
//   sections 2, 3, 4 (derived from Epelbaum, Nogga, Glöckle et al.
//   PRC 66 (2002) 064001 App. A-1/A-4, Golak et al. EPJA 43 (2010) 241,
//   Hebeler et al. PRC 91 (2015) 044001 eqs. 13-15).
//
// All angular momenta use the integer convention (NOT twice-integers)
// except for half-integer quantities (`two_J_1N`, `two_J_3N`, `two_T_3N`)
// which follow the code-wide `two_*` convention. Internally the helpers
// convert to twice-integers for GSL coupling functions.
//
// / [CN] 手征 N2LO 三核子力分波矩阵元的自旋-同位旋-空间重耦合系数。
// 这些系数是解析的、通道依赖且动量无关的，因此可以从径向热循环中分离出来。

// Rank-0 (scalar pair operator) 3NF recoupling coefficient.
//
// Applies to operators of the form
//   O_scalar = (spatial scalar) * (spin scalar) * (isospin scalar)
// acting on the pair (2,3) with the spectator (1) unchanged. The three
// scalar factors include:
//   - spin:    sigma_2 . sigma_3      (eigenvalue 2 S_2N (S_2N+1) - 3)
//   - isospin: tau_2 . tau_3          (eigenvalue 2 T_2N (T_2N+1) - 3)
//   - spatial: any channel-diagonal scalar (pair orbital unchanged)
//
// The returned coefficient is the Kronecker product of selection rules
// times the pair spin and isospin eigenvalues. Since the operator is
// scalar in every subspace, the recoupling reduces to a product of
// delta-functions times the algebraic eigenvalues (all other recoupling
// 6j's collapse to 1/(2J+1) normalisations absorbed into the convention).
//
// This form matches A_c1_rank0 in formula_reference.md §3.4 and is also
// the limiting case of A_cE in §1.4 and A_c3_rank0 in §4 (up to the
// overall LEC- and momentum-kernel dependent prefactors handled by the
// caller).
//
// Arguments (row = bra "r", column = ket "c"):
//   L_2N_*, S_2N_*, J_2N_*, T_2N_* : pair partial-wave integers
//   L_1N_*, two_J_1N_*             : spectator orbital integer, 2*spectator total (half-integer)
//   two_J_3N, two_T_3N             : 2*3N conserved quantum numbers
//
// Returns: real-valued recoupling coefficient (no regulator, no LEC, no
// momentum kernel).
double recoupling_3nf_scalar(
    int L_2N_r, int S_2N_r, int J_2N_r, int T_2N_r,
    int L_1N_r, int two_J_1N_r, int two_J_3N,
    int L_2N_c, int S_2N_c, int J_2N_c, int T_2N_c,
    int L_1N_c, int two_J_1N_c,
    int two_T_3N);

// Rank-2 pair-operator 3NF recoupling coefficient.
//
// Applies to operators of the form
//   O_rank2 = (isospin scalar) * [T^2(pair spin) . T^2(pair spatial)]
// i.e. a rank-2 spin tensor on the pair contracted with a rank-2 spatial
// tensor on the pair. This is the piece responsible for opening
// L_2N = 0 <-> L_2N' = 2 transitions (e.g. 3S1 <-> 3D1) within a fixed
// J_2N and fixed J_3N.
//
// Formula (formula_reference.md §3.4, Hebeler et al. PRC 91 044001 eq. 15):
//   A_rank2 = (2 T_2N (T_2N+1) - 3)   * delta_{T_2N T_2N'}
//           * delta_{L_1N L_1N'} * delta_{two_J_1N two_J_1N'}
//           * sqrt(30 * (2 S_2N+1)(2 S_2N'+1)(2 J_2N+1)(2 J_2N'+1))
//           * NineJ(L_2N , S_2N , J_2N ; L_2N', S_2N', J_2N'; 2, 2, 0)
//           * CG(L_2N, 0; 2, 0 | L_2N', 0)
//
// The factor sqrt(30) is the reduced matrix element
// <1 || [sigma_2 x sigma_3]_2 || 1> (standard Racah result). The 9j
// reduces the rank-2 on pair-spin x rank-2 on pair-orbital product to
// overall J_2N, and the Clebsch-Gordan carries the spatial rank-2 Y_2
// reduced matrix element (factorised up to the radial integral which
// supplies the remaining sqrt((2L+1)/(4pi(2L'+1))) normalisation in the
// caller).
//
// Selection rules enforced: J_2N = J_2N' (rank-2 scalar in 3N coupling),
// T_2N = T_2N' (pair-diagonal tau2.tau3), l_1N = l_1N', j_1N = j_1N'
// (spectator untouched by the rank-2 pair operator).
//
// Returns: real-valued recoupling coefficient. Returns 0 when any
// selection rule is violated.
double recoupling_3nf_rank2(
    int L_2N_r, int S_2N_r, int J_2N_r, int T_2N_r,
    int L_1N_r, int two_J_1N_r, int two_J_3N,
    int L_2N_c, int S_2N_c, int J_2N_c, int T_2N_c,
    int L_1N_c, int two_J_1N_c,
    int two_T_3N);


// ============================================================================
// c_D 1PE-contact spin-isospin recoupling coefficients.
//
// The c_D operator (E2002 eq. 2.10) has the spin structure
//   (sigma_1 . q_3)(sigma_3 . q_3) / (q_3^2 + mpi^2)
// which decomposes on spherical-tensor grounds (formula_reference.md §2.2):
//   (sigma_1 . q_hat)(sigma_3 . q_hat) = (1/3)(sigma_1.sigma_3)   [rank-0]
//                       + [sigma_1 ⊗ sigma_3]_2 . [q_hat ⊗ q_hat]_2  [rank-2]
//
// Isospin structure: (tau_1 . tau_3) — connects pair isospin T_2N = 0 <-> 1.
//
// Selection rules from E2002 eq. A-1 (spectator-1 picture):
//   L_2N = L_2N' = 0    (contact pair vertex: only S-wave pair contributes)
//   L_1N = L_1N'         (sigma operators do not change spectator orbital)
//
// ============================================================================

// Rank-0 (1/3 sigma_1.sigma_3) x (tau_1.tau_3) recoupling coefficient.
//
// Returns the scalar (rank-0) spin-isospin matrix element for the c_D 1PE-CT
// term in the 3N Jj-coupled partial-wave basis:
//
//   (1/3) * <(L_r,S_r)J_r, (l_r, s_1=1/2)j_r; J_3N | sigma_1.sigma_3 |
//            (L_c,S_c)J_c, (l_c, s_1=1/2)j_c; J_3N >
//   * <(T_r, t_1=1/2)T_3N | tau_1.tau_3 | (T_c, t_1=1/2)T_3N >
//
// Selection rules enforced:
//   L_2N_r = L_2N_c = 0   (contact pair vertex, E2002 eq. A-1)
//   L_1N_r = L_1N_c       (orbital angular momenta preserved by sigma operators)
//
// The (1/3) factor is the rank-0 coefficient in the
//   (sigma_1.q_hat)(sigma_3.q_hat) = (1/3)(sigma_1.sigma_3) + rank2
// decomposition.
//
// Arguments:
//   (as per recoupling_3nf_scalar — same channel labelling, see that function)
//
// Returns: real-valued recoupling coefficient. Returns 0 when any selection
//          rule is violated (including tau1.tau3 = 0 for diagonal T channels).
double recoupling_3nf_1pe_ct_scalar(
    int L_2N_r, int S_2N_r, int J_2N_r, int T_2N_r,
    int L_1N_r, int two_J_1N_r, int two_J_3N,
    int L_2N_c, int S_2N_c, int J_2N_c, int T_2N_c,
    int L_1N_c, int two_J_1N_c,
    int two_T_3N);

// Rank-2 [sigma_1 ⊗ sigma_3]_2 . Y_2(q_hat) recoupling coefficient for c_D.
//
// The rank-2 spatial tensor [q_hat ⊗ q_hat]_2 ~ Y_2(q_hat) acts on the
// SPECTATOR direction (q_hat = q/|q|, the spectator Jacobi momentum).
// After partial-wave projection onto (L_2N, l_1N) states, the Y_2 on the
// spectator line requires:
//   CG(l_1N, 0; 2, 0 | l_1N', 0) != 0
// which by the triangle rule demands |l_1N - l_1N'| <= 2 and
// l_1N + l_1N' + 2 even.  For l_1N = l_1N' = 0 the CG is zero
// (since the triangle rule 0+2 >= 0 fails).
//
// This function returns the full rank-2 recoupling coefficient including
// the 9j and CG from the spectator-line Y_2 projection.  For l_1N = 0
// it returns 0 exactly.  For l_1N >= 2 (where it is first non-zero) it
// returns the full coefficient.
//
// Returns: real-valued recoupling coefficient. Returns 0 when any selection
//          rule is violated.
double recoupling_3nf_1pe_ct_rank2(
    int L_2N_r, int S_2N_r, int J_2N_r, int T_2N_r,
    int L_1N_r, int two_J_1N_r, int two_J_3N,
    int L_2N_c, int S_2N_c, int J_2N_c, int T_2N_c,
    int L_1N_c, int two_J_1N_c,
    int two_T_3N);

#endif // CHIRAL_3NF_RECOUPLING_H
