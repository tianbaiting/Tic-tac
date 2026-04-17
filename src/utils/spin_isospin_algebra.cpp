#include "spin_isospin_algebra.h"
#include <cmath>

double reduced_me_sigma_dot_qhat(int L_prime, int S_prime, int J_2N,
                                  int L, int S) {
    // Only nonzero for S=S'=1 (spin-triplet pair)
    if (S != 1 || S_prime != 1) return 0.0;

    // Selection rule: |L' - L| = 1
    int dL = std::abs(L_prime - L);
    if (dL != 1) return 0.0;

    // Triangle inequality for (L', S'=1, J) and (L, S=1, J)
    if (J_2N < std::abs(L_prime - 1) || J_2N > L_prime + 1) return 0.0;
    if (J_2N < std::abs(L - 1) || J_2N > L + 1) return 0.0;

    // <L' 1 J || sigma.qhat || L 1 J>
    // = (-1)^{L'+1+J} * sqrt(3*(2L'+1)*(2L+1)) * {L' 1 J; 1 L 1} * CG(L,0;1,0|L',0)
    int phase_exp = L_prime + 1 + J_2N;
    double phase = (phase_exp % 2 == 0) ? 1.0 : -1.0;

    // CG(L,0;1,0|L',0) — using twice-integer convention
    double cg = clebsch_gordan(2*L, 2*1, 2*L_prime, 0, 0, 0);

    // 6j symbol {L' 1 J; 1 L 1} — using twice-integer convention
    double w6j = wigner_6j(2*L_prime, 2*1, 2*J_2N, 2*1, 2*L, 2*1);

    return phase * std::sqrt(3.0 * (2*L_prime + 1.0) * (2*L + 1.0)) * w6j * cg;
}

double reduced_me_sigma1_dot_sigma3(int L_prime, int S_prime, int J_2N_prime,
                                     int l_prime, int two_j_1N_prime, int two_J_3N,
                                     int L, int S, int J_2N,
                                     int l, int two_j_1N) {
    // Placeholder — full implementation in Task 4
    return 0.0;
}

double reduced_me_sigma_q_tensor(int L_prime, int S_prime, int J_2N,
                                  int L, int S, int rank) {
    // Placeholder — full implementation in Task 6
    return 0.0;
}
