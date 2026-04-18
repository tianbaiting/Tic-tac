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
    // [EN] sigma_1 . sigma_3 matrix element in the 3N coupled basis:
    //   <((L', (s2 s3)S') J_2N', (l', s1) j_1N') J_3N | sigma_1.sigma_3
    //    | ((L, (s2 s3)S) J_2N, (l, s1) j_1N) J_3N >
    //
    // sigma_1 acts on spectator (s1), sigma_3 on pair particle (s3).
    // Scalar (rank-0) operator in J_3N space built from rank-1 operators
    // in different subspaces, requiring 9j recoupling (Varshalovich 12.1.5.1).
    //
    // Formula:
    //   ME = -sqrt(3) * sqrt(2*J_3N+1) * 9j * rme_pair * rme_spec
    //
    // where:
    //   9j = {J_2N, J_2N', 1; j_1N, j_1N', 1; J_3N, J_3N, 0}  (Varshalovich rows)
    //   rme_pair = <(L,S')J_2N' || sigma_3 || (L,S) J_2N>
    //   rme_spec = <(l,1/2)j_1N' || sigma_1 || (l,1/2) j_1N>
    //
    // / [CN] sigma_1 . sigma_3 算子在 3N 耦合基中的矩阵元，需要 9j 重耦合。

    // Selection rules: sigma doesn't act on orbital angular momenta
    if (L_prime != L) return 0.0;
    if (l_prime != l) return 0.0;

    // --- 6j for pair spin peel (sigma_3 acts on s3 within (s2, s3)S) ---
    // <(1/2,1/2)S' || sigma_3 || (1/2,1/2)S>
    //   = (-1)^{s2+s3'+S+1} hat(S') hat(S) {s3' S' s2; S s3 1} <1/2||sigma||1/2>
    //   = (-1)^{1+S} hat(S') hat(S) {1/2 S' 1/2; S 1/2 1} sqrt(6)
    double w6j_pair_spin = wigner_6j(1, 2*S_prime, 1, 2*S, 1, 2);
    if (std::abs(w6j_pair_spin) < 1e-15) return 0.0;

    // --- 6j for pair orbital peel (sigma_3 acts on S within (L, S)J_2N) ---
    // <(L, S') J_2N' || sigma_3 || (L, S) J_2N>
    //   = (-1)^{L+S'+J_2N+1} hat(J_2N') hat(J_2N) {S' J_2N' L; J_2N S 1} <S'||sigma3||S>
    double w6j_pair_orbital = wigner_6j(2*S_prime, 2*J_2N_prime, 2*L,
                                         2*J_2N, 2*S, 2);
    if (std::abs(w6j_pair_orbital) < 1e-15) return 0.0;

    // --- 6j for spectator peel (sigma_1 acts on s1 within (l, s1)j_1N) ---
    // <(l, 1/2) j_1N' || sigma_1 || (l, 1/2) j_1N>
    //   = (-1)^{l+1/2+j_1N+1} hat(j_1N') hat(j_1N) {1/2 j_1N' l; j_1N 1/2 1} sqrt(6)
    double w6j_spec = wigner_6j(1, two_j_1N_prime, 2*l,
                                 two_j_1N, 1, 2);
    if (std::abs(w6j_spec) < 1e-15) return 0.0;

    // --- 9j for recoupling (J_2N, j_1N) J_3N → Varshalovich convention ---
    // Rows: (J_2N, J_2N', 1), (j_1N, j_1N', 1), (J_3N, J_3N, 0)
    double w9j = wigner_9j(2*J_2N, 2*J_2N_prime, 2,
                            two_j_1N, two_j_1N_prime, 2,
                            two_J_3N, two_J_3N, 0);
    if (std::abs(w9j) < 1e-15) return 0.0;

    // --- Phase from peel operations ---
    // pair orbital: (-1)^{L + S' + J_2N + 1}
    // pair spin:    (-1)^{1/2 + 1/2 + S + 1} = (-1)^{S}   [since (-1)^2 = 1]
    // spectator:    (-1)^{l + 1/2 + j_1N + 1}
    //
    // Total = (-1)^{L + S' + J_2N + 1 + S + l + 1/2 + j_1N + 1}
    //       = (-1)^{L + S' + J_2N + S + l + 2 + 1/2 + j_1N}
    //
    // j_1N = two_j_1N/2 half-integer, so 1/2 + j_1N = (1 + two_j_1N)/2 is integer.
    // Use twice-exponent to avoid floating point:
    // 2*exp = 2L + 2S' + 2*J_2N + 2S + 2l + 4 + 1 + two_j_1N
    //       = 2L + 2S' + 2*J_2N + 2S + 2l + 5 + two_j_1N
    int two_exp = 2*L + 2*S_prime + 2*J_2N + 2*S + 2*l + 5 + two_j_1N;
    double phase = ((two_exp / 2) % 2 == 0) ? 1.0 : -1.0;

    // --- Prefactor ---
    // -sqrt(3) from sigma1.sigma3 = -sqrt(3) [sigma1 x sigma3]^0
    // sqrt(2*J_3N+1) from reduced ME → matrix element conversion
    // 6 from <1/2||sigma||1/2>^2 = sqrt(6)^2
    double prefactor = -std::sqrt(3.0) * std::sqrt(two_J_3N + 1.0) * 6.0;

    // --- Angular momentum hat factors ---
    double hat_J2N = std::sqrt((2*J_2N_prime + 1.0) * (2*J_2N + 1.0));
    double hat_S   = std::sqrt((2*S_prime + 1.0) * (2*S + 1.0));
    double hat_j1N = std::sqrt((two_j_1N_prime + 1.0) * (two_j_1N + 1.0));

    return prefactor * phase * hat_J2N * hat_S * hat_j1N
           * w9j * w6j_pair_orbital * w6j_pair_spin * w6j_spec;
}

double tau1_dot_tau3(int T_prime, int T, int two_T_3N) {
    // [EN] tau_1 . tau_3 matrix element in the pair isospin basis:
    //   <((1/2,1/2)T', 1/2) T_3N | tau_1.tau_3 | ((1/2,1/2)T, 1/2) T_3N>
    //
    // Same algebraic structure as sigma_1.sigma_3 but without orbital angular
    // momenta. The coupling is ((tau_2, tau_3)T, tau_1=1/2) T_3N.
    //
    // Formula (Varshalovich 12.1.5.1):
    //   ME = -sqrt(3) * sqrt(2*T_3N+1) * 9j * rme_pair * rme_spec
    //   where rme_pair = (-1)^{1+T} hat(T') hat(T) {1/2 T' 1/2; T 1/2 1} sqrt(6)
    //   and   rme_spec = sqrt(6)  (spectator tau_1=1/2, no orbital)
    //
    // / [CN] tau_1 . tau_3 同位旋算子矩阵元，与 sigma_1.sigma_3 结构相同但无轨道角动量。

    // 6j for pair isospin peel: {1/2 T' 1/2; T 1/2 1}
    double w6j_pair = wigner_6j(1, 2*T_prime, 1, 2*T, 1, 2);
    if (std::abs(w6j_pair) < 1e-15) return 0.0;

    // 9j Varshalovich rows: (T, T', 1), (1/2, 1/2, 1), (T_3N, T_3N, 0)
    double w9j = wigner_9j(2*T, 2*T_prime, 2,
                            1, 1, 2,
                            two_T_3N, two_T_3N, 0);
    if (std::abs(w9j) < 1e-15) return 0.0;

    // Phase from pair isospin peel: (-1)^{1/2+1/2+T+1} = (-1)^{T}
    // (spectator has no peel phase — it's just tau_1=1/2 directly)
    double phase = (T % 2 == 0) ? 1.0 : -1.0;

    // Prefactor: -sqrt(3) * sqrt(2*T_3N+1) * <1/2||tau||1/2>^2
    //          = -sqrt(3) * sqrt(two_T_3N+1) * 6
    double prefactor = -std::sqrt(3.0) * std::sqrt(two_T_3N + 1.0) * 6.0;

    // hat(T') hat(T)
    double hat_T = std::sqrt((2*T_prime + 1.0) * (2*T + 1.0));

    return prefactor * phase * hat_T * w9j * w6j_pair;
}

double reduced_me_sigma_q_tensor(int L_prime, int S_prime, int J_2N,
                                  int L, int S, int rank) {
    // Placeholder — full implementation in Task 6
    return 0.0;
}
