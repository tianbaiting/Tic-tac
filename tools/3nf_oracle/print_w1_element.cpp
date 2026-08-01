// tools/3nf_oracle/print_w1_element.cpp
//
// Phase 3 oracle driver. Constructs a chiral_N2LO_c1c3cDcE_approx 3NF model
// and prints W1_element(alpha_r, alpha_c, p_r, q_r, p_c, q_c) for the channels
// and momenta requested on stdin. The Python angular oracle
// (angular_oracle.py) reads these values and compares them to an INDEPENDENT
// full multi-dimensional angular quadrature of the un-reduced momentum-space
// operator — NOT a re-transcription of the C++ kernel formula.
//
// Usage:
//   print_w1_element cE cD c1 c3 Lambda_MeV
//   then feed lines:  alpha_r alpha_c p_r q_r p_c q_c
//   (momenta in fm^-1). End with EOF (Ctrl-D).
//
// The driver writes "W1 alpha_r alpha_c p_r q_r p_c q_c value" per line.

#include <cstdio>
#include <iostream>
#include "chiral_N2LO_3NF.h"
#include "make_pw_symm_states.h"

static pw_3N_statespace make_test_pw_states() {
    run_params rp = {};
    rp.J_2N_max = 1;
    rp.two_J_3N_max = 1;
    rp.tensor_force = true;
    rp.isospin_breaking_1S0 = false;
    pw_3N_statespace pw = {};
    construct_symmetric_pw_states(pw, rp);
    return pw;
}

static int find_alpha(const pw_3N_statespace& pw,
                      int L2, int S2, int J2, int T2,
                      int l1, int two_j1, int two_J3, int two_T3) {
    for (int a = 0; a < pw.Nalpha; ++a) {
        if (pw.L_2N_array[a]==L2 && pw.S_2N_array[a]==S2 && pw.J_2N_array[a]==J2 &&
            pw.T_2N_array[a]==T2 && pw.L_1N_array[a]==l1 && pw.two_J_1N_array[a]==two_j1 &&
            pw.two_J_3N_array[a]==two_J3 && pw.two_T_3N_array[a]==two_T3) return a;
    }
    return -1;
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::fprintf(stderr, "usage: %s cE cD c1 c3 Lambda_MeV\n", argv[0]);
        return 2;
    }
    double cE = std::stod(argv[1]);
    double cD = std::stod(argv[2]);
    double c1 = std::stod(argv[3]);
    double c3 = std::stod(argv[4]);
    double Lambda = std::stod(argv[5]);

    chiral_N2LO_3NF tnf(cD, cE, Lambda, c1, c3, /*c4=*/0.0);
    pw_3N_statespace pw = make_test_pw_states();

    // emit the channel table (one line per alpha) so the Python oracle can
    // resolve quantum numbers without re-deriving the state space.
    std::printf("# channels Nalpha=%d\n", pw.Nalpha);
    for (int a = 0; a < pw.Nalpha; ++a) {
        std::printf("# alpha %d  L2=%d S2=%d J2=%d T2=%d l1=%d 2j1=%d 2J3=%d 2T3=%d P3=%d\n",
            a, pw.L_2N_array[a], pw.S_2N_array[a], pw.J_2N_array[a], pw.T_2N_array[a],
            pw.L_1N_array[a], pw.two_J_1N_array[a], pw.two_J_3N_array[a],
            pw.two_T_3N_array[a], pw.P_3N_array[a]);
    }

    double ar, ac, pr, qr, pc, qc;
    while (std::cin >> ar >> ac >> pr >> qr >> pc >> qc) {
        int ar_i = (int)ar, ac_i = (int)ac;
        if (ar_i < 0 || ac_i < 0) continue;
        double v = tnf.W1_element(ar_i, ac_i, pr, qr, pc, qc, pw);
        std::printf("W1 %d %d %.10f %.10f %.10f %.10f %.12e\n",
                    ar_i, ac_i, pr, qr, pc, qc, v);
    }
    return 0;
}
