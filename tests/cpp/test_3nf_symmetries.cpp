#include <cstdio>
#include <cmath>
#include "chiral_N2LO_3NF.h"
#include "make_pw_symm_states.h"

static int g_failures = 0;
static int g_passes = 0;

static void check_close(const char* label, double got, double expected, double tol = 1e-12) {
    if (std::abs(got - expected) > tol) {
        std::printf("FAIL %s: got %.15e, expected %.15e\n", label, got, expected);
        g_failures++;
    } else {
        g_passes++;
    }
}

static pw_3N_statespace make_test_pw_states() {
    run_params rp = {};
    rp.J_2N_max = 2;
    rp.two_J_3N_max = 3;
    rp.tensor_force = true;
    rp.isospin_breaking_1S0 = false;
    pw_3N_statespace pw = {};
    construct_symmetric_pw_states(pw, rp);
    return pw;
}

void test_j3n_conservation() {
    chiral_N2LO_3NF tnf(-0.2, -0.398, 500.0, -0.81, -3.2, 0.0);
    pw_3N_statespace pw = make_test_pw_states();
    int count = 0;
    for (int a = 0; a < pw.Nalpha; a++) {
        for (int b = 0; b < pw.Nalpha; b++) {
            if (pw.two_J_3N_array[a] == pw.two_J_3N_array[b]) continue;
            double val = tnf.W1_element(a, b, 1.0, 1.0, 1.0, 1.0, pw);
            check_close("J3N_conserve", val, 0.0);
            count++;
        }
    }
    std::printf("  Checked %d J_3N-violating pairs\n", count);
}

void test_t3n_conservation() {
    chiral_N2LO_3NF tnf(-0.2, -0.398, 500.0, -0.81, -3.2, 0.0);
    pw_3N_statespace pw = make_test_pw_states();
    int count = 0;
    for (int a = 0; a < pw.Nalpha; a++) {
        for (int b = 0; b < pw.Nalpha; b++) {
            if (pw.two_T_3N_array[a] == pw.two_T_3N_array[b]) continue;
            double val = tnf.W1_element(a, b, 1.0, 1.0, 1.0, 1.0, pw);
            check_close("T3N_conserve", val, 0.0);
            count++;
        }
    }
    std::printf("  Checked %d T_3N-violating pairs\n", count);
}

void test_parity_conservation() {
    chiral_N2LO_3NF tnf(-0.2, -0.398, 500.0, -0.81, -3.2, 0.0);
    pw_3N_statespace pw = make_test_pw_states();
    int count = 0;
    for (int a = 0; a < pw.Nalpha; a++) {
        for (int b = 0; b < pw.Nalpha; b++) {
            if (pw.P_3N_array[a] == pw.P_3N_array[b]) continue;
            double val = tnf.W1_element(a, b, 1.0, 1.0, 1.0, 1.0, pw);
            check_close("P3N_conserve", val, 0.0);
            count++;
        }
    }
    std::printf("  Checked %d parity-violating pairs\n", count);
}

void test_hermiticity() {
    chiral_N2LO_3NF tnf(-0.2, -0.398, 500.0, -0.81, -3.2, 0.0);
    pw_3N_statespace pw = make_test_pw_states();
    double p1 = 1.0, q1 = 0.8, p2 = 1.2, q2 = 0.6;
    int count = 0;
    for (int a = 0; a < pw.Nalpha; a++) {
        for (int b = a; b < pw.Nalpha; b++) {
            if (pw.two_J_3N_array[a] != pw.two_J_3N_array[b]) continue;
            if (pw.two_T_3N_array[a] != pw.two_T_3N_array[b]) continue;
            if (pw.P_3N_array[a] != pw.P_3N_array[b]) continue;

            double w_ab = tnf.W1_element(a, b, p1, q1, p2, q2, pw);
            double w_ba = tnf.W1_element(b, a, p2, q2, p1, q1, pw);
            check_close("Hermitian", w_ab, w_ba, 1e-12);
            count++;
        }
    }
    std::printf("  Checked %d Hermiticity pairs\n", count);
}

int main() {
    std::printf("=== L3: 3NF Symmetry Tests ===\n\n");

    test_j3n_conservation();
    test_t3n_conservation();
    test_parity_conservation();
    test_hermiticity();

    std::printf("\n%d passed, %d failed\n", g_passes, g_failures);
    return g_failures > 0 ? 1 : 0;
}
