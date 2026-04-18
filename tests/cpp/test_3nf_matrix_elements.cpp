#include <cstdio>
#include <cmath>
#include "chiral_N2LO_3NF.h"
#include "make_pw_symm_states.h"

static int g_failures = 0;
static int g_passes = 0;

static void check_close(const char* label, double got, double expected, double tol = 1e-10) {
    if (std::abs(got - expected) > tol) {
        std::printf("FAIL %s: got %.15e, expected %.15e\n", label, got, expected);
        g_failures++;
    } else {
        g_passes++;
    }
}

// Build a minimal pw_3N_statespace for testing
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

void test_contact_zero_when_cE_zero() {
    chiral_N2LO_3NF tnf(0.0, 0.0, 500.0);  // c_D=0, c_E=0
    pw_3N_statespace pw = make_test_pw_states();
    for (int a = 0; a < pw.Nalpha; a++) {
        double val = tnf.W1_element(a, a, 1.0, 1.0, 1.0, 1.0, pw);
        check_close("CT_zero_cE0", val, 0.0);
    }
}

void test_contact_diagonal() {
    chiral_N2LO_3NF tnf(0.0, -0.398, 500.0);  // c_E only
    pw_3N_statespace pw = make_test_pw_states();
    for (int a = 0; a < pw.Nalpha; a++) {
        for (int b = 0; b < pw.Nalpha; b++) {
            if (a == b) continue;
            if (pw.two_J_3N_array[a] != pw.two_J_3N_array[b]) continue;
            if (pw.two_T_3N_array[a] != pw.two_T_3N_array[b]) continue;
            if (pw.P_3N_array[a] != pw.P_3N_array[b]) continue;
            double val = tnf.W1_element(a, b, 1.0, 1.0, 1.0, 1.0, pw);
            check_close("CT_offdiag_zero", val, 0.0);
        }
    }
}

void test_1pe_ct_zero_when_cD_zero() {
    chiral_N2LO_3NF tnf(0.0, -0.398, 500.0);  // c_D=0
    pw_3N_statespace pw = make_test_pw_states();
    for (int a = 0; a < pw.Nalpha; a++) {
        for (int b = 0; b < pw.Nalpha; b++) {
            if (a == b) continue;
            if (pw.two_J_3N_array[a] != pw.two_J_3N_array[b]) continue;
            if (pw.two_T_3N_array[a] != pw.two_T_3N_array[b]) continue;
            if (pw.P_3N_array[a] != pw.P_3N_array[b]) continue;
            double val = tnf.W1_element(a, b, 1.0, 1.0, 1.0, 1.0, pw);
            check_close("1PE_zero_cD0", val, 0.0);
        }
    }
}

void test_1pe_ct_nonzero_when_cD_set() {
    chiral_N2LO_3NF tnf(-0.2, 0.0, 500.0);  // c_D=-0.2, c_E=0
    pw_3N_statespace pw = make_test_pw_states();
    bool found_nonzero = false;
    for (int a = 0; a < pw.Nalpha && !found_nonzero; a++) {
        for (int b = 0; b < pw.Nalpha && !found_nonzero; b++) {
            if (a == b) continue;
            if (pw.two_J_3N_array[a] != pw.two_J_3N_array[b]) continue;
            if (pw.two_T_3N_array[a] != pw.two_T_3N_array[b]) continue;
            if (pw.P_3N_array[a] != pw.P_3N_array[b]) continue;
            double val = tnf.W1_element(a, b, 1.0, 1.0, 1.0, 1.0, pw);
            if (std::abs(val) > 1e-15) found_nonzero = true;
        }
    }
    if (!found_nonzero) {
        std::printf("FAIL 1PE_nonzero: no off-diagonal elements found with c_D=-0.2\n");
        g_failures++;
    } else {
        g_passes++;
    }
}

int main() {
    std::printf("=== L2: 3NF Matrix Element Tests ===\n\n");

    test_contact_zero_when_cE_zero();
    test_contact_diagonal();
    test_1pe_ct_zero_when_cD_zero();
    test_1pe_ct_nonzero_when_cD_set();

    std::printf("\n%d passed, %d failed\n", g_passes, g_failures);
    return g_failures > 0 ? 1 : 0;
}
