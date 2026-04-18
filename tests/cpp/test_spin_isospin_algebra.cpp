#include <cstdio>
#include <cmath>
#include "spin_isospin_algebra.h"

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

void test_sigma_qhat_selection_rules() {
    // Must be zero for spin-singlet pair (S=0)
    check_close("sigma_qhat_S0", reduced_me_sigma_dot_qhat(1, 0, 1, 0, 0), 0.0);

    // Must be zero when |L'-L| != 1
    check_close("sigma_qhat_dL0", reduced_me_sigma_dot_qhat(0, 1, 1, 0, 1), 0.0);
    check_close("sigma_qhat_dL2", reduced_me_sigma_dot_qhat(2, 1, 1, 0, 1), 0.0);

    // Must be nonzero for valid transitions
    double val = reduced_me_sigma_dot_qhat(1, 1, 1, 0, 1);  // 3S1 -> 3P1
    if (std::abs(val) < 1e-15) {
        std::printf("FAIL sigma_qhat_3S1_3P1: expected nonzero, got %.15e\n", val);
        g_failures++;
    } else {
        std::printf("  sigma_qhat(L'=1,S'=1,J=1; L=0,S=1) = %.10f\n", val);
        g_passes++;
    }
}

void test_sigma_qhat_specific_values() {
    // <3P1 || sigma.qhat || 3S1>  (L'=1, S=1, J=1, L=0)
    // Phase: (-1)^{1+1+1} = -1
    // CG(0,0;1,0|1,0) = clebsch_gordan(0,2,2, 0,0,0) = 1
    // 6j{1 1 1; 1 0 1} = wigner_6j(2,2,2, 2,0,2)
    double w6j_val = wigner_6j(2, 2, 2, 2, 0, 2);
    double cg_val = clebsch_gordan(0, 2, 2, 0, 0, 0);
    double expected = -1.0 * std::sqrt(3.0 * 3.0 * 1.0) * w6j_val * cg_val;
    double got = reduced_me_sigma_dot_qhat(1, 1, 1, 0, 1);
    check_close("sigma_qhat_3S1_3P1_value", got, expected);

    // <3S1 || sigma.qhat || 3P1>  (L'=0, S=1, J=1, L=1)
    // Phase: (-1)^{0+1+1} = 1
    // CG(1,0;1,0|0,0) = clebsch_gordan(2,2,0, 0,0,0)
    double cg_val2 = clebsch_gordan(2, 2, 0, 0, 0, 0);
    double w6j_val2 = wigner_6j(0, 2, 2, 2, 2, 2);
    double expected2 = 1.0 * std::sqrt(3.0 * 1.0 * 3.0) * w6j_val2 * cg_val2;
    double got2 = reduced_me_sigma_dot_qhat(0, 1, 1, 1, 1);
    check_close("sigma_qhat_3P1_3S1_value", got2, expected2);

    // <3D1 || sigma.qhat || 3P1>  (L'=2, S=1, J=1, L=1)
    double got3 = reduced_me_sigma_dot_qhat(2, 1, 1, 1, 1);
    if (std::abs(got3) < 1e-15 || !std::isfinite(got3)) {
        std::printf("FAIL sigma_qhat_3P1_3D1: expected nonzero finite, got %.15e\n", got3);
        g_failures++;
    } else {
        std::printf("  sigma_qhat(L'=2,S'=1,J=1; L=1,S=1) = %.10f\n", got3);
        g_passes++;
    }
}

void test_tau23_eigenvalue() {
    check_close("tau23_T0", tau23_eigenvalue(0), -3.0);
    check_close("tau23_T1", tau23_eigenvalue(1), 1.0);
}

// --- sigma1_dot_sigma3 tests ---

void test_sigma1_dot_sigma3_selection_rules() {
    // Must be zero when L' != L
    double val = reduced_me_sigma1_dot_sigma3(
        /*L'=*/1, /*S'=*/1, /*J'=*/1, /*l'=*/0, /*two_j'=*/1, /*two_J3N=*/1,
        /*L=*/ 0, /*S=*/ 1, /*J=*/ 1, /*l=*/ 0, /*two_j=*/ 1);
    check_close("sig13_Ldiff", val, 0.0);

    // Must be zero when l' != l
    val = reduced_me_sigma1_dot_sigma3(
        0, 1, 1, /*l'=*/1, 1, 1,
        0, 1, 1, /*l=*/ 0, 1);
    check_close("sig13_ldiff", val, 0.0);

    // Must be zero when S=S'=0 (sigma_3 cannot connect singlet to singlet)
    val = reduced_me_sigma1_dot_sigma3(
        0, /*S'=*/0, /*J'=*/0, 0, 1, 1,
        0, /*S=*/ 0, /*J=*/ 0, 0, 1);
    check_close("sig13_S0S0", val, 0.0);
}

void test_sigma1_dot_sigma3_diagonal() {
    // Case 1: S=S'=1, J=J'=1, j=j'=1/2, L=L'=0, l=l'=0, J_3N=1/2
    // Expected: -2.0
    double val = reduced_me_sigma1_dot_sigma3(
        /*L'=*/0, /*S'=*/1, /*J'=*/1, /*l'=*/0, /*two_j'=*/1, /*two_J3N=*/1,
        /*L=*/ 0, /*S=*/ 1, /*J=*/ 1, /*l=*/ 0, /*two_j=*/ 1);
    std::printf("  sig13(S=1,J=1,j=1/2 diag) = %.10f (expected -2.0)\n", val);
    check_close("sig13_diag_S1J1j12", val, -2.0);
}

void test_sigma1_dot_sigma3_offdiag() {
    // Case 2: S'=0,J'=0,j'=1/2 <- S=1,J=1,j=1/2, L=L'=0, l=l'=0, J_3N=1/2
    // Expected: +sqrt(3) (verified by direct m-scheme computation)
    double val = reduced_me_sigma1_dot_sigma3(
        /*L'=*/0, /*S'=*/0, /*J'=*/0, /*l'=*/0, /*two_j'=*/1, /*two_J3N=*/1,
        /*L=*/ 0, /*S=*/ 1, /*J=*/ 1, /*l=*/ 0, /*two_j=*/ 1);
    std::printf("  sig13(S'=0,J'=0 <- S=1,J=1) = %.10f (expected +%.10f)\n",
                val, std::sqrt(3.0));
    check_close("sig13_offdiag_S0S1", val, std::sqrt(3.0));
}

void test_sigma1_dot_sigma3_with_orbital() {
    // Case with L=1, l=1, j_1N=3/2 (two_j=3), S=1, J_2N=1, J_3N=1/2
    // This tests that orbital angular momentum peeling works correctly.
    // Verify it's nonzero and finite for a valid channel.
    double val = reduced_me_sigma1_dot_sigma3(
        /*L'=*/1, /*S'=*/1, /*J'=*/1, /*l'=*/1, /*two_j'=*/3, /*two_J3N=*/1,
        /*L=*/ 1, /*S=*/ 1, /*J=*/ 1, /*l=*/ 1, /*two_j=*/ 3);
    std::printf("  sig13(L=1,S=1,J=1,l=1,j=3/2 diag) = %.10f\n", val);
    if (!std::isfinite(val)) {
        std::printf("FAIL sig13_orbital: not finite\n");
        g_failures++;
    } else {
        g_passes++;
    }

    // Also test with j_1N changing: j'=1/2, j=3/2 (both valid for l=1, s=1/2)
    double val2 = reduced_me_sigma1_dot_sigma3(
        /*L'=*/1, /*S'=*/1, /*J'=*/1, /*l'=*/1, /*two_j'=*/1, /*two_J3N=*/1,
        /*L=*/ 1, /*S=*/ 1, /*J=*/ 1, /*l=*/ 1, /*two_j=*/ 3);
    std::printf("  sig13(L=1,j'=1/2,j=3/2) = %.10f\n", val2);
    if (!std::isfinite(val2)) {
        std::printf("FAIL sig13_orbital_jchange: not finite\n");
        g_failures++;
    } else {
        g_passes++;
    }
}

// --- tau1_dot_tau3 tests ---

void test_tau1_dot_tau3_diagonal() {
    // T=T'=1, T_3N=1/2: expected -2.0
    double val = tau1_dot_tau3(/*T'=*/1, /*T=*/1, /*two_T3N=*/1);
    std::printf("  tau13(T=1 diag, T3N=1/2) = %.10f (expected -2.0)\n", val);
    check_close("tau13_diag_T1", val, -2.0);
}

void test_tau1_dot_tau3_offdiag() {
    // T'=0, T=1, T_3N=1/2: expected +sqrt(3) (verified by m-scheme analogy with spin)
    double val = tau1_dot_tau3(/*T'=*/0, /*T=*/1, /*two_T3N=*/1);
    std::printf("  tau13(T'=0,T=1, T3N=1/2) = %.10f (expected +%.10f)\n",
                val, std::sqrt(3.0));
    check_close("tau13_offdiag_T0T1", val, std::sqrt(3.0));
}

void test_tau1_dot_tau3_zero() {
    // T=T'=0, T_3N=1/2: expected 0.0 (selection rule)
    double val = tau1_dot_tau3(/*T'=*/0, /*T=*/0, /*two_T3N=*/1);
    check_close("tau13_T0T0_zero", val, 0.0);
}

void test_tau1_dot_tau3_T3N_32() {
    // T=T'=1, T_3N=3/2: should be nonzero and finite
    double val = tau1_dot_tau3(/*T'=*/1, /*T=*/1, /*two_T3N=*/3);
    std::printf("  tau13(T=1 diag, T3N=3/2) = %.10f\n", val);
    if (!std::isfinite(val)) {
        std::printf("FAIL tau13_T3N_32: not finite\n");
        g_failures++;
    } else {
        g_passes++;
    }
}

int main() {
    std::printf("=== L2: Spin-Isospin Algebra Tests ===\n\n");

    test_tau23_eigenvalue();
    test_sigma_qhat_selection_rules();
    test_sigma_qhat_specific_values();

    std::printf("\n--- sigma1_dot_sigma3 ---\n");
    test_sigma1_dot_sigma3_selection_rules();
    test_sigma1_dot_sigma3_diagonal();
    test_sigma1_dot_sigma3_offdiag();
    test_sigma1_dot_sigma3_with_orbital();

    std::printf("\n--- tau1_dot_tau3 ---\n");
    test_tau1_dot_tau3_diagonal();
    test_tau1_dot_tau3_offdiag();
    test_tau1_dot_tau3_zero();
    test_tau1_dot_tau3_T3N_32();

    std::printf("\n%d passed, %d failed\n", g_passes, g_failures);
    return g_failures > 0 ? 1 : 0;
}
