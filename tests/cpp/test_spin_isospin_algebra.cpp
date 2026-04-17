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

int main() {
    std::printf("=== L2: Spin-Isospin Algebra Tests ===\n\n");

    test_tau23_eigenvalue();
    test_sigma_qhat_selection_rules();
    test_sigma_qhat_specific_values();

    std::printf("\n%d passed, %d failed\n", g_passes, g_failures);
    return g_failures > 0 ? 1 : 0;
}
