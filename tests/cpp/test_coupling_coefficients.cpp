#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "coupling_coefficients.h"

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

void test_cg_known_values() {
    // All arguments are twice-integers (GSL convention)

    // CG(1/2,1/2; 1/2,1/2 | 1,1) = 1
    check_close("CG(1/2,1/2;1/2,1/2|1,1)", clebsch_gordan(1, 1, 2, 1, 1, 2), 1.0);

    // CG(1/2,1/2; 1/2,-1/2 | 1,0) = 1/sqrt(2)
    check_close("CG(1/2,1/2;1/2,-1/2|1,0)", clebsch_gordan(1, 1, 2, 1, -1, 0), 1.0/std::sqrt(2.0));

    // CG(1/2,1/2; 1/2,-1/2 | 0,0) = 1/sqrt(2)
    check_close("CG(1/2,1/2;1/2,-1/2|0,0)", clebsch_gordan(1, 1, 0, 1, -1, 0), 1.0/std::sqrt(2.0));

    // CG(1,1; 1,0 | 2,1) = 1/sqrt(2)
    check_close("CG(1,1;1,0|2,1)", clebsch_gordan(2, 2, 4, 2, 0, 2), 1.0/std::sqrt(2.0));

    // CG(3/2,3/2; 1/2,1/2 | 2,2) = 1
    check_close("CG(3/2,3/2;1/2,1/2|2,2)", clebsch_gordan(3, 1, 4, 3, 1, 4), 1.0);
}

void test_6j_known_values() {
    // {1 1 1; 1 1 1} = 1/6
    check_close("6j{1,1,1;1,1,1}", wigner_6j(2, 2, 2, 2, 2, 2), 1.0/6.0);

    // {1/2 1/2 1; 1/2 1/2 0} = 1/2
    check_close("6j{1/2,1/2,1;1/2,1/2,0}", wigner_6j(1, 1, 2, 1, 1, 0), 0.5);

    // {1/2 1/2 1; 1/2 1/2 1} = 1/6  (sympy oracle: +0.16667)
    check_close("6j{1/2,1/2,1;1/2,1/2,1}", wigner_6j(1, 1, 2, 1, 1, 2), 1.0/6.0);
}

void test_9j_known_values() {
    // {1 1 0; 1 1 0; 0 0 0} = 1/3
    check_close("9j{1,1,0;1,1,0;0,0,0}", wigner_9j(2,2,0, 2,2,0, 0,0,0), 1.0/3.0);
}

void test_triangle_violation() {
    // Triangle violation: |j1-j2| > j3 must return 0
    check_close("CG_triangle_violate", clebsch_gordan(6, 0, 2, 0, 0, 0), 0.0);
    check_close("6j_triangle_violate", wigner_6j(10, 0, 2, 0, 0, 0), 0.0);
}

void test_cg_symmetry() {
    // CG(j1,m1;j2,m2|j,m) = (-1)^(j1+j2-j) CG(j2,m2;j1,m1|j,m)
    // C++ args: cg1 = CG(j1=1,j2=2;m1=1,m2=0|j3=2,m3=1) -> twice-int (2,4,4,2,0,2)
    //           cg2 = CG(j1=2,j2=1;m1=0,m2=1|j3=2,m3=1) -> twice-int (4,2,4,0,2,2)
    // phase = (-1)^(j1+j2-j3) = (-1)^(1+2-2) = -1
    double cg1 = clebsch_gordan(2, 4, 4, 2, 0, 2);
    double cg2 = clebsch_gordan(4, 2, 4, 0, 2, 2);
    int phase_exp = (1 + 2 - 2);  // j1+j2-j = 1
    double phase = (phase_exp % 2 == 0) ? 1.0 : -1.0;
    check_close("CG_exchange_symmetry", cg1, phase * cg2);
}

void test_6j_orthogonality() {
    // sum_J (2J+1) {a b J; c d e} {a b J; c d e} = delta_{ee} / (2e+1) = 1/(2e+1)
    // Using a=b=c=d=e=1 (twice-int=2), summing J in {0,1,2} (twice-int {0,2,4})
    int a = 2, b = 2, c = 2, d = 2, e = 2;  // all j=1
    double sum = 0.0;
    for (int twoJ = 0; twoJ <= 4; twoJ += 2) {
        double w = wigner_6j(a, b, twoJ, c, d, e);
        sum += (twoJ + 1.0) * w * w;
    }
    // 1/(2e+1) where e=1 (actual), so = 1/3; in twice-int e=2, so 1/(e+1)=1/3
    check_close("6j_orthogonality", sum, 1.0 / (e + 1.0));
}

int main() {
    std::printf("=== L1: Coupling Coefficient Tests ===\n\n");

    test_cg_known_values();
    test_6j_known_values();
    test_9j_known_values();
    test_triangle_violation();
    test_cg_symmetry();
    test_6j_orthogonality();

    std::printf("\n%d passed, %d failed\n", g_passes, g_failures);
    return g_failures > 0 ? 1 : 0;
}
