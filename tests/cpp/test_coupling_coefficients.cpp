#include <cstdio>
#include <cmath>
#include <cstdlib>
#include "coupling_coefficients.h"

static int g_failures = 0;

static void check_close(const char* label, double got, double expected, double tol = 1e-12) {
    if (std::abs(got - expected) > tol) {
        std::printf("FAIL %s: got %.15e, expected %.15e, diff %.2e\n", label, got, expected, std::abs(got - expected));
        g_failures++;
    } else {
        std::printf("PASS %s\n", label);
    }
}

int main() {
    // Placeholder — real tests added in Task 2
    check_close("placeholder", 1.0, 1.0);

    std::printf("\n%d failure(s)\n", g_failures);
    return g_failures > 0 ? 1 : 0;
}
