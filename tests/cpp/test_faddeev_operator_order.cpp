// Dense, noncommuting verification of the 3NF ordering in the symmetrized AGS
// equation used by the WPCD solver.
//
// Deltuva, Phys. Rev. C 80, 064002 (2009), Eq. (7a), gives
//
//   U = P G0^-1 + (1+P)u + P t G0 U
//       + (1+P)u G0(1+tG0)U .
//
// With G1=G0+G0 t G0, tG0=vG1, and G0(1+tG0)=G1, its
// iteration kernel is exactly
//
//   [P v + (1+P)u] G1 .
//
// The reversed u(1+P) ordering belongs to a different Faddeev-component
// equation and is not interchangeable when [P,u] != 0.  This test evaluates
// both sides with small Hermitian v and u matrices and a physical identical-
// particle permutation sum satisfying P^2=P+2.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>

namespace {

constexpr int N = 3;
using Matrix = std::array<double, N * N>;

int failures = 0;
int passes = 0;

Matrix identity() {
    Matrix a{};
    for (int i = 0; i < N; ++i) a[i * N + i] = 1.0;
    return a;
}

Matrix add(const Matrix& a, const Matrix& b) {
    Matrix c{};
    for (int i = 0; i < N * N; ++i) c[i] = a[i] + b[i];
    return c;
}

Matrix subtract(const Matrix& a, const Matrix& b) {
    Matrix c{};
    for (int i = 0; i < N * N; ++i) c[i] = a[i] - b[i];
    return c;
}

Matrix multiply(const Matrix& a, const Matrix& b) {
    Matrix c{};
    for (int i = 0; i < N; ++i)
        for (int k = 0; k < N; ++k)
            for (int j = 0; j < N; ++j)
                c[i * N + j] += a[i * N + k] * b[k * N + j];
    return c;
}

Matrix inverse(Matrix a) {
    Matrix inv = identity();
    for (int col = 0; col < N; ++col) {
        int pivot = col;
        for (int row = col + 1; row < N; ++row)
            if (std::abs(a[row * N + col]) > std::abs(a[pivot * N + col])) pivot = row;
        if (pivot != col) {
            for (int j = 0; j < N; ++j) {
                std::swap(a[col * N + j], a[pivot * N + j]);
                std::swap(inv[col * N + j], inv[pivot * N + j]);
            }
        }
        const double diagonal = a[col * N + col];
        if (std::abs(diagonal) < 1e-14) {
            std::printf("FAIL inverse: singular toy matrix\n");
            failures++;
            return {};
        }
        for (int j = 0; j < N; ++j) {
            a[col * N + j] /= diagonal;
            inv[col * N + j] /= diagonal;
        }
        for (int row = 0; row < N; ++row) {
            if (row == col) continue;
            const double factor = a[row * N + col];
            for (int j = 0; j < N; ++j) {
                a[row * N + j] -= factor * a[col * N + j];
                inv[row * N + j] -= factor * inv[col * N + j];
            }
        }
    }
    return inv;
}

double max_abs(const Matrix& a) {
    double result = 0.0;
    for (double value : a) result = std::max(result, std::abs(value));
    return result;
}

double max_diff(const Matrix& a, const Matrix& b) {
    return max_abs(subtract(a, b));
}

void require_close(const char* label, const Matrix& got, const Matrix& expected,
                   double tolerance = 2e-13) {
    const double difference = max_diff(got, expected);
    if (difference > tolerance) {
        std::printf("FAIL %s: max|delta|=%.6e\n", label, difference);
        failures++;
    } else {
        std::printf("  PASS %s: max|delta|=%.6e\n", label, difference);
        passes++;
    }
}

void require_distinct(const char* label, const Matrix& a, const Matrix& b,
                      double threshold = 1e-6) {
    const double difference = max_diff(a, b);
    if (difference <= threshold) {
        std::printf("FAIL %s: max|delta|=%.6e\n", label, difference);
        failures++;
    } else {
        std::printf("  PASS %s: max|delta|=%.6e\n", label, difference);
        passes++;
    }
}

}  // namespace

int main() {
    const Matrix I = identity();

    // P=C_123+C_132 in the three-dimensional regular representation.
    const Matrix P = {
        0.0, 1.0, 1.0,
        1.0, 0.0, 1.0,
        1.0, 1.0, 0.0,
    };
    const Matrix two_I = {
        2.0, 0.0, 0.0,
        0.0, 2.0, 0.0,
        0.0, 0.0, 2.0,
    };
    require_close("identical-particle P^2=P+2I",
                  multiply(P, P), add(P, two_I));

    // Hermitian pair potential and bare spectator component.  W is chosen so
    // that it does not commute with P; no unphysical matrix asymmetry is needed
    // to expose the ordering.
    const Matrix V = {
         0.20,  0.03,  0.00,
         0.03, -0.15,  0.04,
         0.00,  0.04,  0.10,
    };
    const Matrix W = {
         0.11,  0.02, -0.01,
         0.02, -0.07,  0.05,
        -0.01,  0.05,  0.16,
    };
    require_distinct("[P,W] is nonzero", multiply(P, W), multiply(W, P));

    // G0 commutes with P, as the free resolvent must.  Build t from the exact
    // finite-dimensional Lippmann-Schwinger equation t=v+vG0t.
    Matrix G0{};
    for (int i = 0; i < N; ++i) G0[i * N + i] = -0.4;
    const Matrix t = multiply(inverse(subtract(I, multiply(V, G0))), V);
    const Matrix G1 = add(G0, multiply(multiply(G0, t), G0));

    require_close("tG0=vG1 identity", multiply(t, G0), multiply(V, G1));
    require_close("G0(1+tG0)=G1 identity",
                  multiply(G0, add(I, multiply(t, G0))), G1);

    const Matrix one_plus_P = add(I, P);
    const Matrix deltuva_kernel = add(
        multiply(multiply(P, t), G0),
        multiply(multiply(multiply(one_plus_P, W), G0),
                 add(I, multiply(t, G0))));

    const Matrix reduced_kernel = multiply(
        add(multiply(P, V), multiply(one_plus_P, W)), G1);
    require_close("Deltuva Eq.(7a) reduces to [Pv+(1+P)W]G1",
                  reduced_kernel, deltuva_kernel);

    const Matrix reversed_kernel = multiply(
        add(multiply(P, V), multiply(W, one_plus_P)), G1);
    require_distinct("reversed [Pv+W(1+P)]G1 is rejected",
                     reversed_kernel, deltuva_kernel);

    std::printf("\n%d passed, %d failed\n", passes, failures);
    return failures == 0 ? 0 : 1;
}
