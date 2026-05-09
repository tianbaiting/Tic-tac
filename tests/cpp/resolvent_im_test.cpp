// tests/cpp/resolvent_im_test.cpp
//
// [EN] Unit-test the analytic Im parts of the cell-averaged resolvent.
// Reference values are derived from the Heaviside-step structure documented in
// the WPCD literature (Kukulin 2007, Eqs. 32, 34-35) and replicated in
// docs/im_path_diagnosis_2026-05-09.md.
// [CN] 单元测试 cell-averaged resolvent 的解析虚部，期望值来自
// Heaviside 阶跃结构。

#include "core/resolvent/make_resolvent.h"
#include <cassert>
#include <cmath>
#include <cstdio>

static bool close(double a, double b, double tol=1e-10) {
    return std::abs(a - b) <= tol * (1.0 + std::abs(b));
}

int main() {
    // Test 1: Bound-continuum on-shell q-bin.
    // Eb = -2.224 MeV, q-bin [50, 60] MeV/c (linear momentum, then converted).
    // E chosen so that Eq_lower + Eb < E < Eq_upper + Eb.
    {
        double Eb = -2.224;
        double q_lo = 50.0, q_hi = 60.0;
        // mu1 from constants.h matches the formula in resolvent_bound_continuum
        double mu1 = 939.565 * (939.565 + 938.272 + Eb) /
                     (939.565 + 939.565 + 938.272 + Eb);
        double Eq_lo = q_lo*q_lo / (2*mu1);
        double Eq_hi = q_hi*q_hi / (2*mu1);
        double Dq = Eq_hi - Eq_lo;
        double E = Eb + 0.5*(Eq_lo + Eq_hi);   // mid-bin on-shell

        cdouble R = resolvent_bound_continuum(E, Eb, q_hi, q_lo);

        // Expected Im_R = -π / Δq when E falls strictly inside the bin.
        double expected_im = -M_PI / Dq;
        if (!close(R.imag(), expected_im, 1e-9)) {
            std::printf("FAIL Test1: Im_R = %.10e, expected %.10e\n",
                        R.imag(), expected_im);
            return 1;
        }
        std::printf("PASS Test1: Im_R = %.10e ~= -pi/Dq\n", R.imag());
    }

    // Test 2: Continuum-continuum cell that strictly straddles E.
    // Im_Q should be non-zero with the documented sign.
    {
        double Eb = -2.224;
        double q_lo = 50.0, q_hi = 60.0;
        double e_lo = 5.0,  e_hi = 15.0;     // p-bin in MeV
        double mu1 = 939.565 * (939.565 + 938.272 + Eb) /
                     (939.565 + 939.565 + 938.272 + Eb);
        double Eq_lo = q_lo*q_lo / (2*mu1);
        double Eq_hi = q_hi*q_hi / (2*mu1);
        double E = e_lo + Eq_lo + 0.5 * ((e_hi - e_lo) + (Eq_hi - Eq_lo));

        cdouble Q = resolvent_continuum_continuum(E, Eb, q_hi, q_lo, e_hi, e_lo);

        if (std::abs(Q.imag()) < 1e-9) {
            std::printf("FAIL Test2: Im_Q = %.10e (expected non-zero straddle)\n", Q.imag());
            return 1;
        }
        std::printf("PASS Test2: Im_Q = %.10e (non-zero)\n", Q.imag());
    }

    // Test 3: Cell well below E -- Im_Q must be zero.
    {
        double Eb = -2.224;
        double E = 1000.0;  // much higher than any cell energy
        cdouble Q = resolvent_continuum_continuum(E, Eb, 60, 50, 15, 5);
        if (std::abs(Q.imag()) > 1e-9) {
            std::printf("FAIL Test3: Im_Q = %.10e (expected ~0 for cell << E)\n", Q.imag());
            return 1;
        }
        std::printf("PASS Test3: Im_Q ~= 0 below threshold\n");
    }

    std::printf("\nAll resolvent_im_test cases passed.\n");
    return 0;
}
