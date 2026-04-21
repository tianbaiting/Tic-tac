#include <cstdio>
#include <cmath>
#include "chiral_3nf_recoupling.h"
#include "coupling_coefficients.h"

static int g_failures = 0;
static int g_passes = 0;

static void check_close(const char* label, double got, double expected, double tol = 1e-10) {
    if (std::abs(got - expected) > tol) {
        std::printf("FAIL %s: got %.15e, expected %.15e (diff %.3e)\n",
                    label, got, expected, std::abs(got - expected));
        g_failures++;
    } else {
        std::printf("  PASS %s = %.10f\n", label, got);
        g_passes++;
    }
}

static void check_nonzero(const char* label, double got) {
    if (std::abs(got) < 1e-15 || !std::isfinite(got)) {
        std::printf("FAIL %s: expected finite nonzero, got %.15e\n", label, got);
        g_failures++;
    } else {
        std::printf("  PASS %s = %.10f (nonzero)\n", label, got);
        g_passes++;
    }
}

// ---------------------------------------------------------------------------
// Scalar recoupling: 3S1 diagonal channel
// ---------------------------------------------------------------------------
// alpha = { L_2N=0, S_2N=1, J_2N=1, T_2N=0, l_1N=0, 2j_1N=1, 2J_3N=1, 2T_3N=1 }
//   sigma_2.sigma_3 (S=1)  = +1
//   tau_2.tau_3     (T=0)  = -3
//   scalar recoupling      = (+1) * (-3) = -3
void test_scalar_3S1_diagonal() {
    double val = recoupling_3nf_scalar(
        /*L_2N_r=*/0, /*S_2N_r=*/1, /*J_2N_r=*/1, /*T_2N_r=*/0,
        /*L_1N_r=*/0, /*two_J_1N_r=*/1, /*two_J_3N=*/1,
        /*L_2N_c=*/0, /*S_2N_c=*/1, /*J_2N_c=*/1, /*T_2N_c=*/0,
        /*L_1N_c=*/0, /*two_J_1N_c=*/1,
        /*two_T_3N=*/1);
    check_close("scalar 3S1 diagonal (sigma.sigma * tau.tau)", val, -3.0);
}

void test_scalar_1S0_diagonal() {
    // 1S0: L_2N=0, S_2N=0, J_2N=0, T_2N=1, l_1N=0, 2j_1N=1, 2J_3N=1, 2T_3N=1
    //   sigma_2.sigma_3 (S=0) = -3
    //   tau_2.tau_3     (T=1) = +1
    //   product               = -3
    double val = recoupling_3nf_scalar(
        0, 0, 0, 1, 0, 1, 1,
        0, 0, 0, 1, 0, 1,
        1);
    check_close("scalar 1S0 diagonal", val, -3.0);
}

// ---------------------------------------------------------------------------
// Scalar recoupling: 3S1 <-> 3D1 must be zero (ΔL=2 selection rule violation
// for a pair-scalar operator).
// ---------------------------------------------------------------------------
void test_scalar_3S1_3D1_zero() {
    double val = recoupling_3nf_scalar(
        /*L_2N_r=*/0, /*S_2N_r=*/1, /*J_2N_r=*/1, /*T_2N_r=*/0,
        /*L_1N_r=*/0, /*two_J_1N_r=*/1, /*two_J_3N=*/1,
        /*L_2N_c=*/2, /*S_2N_c=*/1, /*J_2N_c=*/1, /*T_2N_c=*/0,
        /*L_1N_c=*/0, /*two_J_1N_c=*/1,
        /*two_T_3N=*/1);
    check_close("scalar 3S1<->3D1 (should be zero)", val, 0.0);
}

void test_scalar_spectator_mismatch_zero() {
    // Mismatched spectator l_1N should zero the scalar recoupling.
    double val = recoupling_3nf_scalar(
        0, 1, 1, 0, /*l=*/0, 1, 1,
        0, 1, 1, 0, /*l=*/2, 1,
        1);
    check_close("scalar spectator l mismatch", val, 0.0);
}

// ---------------------------------------------------------------------------
// Rank-2 recoupling: 3S1 <-> 3D1 must be nonzero.
// This is the piece that opens L=0 <-> L'=2 through the spatial Y_2 and
// the spin rank-2 [sigma_2 x sigma_3]_2.
// ---------------------------------------------------------------------------
void test_rank2_3S1_3D1_nonzero() {
    double val = recoupling_3nf_rank2(
        /*L_2N_r=*/0, /*S_2N_r=*/1, /*J_2N_r=*/1, /*T_2N_r=*/0,
        /*L_1N_r=*/0, /*two_J_1N_r=*/1, /*two_J_3N=*/1,
        /*L_2N_c=*/2, /*S_2N_c=*/1, /*J_2N_c=*/1, /*T_2N_c=*/0,
        /*L_1N_c=*/0, /*two_J_1N_c=*/1,
        /*two_T_3N=*/1);
    check_nonzero("rank2 3S1<->3D1 (primary c_1/c_3 sign-flip channel)", val);
}

// ---------------------------------------------------------------------------
// Rank-2 recoupling: selection rules.
// ---------------------------------------------------------------------------
void test_rank2_singlet_zero() {
    // S=S'=0 (no rank-2 on singlet): must return zero.
    double val = recoupling_3nf_rank2(
        0, /*S=*/0, 0, 1, 0, 1, 1,
        0, /*S=*/0, 0, 1, 0, 1,
        1);
    check_close("rank2 S=0 singlet (should be zero)", val, 0.0);
}

void test_rank2_T_mismatch_zero() {
    // Pair isospin diagonal: T' != T must be zero.
    double val = recoupling_3nf_rank2(
        0, 1, 1, /*T=*/0, 0, 1, 1,
        2, 1, 1, /*T=*/1, 0, 1,
        1);
    check_close("rank2 T pair mismatch", val, 0.0);
}

void test_rank2_spectator_mismatch_zero() {
    // Spectator untouched -> l' != l and 2j' != 2j both zero the rank-2.
    double val = recoupling_3nf_rank2(
        0, 1, 1, 0, /*l=*/0, 1, 1,
        2, 1, 1, 0, /*l=*/2, 1,
        1);
    check_close("rank2 spectator l mismatch", val, 0.0);

    val = recoupling_3nf_rank2(
        0, 1, 1, 0, 0, /*2j=*/1, 1,
        2, 1, 1, 0, 0, /*2j=*/3,
        1);
    check_close("rank2 spectator j mismatch", val, 0.0);
}

void test_rank2_dL_1_zero() {
    // Rank-2 spatial: dL=1 must be zero (CG parity 2+L+L' even required).
    double val = recoupling_3nf_rank2(
        0, 1, 1, 0, 0, 1, 1,
        1, 1, 1, 0, 0, 1,
        1);
    check_close("rank2 dL=1 parity forbidden", val, 0.0);
}

// ---------------------------------------------------------------------------
// Rank-2 recoupling: explicit numeric value for the primary channel.
// Benchmark by direct formula evaluation:
//   tau.tau (T=0) = -3
//   sqrt(30)
//   hat = sqrt((2*1+1)*(2*1+1)*(2*1+1)*(2*1+1)) = sqrt(81) = 9
//   CG(L=0, 0; 2, 0 | L'=2, 0) = 1       (standard; CG(0,m;l,0|l,m) = 1)
//   9j{0 1 1; 2 1 1; 2 2 0} from GSL
// Final value = -3 * sqrt(30) * 9 * 1 * 9j
// ---------------------------------------------------------------------------
void test_rank2_3S1_3D1_value() {
    // Use code convention: row = bra (primed), column = ket (unprimed).
    // Here bra = 3S1 (L_r=0), ket = 3D1 (L_c=2).
    // Formula: CG(L_ket, 0; 2, 0 | L_bra, 0) = CG(L_c, 0; 2, 0 | L_r, 0)
    //                                        = CG(2, 0; 2, 0 | 0, 0) = 1/sqrt(5)
    double cg_ref = clebsch_gordan(/*2*L_c=*/4, /*2*2=*/4, /*2*L_r=*/0,
                                    0, 0, 0);
    // 9j{L_r S_r J_r; L_c S_c J_c; 2 2 0} = 9j{0 1 1; 2 1 1; 2 2 0}
    double w9j_ref = wigner_9j(0, 2, 2,
                               4, 2, 2,
                               4, 4, 0);
    double expected = -3.0 * std::sqrt(30.0) * 9.0 * cg_ref * w9j_ref;

    double got = recoupling_3nf_rank2(
        0, 1, 1, 0, 0, 1, 1,
        2, 1, 1, 0, 0, 1,
        1);

    std::printf("  [ref] CG(0 0; 2 0 | 2 0) = %.10f\n", cg_ref);
    std::printf("  [ref] 9j{0 1 1; 2 1 1; 2 2 0} = %.10f\n", w9j_ref);
    std::printf("  [ref] expected = %.10f\n", expected);
    check_close("rank2 3S1<->3D1 numeric value", got, expected);
}

// ---------------------------------------------------------------------------
int main() {
    std::printf("=== L2: Chiral 3NF Recoupling (scalar + rank-2) ===\n\n");

    std::printf("--- Scalar recoupling ---\n");
    test_scalar_3S1_diagonal();
    test_scalar_1S0_diagonal();
    test_scalar_3S1_3D1_zero();
    test_scalar_spectator_mismatch_zero();

    std::printf("\n--- Rank-2 recoupling ---\n");
    test_rank2_3S1_3D1_nonzero();
    test_rank2_singlet_zero();
    test_rank2_T_mismatch_zero();
    test_rank2_spectator_mismatch_zero();
    test_rank2_dL_1_zero();
    test_rank2_3S1_3D1_value();

    std::printf("\n%d passed, %d failed\n", g_passes, g_failures);
    return g_failures > 0 ? 1 : 0;
}
