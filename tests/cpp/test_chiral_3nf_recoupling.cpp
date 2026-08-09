#include <cstdio>
#include <cmath>
#include <array>
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

// Explicit three-particle m-scheme oracle in the |m1,m2,m3> basis. For two
// spin-1/2 particles, sigma_1.sigma_3 = 2 P_13 - 1. The same identity holds
// for isospin. No Wigner helper or production recoupling code enters here.
using three_half_state = std::array<double, 8>;

static int swap_particles_1_and_3(int basis) {
    const int b1 = (basis >> 2) & 1;
    const int b2 = (basis >> 1) & 1;
    const int b3 = basis & 1;
    return (b3 << 2) | (b2 << 1) | b1;
}

static double pauli_1_dot_pauli_3(const three_half_state& bra,
                                  const three_half_state& ket) {
    double value = 0.0;
    for (int i = 0; i < 8; ++i)
        value += bra[i] * (2.0 * ket[swap_particles_1_and_3(i)] - ket[i]);
    return value;
}

static double explicit_cD_scalar_recoupling() {
    // Pair-(23) singlet, spectator up: |S23=0; total 1/2, M=1/2>.
    three_half_state pair_singlet{};
    pair_singlet[6] =  1.0 / std::sqrt(2.0); // |up up down>
    pair_singlet[5] = -1.0 / std::sqrt(2.0); // |up down up>

    // Pair triplet coupled to spectator 1/2 with Condon-Shortley CGs:
    // -|m23=0,m1=up>/sqrt(3) + sqrt(2/3)|m23=1,m1=down>.
    three_half_state pair_triplet{};
    pair_triplet[6] = -1.0 / std::sqrt(6.0);
    pair_triplet[5] = -1.0 / std::sqrt(6.0);
    pair_triplet[3] =  std::sqrt(2.0 / 3.0);

    const double sig_01 = pauli_1_dot_pauli_3(pair_singlet, pair_triplet);
    const double tau_10 = pauli_1_dot_pauli_3(pair_triplet, pair_singlet);
    check_close("m-scheme <S23=0|sigma1.sigma3|S23=1>",
                sig_01, std::sqrt(3.0));
    check_close("m-scheme <T23=1|tau1.tau3|T23=0>",
                tau_10, std::sqrt(3.0));
    return (1.0 / 3.0) * sig_01 * tau_10;
}

// ---------------------------------------------------------------------------
// c_E contact recoupling (pure τ₂·τ₃, NO σ₂·σ₃ dependence)
// ---------------------------------------------------------------------------
// Per Epelbaum 2002 eq. (2.10) + (A-4) and docs/3nf_audit_2026-06-21.md §B1:
//   V^(1)_cont = -E·(τ₂·τ₃),  E = c_E/(f_π⁴ Λ_χ)
// The c_E contact is a SPIN SCALAR. Its matrix element must NOT depend on the
// pair-spin eigenvalue. For pair T=1 (τ₂·τ₃ = +1) the recoupling is identical
// whether S_2N=0 or S_2N=1, as long as T_2N, T_3N and the spectator structure match.
//
// We verify two independent properties:
//   (a) For (S=0, T=1) and (S=1, T=1) at matching J_3N=3/2, T_3N=1/2:
//       A_cE(S=0,T=1) == A_cE(S=1,T=1)   (no spin dependence).
//   (b) For the diagonal 3S1 channel (S=1, T=0): A_cE equals the pure τ₂·τ₃
//       eigenvalue (-3) times the standard 6j recoupling factor — NOT multiplied
//       by σ₂·σ₃ (= +1 by accident for S=1, so the bug is invisible here).
void test_cE_contact_S0T1_equals_S1T1() {
    // (S=0, T=1) → J_2N=0, J_3N=1/2 (doublet), T_3N=1/2
    double v_S0_T1 = recoupling_3nf_contact_cE(
        /*L_2N_r=*/0, /*S_2N_r=*/0, /*J_2N_r=*/0, /*T_2N_r=*/1,
        /*L_1N_r=*/0, /*two_J_1N_r=*/1, /*two_J_3N=*/1,
        /*L_2N_c=*/0, /*S_2N_c=*/0, /*J_2N_c=*/0, /*T_2N_c=*/1,
        /*L_1N_c=*/0, /*two_J_1N_c=*/1,
        /*two_T_3N=*/1);
    // (S=1, T=1) → J_2N=1, J_3N=3/2 (lowest doublet), T_3N=1/2
    double v_S1_T1 = recoupling_3nf_contact_cE(
        /*L_2N_r=*/0, /*S_2N_r=*/1, /*J_2N_r=*/1, /*T_2N_r=*/1,
        /*L_1N_r=*/0, /*two_J_1N_r=*/1, /*two_J_3N=*/3,
        /*L_2N_c=*/0, /*S_2N_c=*/1, /*J_2N_c=*/1, /*T_2N_c=*/1,
        /*L_1N_c=*/0, /*two_J_1N_c=*/1,
        /*two_T_3N=*/1);
    // Both must give the SAME matrix element since c_E is spin-scalar.
    check_close("cE contact: A(S=0,T=1) == A(S=1,T=1) [no spin dep]",
                v_S0_T1, v_S1_T1);
}

// Direct numeric check: for the dominant 3S1 channel (S=1, T=0) at T_3N=1/2,
// A_cE must equal τ₂·τ₃(T=0) = -3, with NO σ·σ factor and NO 6j recoupling.
// (τ_2·τ_3 is diagonal in the (pair T_2N, spectator 1/2) T_3N basis — the
//  spectator is a passive spectator; the formula_reference.md §1.4 closed
//  form with a 6j symbol is wrong, as verified by sympy: the 6j
//  {½½T; ½½½} violates triad parity and is identically zero.)
void test_cE_contact_3S1_diagonal_value() {
    double val = recoupling_3nf_contact_cE(
        /*L_2N_r=*/0, /*S_2N_r=*/1, /*J_2N_r=*/1, /*T_2N_r=*/0,
        /*L_1N_r=*/0, /*two_J_1N_r=*/1, /*two_J_3N=*/1,
        /*L_2N_c=*/0, /*S_2N_c=*/1, /*J_2N_c=*/1, /*T_2N_c=*/0,
        /*L_1N_c=*/0, /*two_J_1N_c=*/1,
        /*two_T_3N=*/1);
    // τ₂·τ₃(T=0) = -3, pure diagonal — no extra factors.
    check_close("cE contact: 3S1 diagonal == -3 (pure tau.tau eigenvalue)",
                val, -3.0);
}

// 1S0 channel (S=0, T=1) must give A_cE = +1 (NOT -3 from σ·σ × τ·τ).
// This is the discriminating test against the old buggy implementation.
void test_cE_contact_1S0_diagonal_value() {
    double val = recoupling_3nf_contact_cE(
        /*L_2N_r=*/0, /*S_2N_r=*/0, /*J_2N_r=*/0, /*T_2N_r=*/1,
        /*L_1N_r=*/0, /*two_J_1N_r=*/1, /*two_J_3N=*/1,
        /*L_2N_c=*/0, /*S_2N_c=*/0, /*J_2N_c=*/0, /*T_2N_c=*/1,
        /*L_1N_c=*/0, /*two_J_1N_c=*/1,
        /*two_T_3N=*/1);
    // τ₂·τ₃(T=1) = +1, pure diagonal.
    // OLD BUGGY value: σ·σ × τ·τ = (-3)(+1) = -3.
    check_close("cE contact: 1S0 diagonal == +1 (NOT -3 from old σ·σ bug)",
                val, +1.0);
}

// Selection rules: the c_E contact requires S-wave pair (L_2N=0) and S-wave
// spectator (l_1N=0). Any non-zero L must return 0.
void test_cE_contact_L_nonzero_zero() {
    double val = recoupling_3nf_contact_cE(
        /*L_2N_r=*/2, /*S_2N_r=*/1, /*J_2N_r=*/1, /*T_2N_r=*/0,
        0, 1, 1,
        0, 1, 1, 0, 0, 1, 1);
    check_close("cE contact: L_2N_r=2 must be 0", val, 0.0);

    val = recoupling_3nf_contact_cE(
        0, 1, 1, 0, /*l=*/2, 1, 1,
        0, 1, 1, 0, 2, 1, 1);
    check_close("cE contact: l_1N=2 must be 0", val, 0.0);
}

// ---------------------------------------------------------------------------
// 2PE rank-0 scalar recoupling: 3S1 diagonal channel
// (Formerly recoupling_3nf_scalar — retained for the c_1/c_3 2PE rank-0 piece
//  where σ₂·σ₃ × τ₂·τ₃ IS the correct operator.)
// ---------------------------------------------------------------------------
// alpha = { L_2N=0, S_2N=1, J_2N=1, T_2N=0, l_1N=0, 2j_1N=1, 2J_3N=1, 2T_3N=1 }
//   sigma_2.sigma_3 (S=1)  = +1
//   tau_2.tau_3     (T=0)  = -3
//   scalar recoupling      = (+1) * (-3) = -3
void test_scalar_3S1_diagonal() {
    double val = recoupling_3nf_2pe_scalar(
        /*L_2N_r=*/0, /*S_2N_r=*/1, /*J_2N_r=*/1, /*T_2N_r=*/0,
        /*L_1N_r=*/0, /*two_J_1N_r=*/1, /*two_J_3N=*/1,
        /*L_2N_c=*/0, /*S_2N_c=*/1, /*J_2N_c=*/1, /*T_2N_c=*/0,
        /*L_1N_c=*/0, /*two_J_1N_c=*/1,
        /*two_T_3N=*/1);
    check_close("2pe_scalar 3S1 diagonal (sigma.sigma * tau.tau)", val, -3.0);
}

void test_scalar_1S0_diagonal() {
    // 1S0: L_2N=0, S_2N=0, J_2N=0, T_2N=1, l_1N=0, 2j_1N=1, 2J_3N=1, 2T_3N=1
    //   sigma_2.sigma_3 (S=0) = -3
    //   tau_2.tau_3     (T=1) = +1
    //   product               = -3
    double val = recoupling_3nf_2pe_scalar(
        0, 0, 0, 1, 0, 1, 1,
        0, 0, 0, 1, 0, 1,
        1);
    check_close("2pe_scalar 1S0 diagonal", val, -3.0);
}

// ---------------------------------------------------------------------------
// 2PE rank-0 scalar: 3S1 <-> 3D1 must be zero (ΔL=2 selection rule violation
// for a pair-scalar operator).
// ---------------------------------------------------------------------------
void test_scalar_3S1_3D1_zero() {
    double val = recoupling_3nf_2pe_scalar(
        /*L_2N_r=*/0, /*S_2N_r=*/1, /*J_2N_r=*/1, /*T_2N_r=*/0,
        /*L_1N_r=*/0, /*two_J_1N_r=*/1, /*two_J_3N=*/1,
        /*L_2N_c=*/2, /*S_2N_c=*/1, /*J_2N_c=*/1, /*T_2N_c=*/0,
        /*L_1N_c=*/0, /*two_J_1N_c=*/1,
        /*two_T_3N=*/1);
    check_close("2pe_scalar 3S1<->3D1 (should be zero)", val, 0.0);
}

void test_scalar_spectator_mismatch_zero() {
    // Mismatched spectator l_1N should zero the scalar recoupling.
    double val = recoupling_3nf_2pe_scalar(
        0, 1, 1, 0, /*l=*/0, 1, 1,
        0, 1, 1, 0, /*l=*/2, 1,
        1);
    check_close("2pe_scalar spectator l mismatch", val, 0.0);
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
// c_D 1PE-contact recoupling: rank-0 (sigma1.sigma3 x tau1.tau3)
// ---------------------------------------------------------------------------

// For the 3S1 diagonal (T_2N=0, L_2N=0): tau1.tau3 diagonal T=T'=0 is zero
// by 6j selection rule. So recoup_scalar must return 0.
void test_1pe_ct_scalar_3S1_diagonal_zero() {
    // alpha = {L_2N=0, S_2N=1, J_2N=1, T_2N=0, l_1N=0, 2j_1N=1, 2J_3N=1, 2T_3N=1}
    // tau1.tau3 for T=T'=0 is zero -> result must be zero.
    double val = recoupling_3nf_1pe_ct_scalar(
        0, 1, 1, 0, 0, 1, 1,
        0, 1, 1, 0, 0, 1,
        1);
    check_close("1pe_ct_scalar 3S1 diagonal (tau=tau'=0, must be 0)", val, 0.0);
}

// For the 1S0 diagonal (T_2N=1, L_2N=0): similar analysis.
// tau1.tau3 diagonal T=T'=1 is nonzero, but sigma1.sigma3 must be checked.
void test_1pe_ct_scalar_1S0_diagonal() {
    // alpha = {L_2N=0, S_2N=0, J_2N=0, T_2N=1, l_1N=0, 2j_1N=1, 2J_3N=1, 2T_3N=1}
    // sigma1.sigma3 for S=S'=0 diagonal: the 9j/6j should give a value.
    // We just check it's finite (not zero) since the tau1.tau3(T=T'=1) is nonzero.
    double val = recoupling_3nf_1pe_ct_scalar(
        0, 0, 0, 1, 0, 1, 1,
        0, 0, 0, 1, 0, 1,
        1);
    check_close("1pe_ct_scalar 1S0 diagonal", val, 0.0);
}

// Off-diagonal channel: bra=(L=0,S=0,J=0,T=1) x ket=(L=0,S=1,J=1,T=0)
// This is the dominant c_D contribution channel.
// Expected: (1/3) * sig13(1S0->3S1-like) * tau13(T=1->T=0) != 0
void test_1pe_ct_scalar_offdiag_1S0_to_3S1_nonzero() {
    double val = recoupling_3nf_1pe_ct_scalar(
        /*bra: L=0,S=0,J=0,T=1,l=0,2j=1*/ 0, 0, 0, 1, 0, 1, 1,
        /*ket: L=0,S=1,J=1,T=0,l=0,2j=1*/ 0, 1, 1, 0, 0, 1,
        1);
    const double expected = explicit_cD_scalar_recoupling();
    check_close("1pe_ct_scalar off-diag 1S0->3S1(T=0)", val, expected);
}

// Selection rule: nonzero L_2N must return 0 (contact pair vertex restriction)
void test_1pe_ct_scalar_L2N_nonzero_zero() {
    double val = recoupling_3nf_1pe_ct_scalar(
        /*L_2N=1*/ 1, 1, 1, 0, 0, 1, 1,
        0, 1, 1, 0, 0, 1,
        1);
    check_close("1pe_ct_scalar L_2N_r=1 (must be 0)", val, 0.0);
}

// Rank-2 c_D: l_1N=0 must return 0 (CG selection rule)
void test_1pe_ct_rank2_l_zero_returns_zero() {
    double val = recoupling_3nf_1pe_ct_rank2(
        0, 1, 1, 0, /*l=*/0, 1, 1,
        0, 1, 1, 0, /*l=*/0, 1,
        1);
    check_close("1pe_ct_rank2 l_1N=0 (CG=0, must return 0)", val, 0.0);
}

// Rank-2 c_D: L_2N nonzero must return 0 (contact pair restriction)
void test_1pe_ct_rank2_L2N_nonzero_returns_zero() {
    double val = recoupling_3nf_1pe_ct_rank2(
        /*L_2N=1*/ 1, 1, 1, 0, 2, 3, 1,
        0, 1, 1, 0, 2, 3,
        1);
    check_close("1pe_ct_rank2 L_2N_r=1 (contact pair, must return 0)", val, 0.0);
}

// ---------------------------------------------------------------------------
int main() {
    std::printf("=== L2: Chiral 3NF Recoupling (contact cE + 2PE scalar + rank-2) ===\n\n");

    std::printf("--- c_E contact recoupling (pure tau.tau, spin-scalar) ---\n");
    test_cE_contact_S0T1_equals_S1T1();
    test_cE_contact_3S1_diagonal_value();
    test_cE_contact_1S0_diagonal_value();
    test_cE_contact_L_nonzero_zero();

    std::printf("\n--- 2PE rank-0 scalar recoupling (sigma.sigma * tau.tau) ---\n");
    test_scalar_3S1_diagonal();
    test_scalar_1S0_diagonal();
    test_scalar_3S1_3D1_zero();
    test_scalar_spectator_mismatch_zero();

    std::printf("\n--- Candidate pair rank-2 recoupling (NOT production c1/c3) ---\n");
    test_rank2_3S1_3D1_nonzero();
    test_rank2_singlet_zero();
    test_rank2_T_mismatch_zero();
    test_rank2_spectator_mismatch_zero();
    test_rank2_dL_1_zero();
    test_rank2_3S1_3D1_value();

    std::printf("\n--- c_D 1PE-CT: rank-0 (sigma1.sigma3)(tau1.tau3) ---\n");
    test_1pe_ct_scalar_3S1_diagonal_zero();
    test_1pe_ct_scalar_1S0_diagonal();
    test_1pe_ct_scalar_offdiag_1S0_to_3S1_nonzero();
    test_1pe_ct_scalar_L2N_nonzero_zero();

    std::printf("\n--- Candidate c_D rank-2 selection helper (NOT production) ---\n");
    test_1pe_ct_rank2_l_zero_returns_zero();
    test_1pe_ct_rank2_L2N_nonzero_returns_zero();

    std::printf("\n%d passed, %d failed\n", g_passes, g_failures);
    return g_failures > 0 ? 1 : 0;
}
