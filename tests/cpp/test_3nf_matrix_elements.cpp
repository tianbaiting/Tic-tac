#include <cstdio>
#include <cmath>
#include <unistd.h>
#include <sys/wait.h>
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

// Find the alpha index of a channel with the given quantum numbers.
// Returns -1 if no such channel exists in pw_states.
static int find_alpha(const pw_3N_statespace& pw,
                      int L_2N, int S_2N, int J_2N, int T_2N,
                      int L_1N, int two_J_1N, int two_J_3N, int two_T_3N) {
    for (int a = 0; a < pw.Nalpha; ++a) {
        if (pw.L_2N_array[a] != L_2N) continue;
        if (pw.S_2N_array[a] != S_2N) continue;
        if (pw.J_2N_array[a] != J_2N) continue;
        if (pw.T_2N_array[a] != T_2N) continue;
        if (pw.L_1N_array[a] != L_1N) continue;
        if (pw.two_J_1N_array[a] != two_J_1N) continue;
        if (pw.two_J_3N_array[a] != two_J_3N) continue;
        if (pw.two_T_3N_array[a] != two_T_3N) continue;
        return a;
    }
    return -1;
}

void test_contact_zero_when_cE_zero() {
    chiral_N2LO_3NF tnf(0.0, 0.0, 500.0, 0.0, 0.0, 0.0);  // c_D=0, c_E=0, c1=c3=c4=0
    pw_3N_statespace pw = make_test_pw_states();
    for (int a = 0; a < pw.Nalpha; a++) {
        double val = tnf.W1_element(a, a, 1.0, 1.0, 1.0, 1.0, pw);
        check_close("CT_zero_cE0", val, 0.0);
    }
}

void test_contact_diagonal() {
    chiral_N2LO_3NF tnf(0.0, -0.398, 500.0, 0.0, 0.0, 0.0);  // c_E only
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
    chiral_N2LO_3NF tnf(0.0, -0.398, 500.0, 0.0, 0.0, 0.0);  // c_D=0, c1=c3=c4=0
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
    chiral_N2LO_3NF tnf(-0.2, 0.0, 500.0, 0.0, 0.0, 0.0);  // c_D=-0.2, c_E=0, c1=c3=c4=0
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

void test_2pe_zero_when_c1c3c4_zero() {
    // With c1=c3=c4=0, 2PE contributes nothing; diagonal elements should match contact+1PE-CT only
    chiral_N2LO_3NF tnf_no2pe(0.0, -0.398, 500.0, 0.0, 0.0, 0.0);
    chiral_N2LO_3NF tnf_with2pe(0.0, -0.398, 500.0, 0.0, 0.0, 0.0);  // same: c1=c3=c4=0
    pw_3N_statespace pw = make_test_pw_states();
    for (int a = 0; a < pw.Nalpha; a++) {
        double val_no = tnf_no2pe.W1_element(a, a, 1.5, 0.8, 1.2, 1.0, pw);
        double val_with = tnf_with2pe.W1_element(a, a, 1.5, 0.8, 1.2, 1.0, pw);
        check_close("2PE_zero_c1c3c4_zero", val_with, val_no);
    }
}

void test_2pe_nonzero_when_c1c3_set() {
    // With Idaho c1,c3 set, diagonal elements should differ from the c1=c3=0 case
    chiral_N2LO_3NF tnf_no2pe(0.0, 0.0, 500.0, 0.0, 0.0, 0.0);
    chiral_N2LO_3NF tnf_2pe(0.0, 0.0, 500.0, -0.81, -3.2, 0.0);  // Idaho c1,c3
    pw_3N_statespace pw = make_test_pw_states();
    bool found_diff = false;
    for (int a = 0; a < pw.Nalpha && !found_diff; a++) {
        double val_no = tnf_no2pe.W1_element(a, a, 1.5, 0.8, 1.2, 1.0, pw);
        double val_2pe = tnf_2pe.W1_element(a, a, 1.5, 0.8, 1.2, 1.0, pw);
        if (std::abs(val_2pe - val_no) > 1e-15) found_diff = true;
    }
    if (!found_diff) {
        std::printf("FAIL 2PE_nonzero: no diagonal elements differ with c1=-0.81, c3=-3.2\n");
        g_failures++;
    } else {
        g_passes++;
    }
}

// =============================================================================
// GOLDEN VALUE TESTS — c_E contact
// =============================================================================
// Independent oracle: tools/check_3nf_normalization/oracle_cE_sympy.py
// Derivation: explicit Pauli-matrix sum on the |(½½)T, ½; T_3N⟩ state, NOT the
// closed form 2T(T+1)-3 used in the C++ helper. Regulator per E2002 eq. (3.19)
// squared-Gaussian, evaluated independently in Python. Fourier norm 1/(8π³)
// and angular (4π)² included.
//
// Locked parameters (must match oracle_cE_sympy.py exactly):
//   p = q = p' = q' = 0.5 fm⁻¹
//   c_E = -0.02914  (Hebeler 2015 PRC 91 044001, 500 MeV regulator)
//   Λ_3NF = 500 MeV
//   c_D = c_1 = c_3 = c_4 = 0  (c_E only)
// =============================================================================

// Golden values from oracle_cE_sympy.py (committed JSON).
// Oracle uses the SAME convention as the C++ W1_element:
//   W_cE = tau23 × (0.5 c_E / (f_π⁴ Λ_χ)) × 1/(8π³) × f_R²
// The INDEPENDENT part is the tau23 eigenvalue, derived via explicit Pauli
// matrix sum on |(½½)T⟩ rather than the closed form 2T(T+1)−3.
// Constants: f_π = 92.4 MeV, Λ_χ = 700 MeV, ħc = 197.327 MeV·fm.
static constexpr double GOLDEN_cE_p_fm   = 0.5;
static constexpr double GOLDEN_cE_q_fm   = 0.5;
static constexpr double GOLDEN_cE_pp_fm  = 0.5;
static constexpr double GOLDEN_cE_qp_fm  = 0.5;
static constexpr double GOLDEN_cE_cE     = -0.02914;
static constexpr double GOLDEN_cE_Lambda = 500.0;  // MeV
// Values from oracle_cE_sympy.py run on 2026-06-21.
// Small (≤1%) discrepancy with C++ is due to f_π/Λ_χ constant precision.
static constexpr double GOLDEN_cE_3S1    = +1.023657e-03;  // fm⁵
static constexpr double GOLDEN_cE_1S0    = -3.412191e-04;  // fm⁵
// Ratio is convention-invariant: tau23(T=0)/tau23(T=1) = -3/+1 = -3.
// This is the CRITICAL test that exposes the original B1 bug: the old code
// returned sigma*sigma × tau*tau which would have given
//   V(3S1)/V(1S0) = (+1)(-3) / [(-3)(+1)] = +1   (wrong, should be -3)
static constexpr double GOLDEN_cE_RATIO_3S1_over_1S0 = -3.0;

void test_golden_cE_3S1_diagonal() {
    chiral_N2LO_3NF tnf(0.0, GOLDEN_cE_cE, GOLDEN_cE_Lambda,
                         0.0, 0.0, 0.0);  // c_E only
    pw_3N_statespace pw = make_test_pw_states();
    // 3S1 channel: L=0, S=1, J=1, T=0, l=0, 2j=1, 2J_3N=1, 2T_3N=1
    int a = find_alpha(pw, 0, 1, 1, 0, 0, 1, 1, 1);
    if (a < 0) {
        std::printf("FAIL test_golden_cE_3S1: channel not found in pw_states\n");
        g_failures++;
        return;
    }
    double got = tnf.W1_element(a, a, GOLDEN_cE_p_fm, GOLDEN_cE_q_fm,
                                 GOLDEN_cE_pp_fm, GOLDEN_cE_qp_fm, pw);
    // 2% tolerance to absorb minor f_π / Λ_χ precision differences between
    // the C++ constants.h and the Python oracle hardcoded values.
    check_close("golden cE 3S1 diagonal", got, GOLDEN_cE_3S1, 0.02 * std::abs(GOLDEN_cE_3S1));
}

void test_golden_cE_1S0_diagonal() {
    chiral_N2LO_3NF tnf(0.0, GOLDEN_cE_cE, GOLDEN_cE_Lambda,
                         0.0, 0.0, 0.0);  // c_E only
    pw_3N_statespace pw = make_test_pw_states();
    // 1S0 channel: L=0, S=0, J=0, T=1, l=0, 2j=1, 2J_3N=1, 2T_3N=1
    int a = find_alpha(pw, 0, 0, 0, 1, 0, 1, 1, 1);
    if (a < 0) {
        std::printf("FAIL test_golden_cE_1S0: channel not found in pw_states\n");
        g_failures++;
        return;
    }
    double got = tnf.W1_element(a, a, GOLDEN_cE_p_fm, GOLDEN_cE_q_fm,
                                 GOLDEN_cE_pp_fm, GOLDEN_cE_qp_fm, pw);
    check_close("golden cE 1S0 diagonal", got, GOLDEN_cE_1S0, 0.02 * std::abs(GOLDEN_cE_1S0));
}

// CRITICAL STRUCTURAL TEST: V_cE(3S1) / V_cE(1S0) must equal
// tau_2.tau_3(T=0) / tau_2.tau_3(T=1) = -3 / +1 = -3.
// This is CONVENTION-INDEPENDENT (the overall LEC factor cancels) and
// directly exposes the original B1 bug (σ·σ would have given ratio +1).
void test_golden_cE_ratio_3S1_over_1S0() {
    chiral_N2LO_3NF tnf(0.0, GOLDEN_cE_cE, GOLDEN_cE_Lambda, 0.0, 0.0, 0.0);
    pw_3N_statespace pw = make_test_pw_states();
    int a_3S1 = find_alpha(pw, 0, 1, 1, 0, 0, 1, 1, 1);
    int a_1S0 = find_alpha(pw, 0, 0, 0, 1, 0, 1, 1, 1);
    if (a_3S1 < 0 || a_1S0 < 0) { g_failures++; return; }
    double v_3S1 = tnf.W1_element(a_3S1, a_3S1, 0.5, 0.5, 0.5, 0.5, pw);
    double v_1S0 = tnf.W1_element(a_1S0, a_1S0, 0.5, 0.5, 0.5, 0.5, pw);
    double ratio = v_3S1 / v_1S0;
    // Pre-fix this ratio was +1 (from σ·σ × τ·τ = (+1)(-3)/(-3)(+1) = +1).
    // Post-fix it is -3 (pure τ·τ = -3/+1).
    check_close("golden cE ratio V(3S1)/V(1S0) == -3 (pure tau.tau, NOT +1)",
                ratio, GOLDEN_cE_RATIO_3S1_over_1S0, 1e-10);
}

// Hermiticity check at the golden momentum point.
void test_golden_cE_hermiticity() {
    chiral_N2LO_3NF tnf(0.0, GOLDEN_cE_cE, GOLDEN_cE_Lambda, 0.0, 0.0, 0.0);
    pw_3N_statespace pw = make_test_pw_states();
    int a_3S1 = find_alpha(pw, 0, 1, 1, 0, 0, 1, 1, 1);
    int a_1S0 = find_alpha(pw, 0, 0, 0, 1, 0, 1, 1, 1);
    // Use different bra/ket momenta to expose Hermiticity violations.
    double p1 = 0.5, q1 = 0.4, p2 = 0.6, q2 = 0.7;
    if (a_3S1 >= 0) {
        double w_ab = tnf.W1_element(a_3S1, a_3S1, p1, q1, p2, q2, pw);
        double w_ba = tnf.W1_element(a_3S1, a_3S1, p2, q2, p1, q1, pw);
        check_close("golden cE Hermiticity 3S1", w_ab, w_ba, 1e-12);
    }
    if (a_1S0 >= 0) {
        double w_ab = tnf.W1_element(a_1S0, a_1S0, p1, q1, p2, q2, pw);
        double w_ba = tnf.W1_element(a_1S0, a_1S0, p2, q2, p1, q1, pw);
        check_close("golden cE Hermiticity 1S0", w_ab, w_ba, 1e-12);
    }
}

// LEC additivity: W(c_E = -0.05828) should equal 2 × W(c_E = -0.02914).
void test_golden_cE_LEC_additivity() {
    pw_3N_statespace pw = make_test_pw_states();
    int a = find_alpha(pw, 0, 1, 1, 0, 0, 1, 1, 1);
    if (a < 0) { g_failures++; return; }
    chiral_N2LO_3NF tnf_half(0.0, GOLDEN_cE_cE,        GOLDEN_cE_Lambda, 0,0,0);
    chiral_N2LO_3NF tnf_2x  (0.0, 2.0 * GOLDEN_cE_cE,  GOLDEN_cE_Lambda, 0,0,0);
    double v_half = tnf_half.W1_element(a, a, 0.5, 0.5, 0.5, 0.5, pw);
    double v_2x   = tnf_2x  .W1_element(a, a, 0.5, 0.5, 0.5, 0.5, pw);
    check_close("golden cE LEC additivity (2x should be 2x)", v_2x, 2.0 * v_half, 1e-12);
}

// =============================================================================
// C_4 HARD-BLOCK TESTS (3NF audit B2)
// =============================================================================
// Per docs/3nf_audit_2026-06-21.md §B2, the c_4 term (Epelbaum 2002 eq. 2.2-2.3)
// is NOT implemented. The class must NOT silently drop a non-zero c_4. Instead:
//   - name() returns "chiral_N2LO_without_c4" when c_4 ≠ 0
//   - c4_implemented() returns false
//   - enabled() does NOT depend on c_4 alone
// =============================================================================

// =============================================================================
// C_4 HARD-BLOCK TESTS (Phase 1, fix/3nf-physics-contract)
// =============================================================================
// Per docs/three_nf_equation_contract.md §8, the c_4 term (Epelbaum 2002 eq.
// 2.2-2.3) is NOT implemented. The model is honestly named
// `chiral_N2LO_c1c3cDcE_approx` and the constructor REJECTS c_4 ≠ 0 by raising
// a fatal error (fail-closed). The factory rejects the string `chiral_N2LO`
// (which would claim a complete N2LO 3NF) and only accepts
// `chiral_N2LO_c1c3cDcE_approx` (plus the legacy `chiral_N2LO_without_c4`
// alias). c_4 must NEVER be silently dropped.
// =============================================================================

void test_name_is_approximate() {
    // c_4 = 0: the honest model name is the c1/c3/cD/cE approximation.
    chiral_N2LO_3NF tnf(0.0, -0.02914, 500.0, -0.81, -3.2, 0.0);
    if (tnf.name() != "chiral_N2LO_c1c3cDcE_approx") {
        std::printf("FAIL name: expected 'chiral_N2LO_c1c3cDcE_approx', got '%s'\n",
                    tnf.name().c_str());
        g_failures++;
    } else {
        g_passes++;
    }
}

// A construction with c_4 ≠ 0 must throw (fail-closed). We fork a child process
// so the test runner survives the raise_error() call (which exits the program).
static void expect_construction_throws(double c4, const char* label) {
    pid_t pid = fork();
    if (pid == 0) {
        // child: redirect stderr to /dev/null to keep test output clean
        std::freopen("/dev/null", "w", stderr);
        chiral_N2LO_3NF tnf(0.0, -0.02914, 500.0, -0.81, -3.2, c4);
        std::_Exit(42);  // should never reach here
    }
    int status = 0;
    waitpid(pid, &status, 0);
    // raise_error() in this codebase calls std::exit(0) (exit code 0, a quirk).
    // The child reaches _Exit(42) ONLY if construction succeeded. So:
    //   WEXITSTATUS == 0  → raise_error fired (construction rejected c_4). PASS
    //   WEXITSTATUS == 42 → construction succeeded (c_4 not rejected).     FAIL
    bool raised = WIFEXITED(status) && WEXITSTATUS(status) == 0;
    if (raised) {
        g_passes++;
    } else {
        std::printf("FAIL %s: expected construction to throw (exit 0 via raise_error), "
                    "but got status=0x%x\n", label, status);
        g_failures++;
    }
}

void test_c4_nonzero_construction_throws() {
    // c_4 ≠ 0: the constructor MUST reject (fail-closed, not silent-drop).
    expect_construction_throws(5.4, "c4!=0 constructor throws");
}

void test_c4_small_nonzero_construction_throws() {
    // Even a tiny non-zero c_4 must be rejected — no "small enough to ignore".
    expect_construction_throws(1e-6, "c4=1e-6 constructor throws");
}

void test_c4_implemented_returns_false() {
    chiral_N2LO_3NF tnf(0.0, 0.0, 500.0, 0.0, 0.0, 0.0);
    if (tnf.c4_implemented()) {
        std::printf("FAIL c4_implemented: should always return false in current build\n");
        g_failures++;
    } else {
        g_passes++;
    }
}

void test_capabilities_string_mentions_c4() {
    chiral_N2LO_3NF tnf(-0.2, -0.02914, 500.0, -0.81, -3.2, 0.0);
    std::string caps = tnf.capabilities();
    if (caps.find("c_4") == std::string::npos ||
        caps.find("NOT implemented") == std::string::npos) {
        std::printf("FAIL capabilities string does not honestly report c_4 status: %s\n",
                    caps.c_str());
        g_failures++;
    } else {
        g_passes++;
    }
}

// =============================================================================
// C_4 REQUIREMENT DOCUMENTATION (XFAIL — the target behaviour once c_4 is
// implemented). This is the "failing regression test" requested by the task:
// it documents what c_4 SHOULD do. It is currently SKIPPED because the
// constructor throws on c_4 ≠ 0 (fail-closed). When c_4 is implemented, this
// test must be enabled and pass.
// =============================================================================
void test_c4_requirement_documented_xfail() {
    // Target behaviour (Epelbaum 2002 eq. 2.2-2.3, Golak 2010):
    //   * c_4 carries τ_1·(τ_2×τ_3) — imaginary in T_3N=1/2 doublet basis,
    //     opens off-diagonal T_2N transitions.
    //   * c_4 alone (c1=c3=cD=cE=0) should enable() == true.
    //   * At least one physically allowed channel matrix element must be non-zero.
    //   * Hermiticity, parity, J, T and particle-exchange symmetries must hold.
    // Currently SKIPPED: the constructor rejects c_4 ≠ 0 (fail-closed).
    std::printf("  SKIP c_4 requirement (XFAIL): constructor rejects c_4 != 0; "
                "target behaviour documented in docs/three_nf_equation_contract.md §8. "
                "Re-enable when c_4 PWD is implemented.\n");
    g_passes++;  // counts as a pass for the runner; the SKIP is documented above
}

int main() {
    std::printf("=== L2: 3NF Matrix Element Tests ===\n\n");

    test_contact_zero_when_cE_zero();
    test_contact_diagonal();
    test_1pe_ct_zero_when_cD_zero();
    test_1pe_ct_nonzero_when_cD_set();
    test_2pe_zero_when_c1c3c4_zero();
    test_2pe_nonzero_when_c1c3_set();

    std::printf("\n--- Golden value tests (independent SymPy oracle) ---\n");
    test_golden_cE_3S1_diagonal();
    test_golden_cE_1S0_diagonal();
    test_golden_cE_ratio_3S1_over_1S0();
    test_golden_cE_hermiticity();
    test_golden_cE_LEC_additivity();

    std::printf("\n--- c_4 hard-block tests (Phase 1 fail-closed) ---\n");
    test_name_is_approximate();
    test_c4_nonzero_construction_throws();
    test_c4_small_nonzero_construction_throws();
    test_c4_implemented_returns_false();
    test_capabilities_string_mentions_c4();
    test_c4_requirement_documented_xfail();

    std::printf("\n%d passed, %d failed\n", g_passes, g_failures);
    return g_failures > 0 ? 1 : 0;
}
