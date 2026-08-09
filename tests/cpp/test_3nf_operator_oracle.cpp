// tests/cpp/test_3nf_operator_oracle.cpp
//
// Phase 0 acceptance test (see docs/three_nf_equation_contract.md §6).
//
// Finite-dimensional, operator-level oracle that exercises the REAL production
// kernel builders `calculate_CPVC_col` and `calculate_all_CPVC_rows` with
// hand-constructed, deliberately NON-symmetric mock P, V, W^(1), C and checks
// that every code path returns the same matrix
//
//     M = C^T · ( P·V  +  W^(1)·(1 + P) ) · C
//
// assembled directly from the dense definitions. The mock operators are
// independent of the production 3NF recoupling / kernels: they are random
// dense matrices with a controlled structure (C and V contain genuine
// alpha-off-diagonal 3S1--3D1-like blocks, P is a generic permutation).
// Because W^(1), C, and the off-diagonal C_01/C_10 blocks are asymmetric,
// coupled-alpha contraction, transposition, and ordering bugs cannot hide
// behind symmetry.
//
// This is the gate for Phase 3 (independent angular oracle): until every check
// here passes, no claim about partial-wave matrix elements is admissible.

#include <cstdio>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

#include "constants.h"
#include "type_defs.h"
#include "make_pw_symm_states.h"
#include "solve_faddeev.h"
#include "coupled_channel_transform.h"
#include "pade_approximant.h"
#include "three_nucleon_force_model.h"
#include "w1_pw_cache.h"

namespace {

constexpr int N_ALPHA = 2;
constexpr int N_P     = 2;
constexpr int N_Q     = 2;
constexpr int DIM     = N_ALPHA * N_P * N_Q;

constexpr double TOL = 1e-10;

static int g_failures = 0;
static int g_passes   = 0;

static void check_close(const char* label, double got, double expected,
                        double tol = TOL) {
    if (std::abs(got - expected) > tol) {
        std::printf("FAIL %s: got %.10e, expected %.10e (diff %.3e)\n",
                    label, got, expected, std::abs(got - expected));
        g_failures++;
    } else {
        g_passes++;
    }
}

// Flattened (α, p, q) → linear index, matching solve_faddeev.cpp's convention
// (α outermost, then q, then p innermost — see e.g. line 362 of solve_faddeev.cpp:
//  idx_row = idx_alpha_r * Nq_WP * Np_WP + idx_q_r * Np_WP + idx_p_r).
inline int idx_apq(int a, int p, int q) {
    return a * N_P * N_Q + q * N_P + p;
}

// -----------------------------------------------------------------------------
// Mock WP grid. The midpoints and widths are used both by the production code
// (which multiplies W^(1) by p_mid·√dp·q_mid·√dq) and by the mock W1_element
// (which divides them back out, so the algebra is tested free of unit noise).
// -----------------------------------------------------------------------------
struct MockGrid {
    double p_WP[N_P + 1];
    double q_WP[N_Q + 1];
    double p_mid[N_P];
    double q_mid[N_Q];
    double dp[N_P];
    double dq[N_Q];
    void init() {
        // Non-uniform bins so midpoints/widths are distinct and the test can
        // detect index-mixing bugs.
        p_WP[0] = 1.0;  p_WP[1] = 3.0;  p_WP[2] = 7.0;
        q_WP[0] = 2.0;  q_WP[1] = 5.0;  q_WP[2] = 9.0;
        for (int i = 0; i < N_P; ++i) {
            p_mid[i] = 0.5 * (p_WP[i] + p_WP[i + 1]);
            dp[i]    = p_WP[i + 1] - p_WP[i];
        }
        for (int i = 0; i < N_Q; ++i) {
            q_mid[i] = 0.5 * (q_WP[i] + q_WP[i + 1]);
            dq[i]    = q_WP[i + 1] - q_WP[i];
        }
    }
    // Map an fm⁻¹ momentum back to its bin index.
    int p_bin(double p_fm) const {
        double p_MeV = p_fm * hbarc;
        for (int i = 0; i < N_P; ++i)
            if (p_MeV >= p_WP[i] && p_MeV < p_WP[i + 1]) return i;
        return N_P - 1;
    }
    int q_bin(double q_fm) const {
        double q_MeV = q_fm * hbarc;
        for (int i = 0; i < N_Q; ++i)
            if (q_MeV >= q_WP[i] && q_MeV < q_WP[i + 1]) return i;
        return N_Q - 1;
    }
};

static MockGrid g_grid;

// -----------------------------------------------------------------------------
// Mock W^(1) table: a fully general 6-index array W1[α_r][α_c][p_r][q_r][p_c][q_c]
// (NOT separable, so rank-1 bugs cannot hide). The mock three_nucleon_force_model
// returns  W1_table / (w1_unit·wp_r·wq_r·wp_c·wq_c)  so that the production
// code's  w1_bin = w1_raw · w1_unit · (wp_r·wq_r·wp_c·wq_c)  == W1_table.
// This isolates the ALGEBRA from the WP-normalisation bookkeeping (Phase 4).
// -----------------------------------------------------------------------------
struct W1Table {
    double v[N_ALPHA][N_ALPHA][N_P][N_Q][N_P][N_Q];
    void fill(std::mt19937_64& rng) {
        std::uniform_real_distribution<double> d(-2.0, 2.0);
        for (int ar = 0; ar < N_ALPHA; ++ar)
        for (int ac = 0; ac < N_ALPHA; ++ac)
        for (int pr = 0; pr < N_P; ++pr)
        for (int qr = 0; qr < N_Q; ++qr)
        for (int pc = 0; pc < N_P; ++pc)
        for (int qc = 0; qc < N_Q; ++qc)
            v[ar][ac][pr][qr][pc][qc] = d(rng) + 0.3;  // biased away from zero
    }
};

static W1Table g_w1;

class MockTNF : public three_nucleon_force_model {
public:
    bool enabled() const override { return true; }
    std::string name() const override { return "mock_oracle_W1"; }

    // Returns W1_table / (w1_unit·wp·wq·wp·wq) so the production code's
    // normalised w1_bin equals W1_table exactly.
    double W1_element(int alpha_r, int alpha_c,
                      double p_r_fm, double q_r_fm,
                      double p_c_fm, double q_c_fm,
                      const pw_3N_statespace& /*pw*/) const override {
        const int pr = g_grid.p_bin(p_r_fm);
        const int qr = g_grid.q_bin(q_r_fm);
        const int pc = g_grid.p_bin(p_c_fm);
        const int qc = g_grid.q_bin(q_c_fm);
        const double wp_r = g_grid.p_mid[pr] * std::sqrt(g_grid.dp[pr]);
        const double wq_r = g_grid.q_mid[qr] * std::sqrt(g_grid.dq[qr]);
        const double wp_c = g_grid.p_mid[pc] * std::sqrt(g_grid.dp[pc]);
        const double wq_c = g_grid.q_mid[qc] * std::sqrt(g_grid.dq[qc]);
        const double inv_hbarc  = 1.0 / hbarc;
        const double inv_hbarc5 = inv_hbarc*inv_hbarc*inv_hbarc*inv_hbarc*inv_hbarc;
        const double w1_unit = inv_hbarc5;  // w1_scale == 1
        return g_w1.v[alpha_r][alpha_c][pr][qr][pc][qc]
               / (w1_unit * wp_r * wq_r * wp_c * wq_c);
    }
};

class ZeroTNF : public three_nucleon_force_model {
public:
    bool enabled() const override { return true; }
    std::string name() const override { return "zero_oracle_W1"; }
    double W1_element(int, int, double, double, double, double,
                      const pw_3N_statespace&) const override { return 0.0; }
};

// -----------------------------------------------------------------------------
// Build the dense mock operators. C and V act only at fixed q, but are full in
// the two alpha channels (3S1/3D1-like Np×Np blocks). P is a full DIM×DIM
// permutation matrix.
// -----------------------------------------------------------------------------
struct DenseMock {
    std::vector<double> C;    // DIM×DIM, q-diagonal and full in alpha
    std::vector<double> V;    // DIM×DIM, q-diagonal and full in alpha
    std::vector<double> P;    // DIM×DIM permutation
    std::vector<double> W1;   // DIM×DIM = g_w1 with normalisation folded in
    std::vector<double> CT;   // C^T
};

static std::vector<double> transpose(const std::vector<double>& M) {
    std::vector<double> T(DIM * DIM);
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            T[i * DIM + j] = M[j * DIM + i];
    return T;
}

static std::vector<double> mat_mat(const std::vector<double>& A,
                                   const std::vector<double>& B) {
    std::vector<double> C(DIM * DIM, 0.0);
    for (int i = 0; i < DIM; ++i)
        for (int k = 0; k < DIM; ++k) {
            const double a = A[i * DIM + k];
            if (a == 0.0) continue;
            for (int j = 0; j < DIM; ++j)
                C[i * DIM + j] += a * B[k * DIM + j];
        }
    return C;
}

static std::vector<double> mat_add(const std::vector<double>& A,
                                   const std::vector<double>& B) {
    std::vector<double> C(DIM * DIM);
    for (int i = 0; i < DIM * DIM; ++i) C[i] = A[i] + B[i];
    return C;
}

static DenseMock build_dense_mock(std::mt19937_64& rng) {
    DenseMock m;
    m.C.assign(DIM * DIM, 0.0);
    m.V.assign(DIM * DIM, 0.0);
    m.P.assign(DIM * DIM, 0.0);
    m.W1.assign(DIM * DIM, 0.0);
    std::uniform_real_distribution<double> d(-1.0, 1.0);

    // C and V are q-independent pair operators.  Unlike the old oracle, every
    // (alpha_r,alpha_c) block is populated.  The alpha-dependent offset makes
    // C_01 and C_10 explicitly nonzero and non-transposes of one another; the
    // random term makes every p block asymmetric as well.
    for (int ar = 0; ar < N_ALPHA; ++ar)
    for (int ac = 0; ac < N_ALPHA; ++ac) {
        std::vector<double> Cblk(N_P * N_P), Vblk(N_P * N_P);
        for (int pr = 0; pr < N_P; ++pr)
        for (int pc = 0; pc < N_P; ++pc) {
            const double diagonal = (ar == ac && pr == pc) ? 1.0 : 0.0;
            Cblk[pr * N_P + pc] = diagonal
                                 + 0.07 * (1 + 2 * ar + 3 * ac)
                                 + 0.15 * d(rng);
            Vblk[pr * N_P + pc] = 0.35
                                 + 0.05 * (2 * ar + ac)
                                 + 0.25 * d(rng);
        }
        // Replicate the pair blocks across q.
        for (int q = 0; q < N_Q; ++q)
        for (int pr = 0; pr < N_P; ++pr)
        for (int pc = 0; pc < N_P; ++pc) {
            int i = idx_apq(ar, pr, q);
            int j = idx_apq(ac, pc, q);
            m.C[i * DIM + j] = Cblk[pr * N_P + pc];
            m.V[i * DIM + j] = Vblk[pr * N_P + pc];
        }
    }

    // P: a non-trivial permutation on the full (α,p,q) index. Use a 4-cycle
    // mixing α and p so P ≠ P^T (to expose transpose bugs) and P ≠ I.
    // 4-cycle: 0 → 1 → 2 → 3 → 0  on the first four (α,p,q) states (q=0).
    {
        int cycle[] = {idx_apq(0,0,0), idx_apq(0,1,0), idx_apq(1,0,0), idx_apq(1,1,0)};
        for (int t = 0; t < 4; ++t) {
            int from = cycle[t];
            int to   = cycle[(t + 1) % 4];
            m.P[to * DIM + from] = 1.0;  // P[to, from] = 1  →  (P·x)[to] = x[from]
        }
        // identity elsewhere
        for (int i = 0; i < DIM; ++i) {
            bool placed = false;
            for (int j = 0; j < DIM; ++j) if (m.P[i*DIM+j] != 0.0) { placed = true; break; }
            if (!placed) m.P[i * DIM + i] = 1.0;
        }
    }

    // W1: full DIM×DIM from g_w1 with the production normalisation folded in.
    // W1_dense[α_r,p_r,q_r ; α_c,p_c,q_c] = g_w1[...] (already post-normalisation,
    // because the mock W1_element divided the normalisation out).
    for (int ar = 0; ar < N_ALPHA; ++ar)
    for (int ac = 0; ac < N_ALPHA; ++ac)
    for (int pr = 0; pr < N_P; ++pr)
    for (int qr = 0; qr < N_Q; ++qr)
    for (int pc = 0; pc < N_P; ++pc)
    for (int qc = 0; qc < N_Q; ++qc)
        m.W1[idx_apq(ar,pr,qr) * DIM + idx_apq(ac,pc,qc)] = g_w1.v[ar][ac][pr][qr][pc][qc];

    m.CT = transpose(m.C);
    return m;
}

// The mathematical oracle:  M = C^T·(P·V + W1·(1+P))·C
static std::vector<double> build_M_math(const DenseMock& m) {
    auto PV  = mat_mat(m.P, m.V);
    auto W1P  = mat_mat(m.W1, m.P);
    auto W1_one_plus_P = mat_add(m.W1, W1P);          // W1·(1+P)
    auto K    = mat_add(PV, W1_one_plus_P);           // P·V + W1·(1+P)
    auto KC  = mat_mat(K, m.C);                        // K·C
    return mat_mat(m.CT, KC);                          // C^T·K·C
}

// -----------------------------------------------------------------------------
// Build the mock WP environment (CT_RM_array, VC_CM_array, P123 sparse, pw_states)
// from the dense mock operators, in the exact format the production code expects.
// -----------------------------------------------------------------------------
struct MockWP {
    pw_3N_statespace pw;
    int L_2N[N_ALPHA], S_2N[N_ALPHA], J_2N[N_ALPHA], T_2N[N_ALPHA];
    int L_1N[N_ALPHA], two_J_1N[N_ALPHA], two_J_3N[N_ALPHA], two_T_3N[N_ALPHA], P_3N[N_ALPHA];
    int chn_3N_idx[N_ALPHA];

    // CT_RM_array[α·Nα + α']  = Np×Np block (row-major), or nullptr.
    // Holds C^T restricted to the pair line within (alpha,alpha').  All four
    // blocks are populated. VC_CM_array similarly stores V*C (column-major).
    double* CT_blocks[N_ALPHA * N_ALPHA];  // owning
    double* VC_blocks[N_ALPHA * N_ALPHA];
    std::vector<double> CT_storage;       // backing store
    std::vector<double> VC_storage;
    double** CT_RM_array;
    double** VC_CM_array;

    // P123 CSC sparse.  P_code = 2 · P123_stored, so we store P123 = P/2.
    std::vector<double> P_val;
    std::vector<int>    P_row;
    std::vector<size_t> P_col;            // size DIM+1
    size_t P_dim;

    MockWP(const DenseMock& m) {
        // pw_states: all α share the same conserved quantum numbers so the
        // production code's conservation checks pass for all (α_r, α_c) pairs.
        pw.Nalpha = N_ALPHA;
        pw.J_2N_max = 1;
        for (int a = 0; a < N_ALPHA; ++a) {
            L_2N[a] = 2 * a;  // alpha 0 = 3S1-like, alpha 1 = 3D1-like
            S_2N[a] = 1; J_2N[a] = 1; T_2N[a] = 0;
            L_1N[a] = 0; two_J_1N[a] = 1; two_J_3N[a] = 1; two_T_3N[a] = 1; P_3N[a] = +1;
            chn_3N_idx[a] = 0;
        }
        pw.L_2N_array = L_2N; pw.S_2N_array = S_2N; pw.J_2N_array = J_2N; pw.T_2N_array = T_2N;
        pw.L_1N_array = L_1N; pw.two_J_1N_array = two_J_1N;
        pw.two_J_3N_array = two_J_3N; pw.two_T_3N_array = two_T_3N; pw.P_3N_array = P_3N;
        pw.chn_3N_idx_array = chn_3N_idx; pw.N_chn_3N = 1;

        // For a requested right column (alpha_c,p_c), the production builders
        // index both pointer tables as [alpha_c*Nalpha + alpha_j]:
        //
        // CT_RM[alpha_c,alpha_j][p_c,p_j]
        //   = (C^T)[alpha_c p_c,alpha_j p_j]
        //   = C[alpha_j p_j,alpha_c p_c].
        //
        // VC_CM[alpha_c,alpha_j][p_c,p_j]
        //   = (V*C)[alpha_j p_j,alpha_c p_c].
        //
        // The first index is therefore the external column channel and the
        // second is the intermediate row channel, despite the pointer-array
        // storage being described as C^T row-major.
        CT_RM_array = new double*[N_ALPHA * N_ALPHA];
        VC_CM_array = new double*[N_ALPHA * N_ALPHA];
        CT_storage.assign((size_t)N_ALPHA * N_ALPHA * N_P * N_P, 0.0);
        VC_storage.assign((size_t)N_ALPHA * N_ALPHA * N_P * N_P, 0.0);
        const int q = 0;
        for (int alpha_c = 0; alpha_c < N_ALPHA; ++alpha_c)
        for (int alpha_j = 0; alpha_j < N_ALPHA; ++alpha_j) {
            const size_t bi = (size_t)alpha_c * N_ALPHA + alpha_j;
            double* ct = &CT_storage[bi * N_P * N_P];
            double* vc = &VC_storage[bi * N_P * N_P];
            for (int p_c = 0; p_c < N_P; ++p_c)
            for (int p_j = 0; p_j < N_P; ++p_j) {
                ct[(size_t)p_c * N_P + p_j]
                    = m.C[idx_apq(alpha_j, p_j, q) * DIM
                          + idx_apq(alpha_c, p_c, q)];

                double vc_element = 0.0;
                for (int alpha_k = 0; alpha_k < N_ALPHA; ++alpha_k)
                for (int p_k = 0; p_k < N_P; ++p_k) {
                    const double V_jk
                        = m.V[idx_apq(alpha_j, p_j, q) * DIM
                              + idx_apq(alpha_k, p_k, q)];
                    const double C_kc
                        = m.C[idx_apq(alpha_k, p_k, q) * DIM
                              + idx_apq(alpha_c, p_c, q)];
                    vc_element += V_jk * C_kc;
                }
                vc[(size_t)p_c * N_P + p_j] = vc_element;
            }
            CT_RM_array[bi] = ct;
            VC_CM_array[bi] = vc;
        }

        // P123 CSC sparse: store P/2 so production 2·P123 = P (mock).
        // P is DIM×DIM; build CSC over columns.
        P_col.assign(DIM + 1, 0);
        std::vector<std::vector<std::pair<int,double>>> rows_per_col(DIM);
        for (int j = 0; j < DIM; ++j)
        for (int i = 0; i < DIM; ++i)
            if (m.P[i * DIM + j] != 0.0)
                rows_per_col[j].push_back({i, m.P[i*DIM+j] / 2.0});
        size_t nnz = 0;
        for (int j = 0; j < DIM; ++j) {
            P_col[j] = nnz;
            for (auto& pr : rows_per_col[j]) {
                P_row.push_back(pr.first);
                P_val.push_back(pr.second);
                nnz++;
            }
        }
        P_col[DIM] = nnz;
        P_dim = nnz;
    }
    ~MockWP() {
        delete[] CT_RM_array;
        delete[] VC_CM_array;
    }
};

}  // namespace

// =============================================================================
// Test 1: pointer-table orientation, independent of the kernel contraction.
//
// CT_RM has one convention, regardless of how its indices are named at a call
// site:
//   CT_RM[a,b][i,j] = (C^T)[a i,b j] = C[b j,a i].
// The right-C application names (a,b,i,j)=(alpha_c,alpha_j,p_c,p_j), while the
// final left-C^T application names them (alpha_r,alpha_i,p_r,p_i).  These are
// two uses of the same layout, not two pointer-table conventions.
// =============================================================================
static void test_pointer_table_layout(const DenseMock& m) {
    MockWP wp(m);
    const auto VC = mat_mat(m.V, m.C);

    int ct_mismatches = 0;
    int right_c_mismatches = 0;
    int vc_mismatches = 0;
    double ct_max_diff = 0.0;
    double right_c_max_diff = 0.0;
    double vc_max_diff = 0.0;
    for (int a = 0; a < N_ALPHA; ++a)
    for (int b = 0; b < N_ALPHA; ++b)
    for (int i = 0; i < N_P; ++i)
    for (int j = 0; j < N_P; ++j) {
        const double stored_ct = wp.CT_RM_array[a * N_ALPHA + b][i * N_P + j];
        const double expected_ct
            = m.CT[idx_apq(a, i, 0) * DIM + idx_apq(b, j, 0)];
        const double expected_right_c
            = m.C[idx_apq(b, j, 0) * DIM + idx_apq(a, i, 0)];
        const double ct_diff = std::abs(stored_ct - expected_ct);
        const double right_c_diff = std::abs(stored_ct - expected_right_c);
        if (ct_diff > TOL) {
            ct_mismatches++;
            ct_max_diff = std::max(ct_max_diff, ct_diff);
        }
        if (right_c_diff > TOL) {
            right_c_mismatches++;
            right_c_max_diff = std::max(right_c_max_diff, right_c_diff);
        }

        // VC_CM[a,b][i,j] is column-major storage of
        // (V*C)[row=(b,j), col=(a,i)].
        const double stored_vc = wp.VC_CM_array[a * N_ALPHA + b][i * N_P + j];
        const double expected_vc
            = VC[idx_apq(b, j, 0) * DIM + idx_apq(a, i, 0)];
        const double vc_diff = std::abs(stored_vc - expected_vc);
        if (vc_diff > TOL) {
            vc_mismatches++;
            vc_max_diff = std::max(vc_max_diff, vc_diff);
        }
    }

    check_close("CT_RM[a,b][i,j] == C^T[ai,bj]: #mismatches",
                (double)ct_mismatches, 0.0, 0.5);
    check_close("CT_RM layout: max|diff|", ct_max_diff, 0.0);
    check_close("right C from CT_RM[a,b][i,j] == C[bj,ai]: #mismatches",
                (double)right_c_mismatches, 0.0, 0.5);
    check_close("right C orientation: max|diff|", right_c_max_diff, 0.0);
    check_close("VC_CM[a,b][i,j] == (V*C)[bj,ai]: #mismatches",
                (double)vc_mismatches, 0.0, 0.5);
    check_close("VC_CM layout: max|diff|", vc_max_diff, 0.0);
}

// =============================================================================
// Test 2: the production coupled-channel block transform for an actual
// 3S1--3D1 pair obeys CT_RM[a,b][i,j] = C[b,j; a,i].  The deliberately
// non-symmetric full C matrix makes the wrong, untransposed orientation fail.
// =============================================================================
static void test_real_3S1_3D1_CT_storage() {
    constexpr int pair_dim = 2 * N_P;
    const int L[N_ALPHA] = {0, 2};
    constexpr int J = 1;

    std::vector<double> C(pair_dim * pair_dim);
    for (int row = 0; row < pair_dim; ++row)
    for (int col = 0; col < pair_dim; ++col)
        C[row * pair_dim + col] = 10.0 * (row + 1) + 0.7 * (col + 1)
                                      + 0.03 * row * col;

    std::vector<double> CT_blocks = C;
    for (int row = 0; row < pair_dim; ++row)
    for (int col = row + 1; col < pair_dim; ++col)
        std::swap(CT_blocks[row * pair_dim + col],
                  CT_blocks[col * pair_dim + row]);
    restructure_coupled_matrix_blocks(CT_blocks.data(), N_P);

    int mismatches = 0;
    int wrong_orientation_coincidences = 0;
    double max_diff = 0.0;
    for (int alpha_r = 0; alpha_r < N_ALPHA; ++alpha_r)
    for (int alpha_c = 0; alpha_c < N_ALPHA; ++alpha_c) {
        const bool row_is_D = L[alpha_r] > J;
        const bool col_is_D = L[alpha_c] > J;
        const std::size_t block
            = coupled_matrix_block_index(row_is_D, col_is_D);
        const double* stored = CT_blocks.data() + block * N_P * N_P;
        for (int p_r = 0; p_r < N_P; ++p_r)
        for (int p_c = 0; p_c < N_P; ++p_c) {
            const int full_r = alpha_r * N_P + p_r;
            const int full_c = alpha_c * N_P + p_c;
            const double got = stored[p_r * N_P + p_c];
            const double expected = C[full_c * pair_dim + full_r];
            const double wrong = C[full_r * pair_dim + full_c];
            const double diff = std::abs(got - expected);
            if (diff > TOL) {
                mismatches++;
                max_diff = std::max(max_diff, diff);
            }
            if (full_r != full_c && std::abs(got - wrong) <= TOL)
                wrong_orientation_coincidences++;
        }
    }

    check_close("real 3S1-3D1 CT storage: #mismatches",
                static_cast<double>(mismatches), 0.0, 0.5);
    check_close("real 3S1-3D1 CT storage: max|diff|", max_diff, 0.0);
    check_close("real 3S1-3D1 rejects C-for-CT orientation",
                static_cast<double>(wrong_orientation_coincidences), 0.0, 0.5);
}

// =============================================================================
// Test 2: isolate C^T*W1*C with coupled alpha blocks.
//
// P and V are set to zero in the production call, so a failure cannot be
// attributed to the sparse permutation convention or to VC construction.  The
// old implementation drops alpha_j != alpha_c from W1*C and must fail here.
// =============================================================================
static void test_coupled_alpha_W1C(const DenseMock& m) {
    MockWP wp(m);
    MockTNF tnf;
    tnf_kernel_context ctx{};
    ctx.tnf = &tnf;
    ctx.pw_states = &wp.pw;
    ctx.p_WP_array = g_grid.p_WP;
    ctx.q_WP_array = g_grid.q_WP;
    ctx.CT_RM_array = wp.CT_RM_array;
    ctx.w1_scale = 1.0;
    ctx.w1_cache = nullptr;

    const auto W1C = mat_mat(m.W1, m.C);
    const auto expected = mat_mat(m.CT, W1C);

    std::vector<double> zero_vc_storage(
        (size_t)N_ALPHA * N_ALPHA * N_P * N_P, 0.0);
    std::vector<double*> zero_vc_blocks(N_ALPHA * N_ALPHA, nullptr);
    for (int b = 0; b < N_ALPHA * N_ALPHA; ++b)
        zero_vc_blocks[b] = &zero_vc_storage[(size_t)b * N_P * N_P];

    // Empty CSC matrix: every column starts and ends at offset zero.
    std::vector<size_t> zero_p_col(DIM + 1, 0);
    double dummy_p_val = 0.0;
    int dummy_p_row = 0;

    int mismatch_count = 0;
    double max_diff = 0.0;
    for (int alpha_c = 0; alpha_c < N_ALPHA; ++alpha_c)
    for (int q_c = 0; q_c < N_Q; ++q_c)
    for (int p_c = 0; p_c < N_P; ++p_c) {
        std::vector<double> col(DIM, 0.0);
        std::vector<int> row_to_nnz(DIM, -1);
        std::vector<int> nnz_to_row(DIM, -1);
        size_t num_nnz = 0;
        calculate_CPVC_col(col.data(), row_to_nnz.data(), nnz_to_row.data(), num_nnz,
                           alpha_c, p_c, q_c, N_ALPHA, N_Q, N_P,
                           wp.CT_RM_array, zero_vc_blocks.data(),
                           &dummy_p_val, &dummy_p_row, zero_p_col.data(), 0, ctx);
        const int col_idx = idx_apq(alpha_c, p_c, q_c);
        for (int row = 0; row < DIM; ++row) {
            const double diff = std::abs(col[row] - expected[row * DIM + col_idx]);
            if (diff > TOL) {
                mismatch_count++;
                max_diff = std::max(max_diff, diff);
            }
        }
    }
    check_close("coupled-alpha C^T*W1*C: #mismatches",
                (double)mismatch_count, 0.0, 0.5);
    check_close("coupled-alpha C^T*W1*C: max|diff|",
                max_diff, 0.0);
}

// =============================================================================
// Test 3: calculate_CPVC_col reproduces M_math column-by-column.
// =============================================================================
static void test_CPVC_col_matches_math(const DenseMock& m, const std::vector<double>& M_math) {
    MockWP wp(m);
    MockTNF tnf;
    tnf_kernel_context ctx{};
    ctx.tnf = &tnf;
    ctx.pw_states = &wp.pw;
    ctx.p_WP_array = g_grid.p_WP;
    ctx.q_WP_array = g_grid.q_WP;
    ctx.CT_RM_array = wp.CT_RM_array;
    ctx.w1_scale = 1.0;
    ctx.w1_cache = nullptr;

    int max_diff_count = 0;
    double max_diff = 0.0;
    for (int alpha_c = 0; alpha_c < N_ALPHA; ++alpha_c)
    for (int q_c = 0; q_c < N_Q; ++q_c)
    for (int p_c = 0; p_c < N_P; ++p_c) {
        std::vector<double> col(DIM, 0.0);
        std::vector<int> row_to_nnz(DIM, -1);
        std::vector<int> nnz_to_row(DIM, -1);
        size_t num_nnz = 0;
        calculate_CPVC_col(col.data(), row_to_nnz.data(), nnz_to_row.data(), num_nnz,
                           alpha_c, p_c, q_c, N_ALPHA, N_Q, N_P,
                           wp.CT_RM_array, wp.VC_CM_array,
                           wp.P_val.data(), wp.P_row.data(), wp.P_col.data(), wp.P_dim, ctx);
        int col_idx = idx_apq(alpha_c, p_c, q_c);
        for (int i = 0; i < DIM; ++i) {
            double diff = std::abs(col[i] - M_math[i * DIM + col_idx]);
            if (diff > TOL) {
                max_diff_count++;
                if (diff > max_diff) max_diff = diff;
            }
        }
    }
    check_close("CPVC_col: #mismatches", (double)max_diff_count, 0.0, 0.5);
    check_close("CPVC_col: max|diff| vs M_math", max_diff, 0.0);
}

// =============================================================================
// Test 4: production row path equals both the production column path and the
// independent dense matrix for every selected elastic row and every column.
// =============================================================================
static void test_CPVC_rows_share_algebra(const DenseMock& m,
                                         const std::vector<double>& M_math) {
    MockWP wp(m);
    MockTNF tnf;
    tnf_kernel_context ctx{};
    ctx.tnf = &tnf;
    ctx.pw_states = &wp.pw;
    ctx.p_WP_array = g_grid.p_WP;
    ctx.q_WP_array = g_grid.q_WP;
    ctx.CT_RM_array = wp.CT_RM_array;
    ctx.w1_scale = 1.0;
    ctx.w1_cache = nullptr;

    int q_com_idx[N_Q] = {0, 1};
    int deuteron_idx[N_ALPHA] = {0, 1};
    constexpr int N_SELECTED_ROWS = N_ALPHA * N_Q;
    std::vector<double> rows(N_SELECTED_ROWS * DIM, 0.0);
    calculate_all_CPVC_rows(rows.data(),
                            q_com_idx, N_Q,
                            deuteron_idx, N_ALPHA,
                            N_ALPHA, N_Q, N_P,
                            wp.CT_RM_array, wp.VC_CM_array,
                            wp.P_val.data(), wp.P_row.data(), wp.P_col.data(), wp.P_dim,
                            ctx);

    int math_mismatches = 0;
    int column_mismatches = 0;
    double max_math_diff = 0.0;
    double max_column_diff = 0.0;
    for (int alpha_r = 0; alpha_r < N_ALPHA; ++alpha_r)
    for (int q_r = 0; q_r < N_Q; ++q_r) {
        const int selected_row = alpha_r * N_Q + q_r;
        const int dense_row = idx_apq(alpha_r, 0, q_r);
        for (int alpha_c = 0; alpha_c < N_ALPHA; ++alpha_c)
        for (int q_c = 0; q_c < N_Q; ++q_c)
        for (int p_c = 0; p_c < N_P; ++p_c) {
            const int col_idx = idx_apq(alpha_c, p_c, q_c);
            const double row_value = rows[selected_row * DIM + col_idx];
            const double math_diff = std::abs(row_value - M_math[dense_row * DIM + col_idx]);
            if (math_diff > TOL) {
                math_mismatches++;
                max_math_diff = std::max(max_math_diff, math_diff);
            }

            std::vector<double> col(DIM, 0.0);
            std::vector<int> row_to_nnz(DIM, -1);
            std::vector<int> nnz_to_row(DIM, -1);
            size_t num_nnz = 0;
            calculate_CPVC_col(col.data(), row_to_nnz.data(), nnz_to_row.data(), num_nnz,
                               alpha_c, p_c, q_c, N_ALPHA, N_Q, N_P,
                               wp.CT_RM_array, wp.VC_CM_array,
                               wp.P_val.data(), wp.P_row.data(), wp.P_col.data(), wp.P_dim,
                               ctx);
            const double column_diff = std::abs(row_value - col[dense_row]);
            if (column_diff > TOL) {
                column_mismatches++;
                max_column_diff = std::max(max_column_diff, column_diff);
            }
        }
    }

    check_close("CPVC_rows: #mismatches vs dense math", (double)math_mismatches, 0.0, 0.5);
    check_close("CPVC_rows: max|diff| vs dense math", max_math_diff, 0.0);
    check_close("CPVC_rows: #mismatches vs CPVC_col", (double)column_mismatches, 0.0, 0.5);
    check_close("CPVC_rows: max|diff| vs CPVC_col", max_column_diff, 0.0);
}

// =============================================================================
// Test 5: all disabled/zero-strength 3NF routes reproduce the pure-2NF kernel.
// This covers a null context, the diagnostic w1_scale=0 path, and an enabled
// model whose W1_element is identically zero.
// =============================================================================
static void test_zero_W1_reproduces_2NF(const DenseMock& m) {
    MockWP wp(m);
    MockTNF mock_tnf;
    ZeroTNF zero_tnf;

    const auto two_nf_math = mat_mat(m.CT, mat_mat(mat_mat(m.P, m.V), m.C));
    tnf_kernel_context contexts[3]{};
    contexts[0].tnf = nullptr;
    contexts[0].w1_scale = 1.0;
    contexts[1].tnf = &mock_tnf;
    contexts[1].w1_scale = 0.0;
    contexts[2].tnf = &zero_tnf;
    contexts[2].w1_scale = 1.0;
    for (auto& ctx : contexts) {
        ctx.pw_states = &wp.pw;
        ctx.p_WP_array = g_grid.p_WP;
        ctx.q_WP_array = g_grid.q_WP;
        ctx.CT_RM_array = wp.CT_RM_array;
        ctx.w1_cache = nullptr;
    }

    int mismatches = 0;
    double max_diff = 0.0;
    for (const auto& ctx : contexts)
    for (int alpha_c = 0; alpha_c < N_ALPHA; ++alpha_c)
    for (int q_c = 0; q_c < N_Q; ++q_c)
    for (int p_c = 0; p_c < N_P; ++p_c) {
        std::vector<double> col(DIM, 0.0);
        std::vector<int> row_to_nnz(DIM, -1);
        std::vector<int> nnz_to_row(DIM, -1);
        size_t num_nnz = 0;
        calculate_CPVC_col(col.data(), row_to_nnz.data(), nnz_to_row.data(), num_nnz,
                           alpha_c, p_c, q_c, N_ALPHA, N_Q, N_P,
                           wp.CT_RM_array, wp.VC_CM_array,
                           wp.P_val.data(), wp.P_row.data(), wp.P_col.data(), wp.P_dim,
                           ctx);
        const int col_idx = idx_apq(alpha_c, p_c, q_c);
        for (int row = 0; row < DIM; ++row) {
            const double diff = std::abs(col[row] - two_nf_math[row * DIM + col_idx]);
            if (diff > TOL) {
                mismatches++;
                max_diff = std::max(max_diff, diff);
            }
        }
    }
    check_close("W1=0/disabled: #mismatches vs 2NF", (double)mismatches, 0.0, 0.5);
    check_close("W1=0/disabled: max|diff| vs 2NF", max_diff, 0.0);
}

// =============================================================================
// Test 6: the real W1_PW_cache consumer path and the direct W1_element fallback
// agree element-by-element. Np=Nq=1 quadrature is intentional here: it makes
// the cache and fallback evaluate the identical midpoint matrix element, so
// this test isolates cache indexing (including the intermediate alpha_j).
// =============================================================================
static void test_W1_cache_matches_fallback(const DenseMock& m) {
    MockWP wp(m);
    MockTNF tnf;
    run_params params{};
    params.Np_WP = N_P;
    params.Nq_WP = N_Q;
    params.Np_per_WP_W1 = 1;
    params.Nq_per_WP_W1 = 1;
    params.potential_model = "mock_2NF";
    params.three_nucleon_force = tnf.name();
    params.c_D = 0.0;
    params.c_E = 0.0;
    params.Lambda_3NF = 500.0;

    W1_PW_cache cache;
    cache.build(tnf, g_grid.p_WP, N_P, g_grid.q_WP, N_Q, wp.pw, params);

    tnf_kernel_context direct_ctx{};
    direct_ctx.tnf = &tnf;
    direct_ctx.pw_states = &wp.pw;
    direct_ctx.p_WP_array = g_grid.p_WP;
    direct_ctx.q_WP_array = g_grid.q_WP;
    direct_ctx.CT_RM_array = wp.CT_RM_array;
    direct_ctx.w1_scale = 0.37;
    direct_ctx.w1_cache = nullptr;
    tnf_kernel_context cache_ctx = direct_ctx;
    cache_ctx.w1_cache = &cache;

    int mismatches = 0;
    double max_diff = 0.0;
    for (int alpha_c = 0; alpha_c < N_ALPHA; ++alpha_c)
    for (int q_c = 0; q_c < N_Q; ++q_c)
    for (int p_c = 0; p_c < N_P; ++p_c) {
        std::vector<double> direct_col(DIM, 0.0), cached_col(DIM, 0.0);
        std::vector<int> row_to_nnz(DIM, -1), nnz_to_row(DIM, -1);
        size_t direct_nnz = 0;
        calculate_CPVC_col(direct_col.data(), row_to_nnz.data(), nnz_to_row.data(), direct_nnz,
                           alpha_c, p_c, q_c, N_ALPHA, N_Q, N_P,
                           wp.CT_RM_array, wp.VC_CM_array,
                           wp.P_val.data(), wp.P_row.data(), wp.P_col.data(), wp.P_dim,
                           direct_ctx);
        std::fill(row_to_nnz.begin(), row_to_nnz.end(), -1);
        std::fill(nnz_to_row.begin(), nnz_to_row.end(), -1);
        size_t cached_nnz = 0;
        calculate_CPVC_col(cached_col.data(), row_to_nnz.data(), nnz_to_row.data(), cached_nnz,
                           alpha_c, p_c, q_c, N_ALPHA, N_Q, N_P,
                           wp.CT_RM_array, wp.VC_CM_array,
                           wp.P_val.data(), wp.P_row.data(), wp.P_col.data(), wp.P_dim,
                           cache_ctx);
        for (int row = 0; row < DIM; ++row) {
            const double diff = std::abs(cached_col[row] - direct_col[row]);
            if (diff > TOL) {
                mismatches++;
                max_diff = std::max(max_diff, diff);
            }
        }
    }
    check_close("W1 cache/fallback: #mismatches", (double)mismatches, 0.0, 0.5);
    check_close("W1 cache/fallback: max|diff|", max_diff, 0.0);
}

// =============================================================================
// Test 7: Operator-ordering guard — wrong order (1+P)·W1 ≠ right order W1·(1+P).
// =============================================================================
static void test_operator_ordering(const DenseMock& m) {
    auto PV      = mat_mat(m.P, m.V);
    auto PW1     = mat_mat(m.P, m.W1);          // P·W1 (wrong order)
    auto W1C     = mat_mat(m.W1, m.C);
    auto wrong_K = mat_add(PV, mat_add(PW1, m.W1));   // P·V + P·W1 + W1  (wrong)
    auto wrong   = mat_mat(m.CT, mat_mat(wrong_K, m.C));

    auto right = build_M_math(m);

    double max_diff = 0.0;
    for (int i = 0; i < DIM * DIM; ++i)
        max_diff = std::max(max_diff, std::abs(wrong[i] - right[i]));
    // wrong ≠ right  →  max_diff must be > 0  (this guards against a regression
    // where the two orderings accidentally coincide).
    if (max_diff > 1e-6) g_passes++;
    else { g_failures++; std::printf("FAIL operator-ordering: wrong==right (max_diff=%.3e)\n", max_diff); }
}

// =============================================================================
// Test 6: C vs C^T guard — using C in place of C^T must NOT match M_math.
// =============================================================================
static void test_CT_distinction(const DenseMock& m, const std::vector<double>& M_math) {
    // Build M_wrongC = C·(P·V + W1·(1+P))·C^T  (C and C^T swapped).
    auto PV      = mat_mat(m.P, m.V);
    auto W1P     = mat_mat(m.W1, m.P);
    auto K       = mat_add(PV, mat_add(m.W1, W1P));
    auto KCT     = mat_mat(K, m.CT);
    auto M_wrong = mat_mat(m.C, KCT);   // C·K·C^T  (swapped)

    double max_diff = 0.0;
    for (int i = 0; i < DIM * DIM; ++i)
        max_diff = std::max(max_diff, std::abs(M_wrong[i] - M_math[i]));
    if (max_diff > 1e-6) g_passes++;
    else { g_failures++; std::printf("FAIL CT-distinction: C==C^T (max_diff=%.3e)\n", max_diff); }
}

// =============================================================================
// Test 7: Dense-solver identity — verify U = (1 − A·G)⁻¹·A satisfies (1-AG)·U = A.
// Uses a clean Gauss-Jordan solve with full pivoting on the augmented matrix.
// =============================================================================
static void test_dense_solver_identity(const DenseMock& m, const std::vector<double>& A) {
    // G = diag(g_i); real g_i for a non-singular (1 − A·G).
    std::vector<double> G(DIM, 0.0);
    for (int i = 0; i < DIM; ++i) G[i] = -0.002 / (1.0 + 0.03 * i);

    // Build M = 1 − A·G  (G diagonal ⇒ (A·G)[i,j] = A[i,j]·G[j]).
    std::vector<double> M(DIM * DIM);
    for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j)
        M[i * DIM + j] = (i == j ? 1.0 : 0.0) - A[i * DIM + j] * G[j];

    // Solve M·X = A by Gauss-Jordan on the augmented [M | A] (DIM × 2·DIM).
    std::vector<double> aug((size_t)DIM * 2 * DIM, 0.0);
    for (int i = 0; i < DIM; ++i) {
        for (int j = 0; j < DIM; ++j) aug[i * (2 * DIM) + j] = M[i * DIM + j];
        for (int j = 0; j < DIM; ++j) aug[i * (2 * DIM) + (DIM + j)] = A[i * DIM + j];
    }
    for (int k = 0; k < DIM; ++k) {
        // partial pivot
        int pivrow = k; double maxv = std::abs(aug[k * (2 * DIM) + k]);
        for (int i = k + 1; i < DIM; ++i)
            if (std::abs(aug[i * (2 * DIM) + k]) > maxv) {
                maxv = std::abs(aug[i * (2 * DIM) + k]); pivrow = i;
            }
        if (pivrow != k)
            for (int j = 0; j < 2 * DIM; ++j)
                std::swap(aug[k * (2 * DIM) + j], aug[pivrow * (2 * DIM) + j]);
        double piv = aug[k * (2 * DIM) + k];
        if (std::abs(piv) < 1e-15) { g_failures++; std::printf("FAIL dense-solver: singular (1-AG)\n"); return; }
        for (int j = 0; j < 2 * DIM; ++j) aug[k * (2 * DIM) + j] /= piv;
        for (int i = 0; i < DIM; ++i) {
            if (i == k) continue;
            double f = aug[i * (2 * DIM) + k];
            if (f == 0.0) continue;
            for (int j = 0; j < 2 * DIM; ++j)
                aug[i * (2 * DIM) + j] -= f * aug[k * (2 * DIM) + j];
        }
    }
    // Extract X (right half of the reduced augmented matrix).
    std::vector<double> X(DIM * DIM);
    for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            X[i * DIM + j] = aug[i * (2 * DIM) + (DIM + j)];

    // Residual: (1 − A·G)·X − A == 0.  G diagonal ⇒ (A·G)[i,j] = A[i,j]·G[j].
    // (1−AG)·X = X − (AG)·X  (the "1" is the identity matrix, so 1·X = X).
    double max_res = 0.0;
    for (int i = 0; i < DIM; ++i)
    for (int col = 0; col < DIM; ++col) {
        double s = X[i * DIM + col];
        for (int j = 0; j < DIM; ++j) s -= A[i * DIM + j] * G[j] * X[j * DIM + col];
        max_res = std::max(max_res, std::abs(s - A[i * DIM + col]));
    }
    check_close("dense-solver: residual (1-AG)X - A", max_res, 0.0);
    (void)m;
}

// =============================================================================
// Test 8: Neumann series matches dense inverse to Padé tolerance.
// a_0 = A, a_{n+1} = a_n·G·A;  U_partial = Σ_{n=0}^{N} a_n  →  U_analytic.
// =============================================================================
static void test_neumann_matches_dense(const std::vector<double>& A) {
    std::vector<double> G(DIM, 0.0);
    // Keep the deliberately dense coupled-alpha toy kernel inside the Neumann
    // convergence disk. This is an algebra test, not a physical resolvent.
    for (int i = 0; i < DIM; ++i) G[i] = -0.002 / (1.0 + 0.03 * i);

    // Analytic U = (1 − A·G)⁻¹·A
    std::vector<double> M(DIM * DIM);
    for (int i = 0; i < DIM; ++i)
    for (int j = 0; j < DIM; ++j)
        M[i * DIM + j] = (i == j ? 1.0 : 0.0) - A[i * DIM + j] * G[j];
    std::vector<double> U = A;
    for (int k = 0; k < DIM; ++k) {
        int pivrow = k; double maxv = std::abs(M[k * DIM + k]);
        for (int i = k + 1; i < DIM; ++i)
            if (std::abs(M[i * DIM + k]) > maxv) { maxv = std::abs(M[i*DIM+k]); pivrow = i; }
        if (pivrow != k) {
            for (int j = 0; j < DIM; ++j) std::swap(M[k*DIM+j], M[pivrow*DIM+j]);
            for (int j = 0; j < DIM; ++j) std::swap(U[k*DIM+j], U[pivrow*DIM+j]);
        }
        double p = M[k * DIM + k]; if (std::abs(p) < 1e-15) continue;
        for (int i = k + 1; i < DIM; ++i) {
            double f = M[i * DIM + k] / p;
            for (int j = k; j < DIM; ++j) M[i * DIM + j] -= f * M[k * DIM + j];
            for (int j = 0; j < DIM; ++j) U[i * DIM + j] -= f * U[k * DIM + j];
        }
    }
    std::vector<double> U_an(DIM * DIM, 0.0);
    for (int col = 0; col < DIM; ++col)
        for (int i = DIM - 1; i >= 0; --i) {
            double s = U[i * DIM + col];
            for (int j = i + 1; j < DIM; ++j) s -= M[i * DIM + j] * U_an[j * DIM + col];
            U_an[i * DIM + col] = s / M[i * DIM + i];
        }

    // Neumann partial sum  Σ_{n=0}^{20} A·(G·A)^n
    std::vector<double> term = A;          // a_0
    std::vector<double> sum  = A;
    for (int n = 1; n <= 20; ++n) {
        // term_{n} = term_{n-1} · G · A
        std::vector<double> GA(DIM * DIM, 0.0);  // = G·A (G diagonal)
        for (int i = 0; i < DIM; ++i)
        for (int j = 0; j < DIM; ++j)
            GA[i * DIM + j] = G[i] * A[i * DIM + j];
        std::vector<double> next(DIM * DIM, 0.0);
        for (int i = 0; i < DIM; ++i)
        for (int k = 0; k < DIM; ++k) {
            double t = term[i * DIM + k]; if (t == 0.0) continue;
            for (int j = 0; j < DIM; ++j) next[i * DIM + j] += t * GA[k * DIM + j];
        }
        term = next;
        for (int i = 0; i < DIM * DIM; ++i) sum[i] += term[i];
    }
    double max_diff = 0.0;
    for (int i = 0; i < DIM * DIM; ++i)
        max_diff = std::max(max_diff, std::abs(sum[i] - U_an[i]));
    // Neumann converges to the dense inverse (looser tolerance: truncated series).
    check_close("Neumann vs dense (1-AG)^-1·A", max_diff, 0.0, 1e-4);
}

// =============================================================================
// Test 11: production Padé routine vs an independent dense direct solve for a
// two-channel non-symmetric toy kernel. Both off-diagonal alpha couplings are
// nonzero. The scalar Neumann coefficients are generated by matrix recursion,
// not from the Padé implementation.
// =============================================================================
static void test_production_pade_matches_dense_toy() {
    constexpr int D = 2;
    const cdouble A[D * D] = {
        {0.70,  0.05}, {-0.23, 0.11},
        {0.41, -0.07}, { 0.52, 0.03}
    };
    const cdouble G[D] = {{0.31, 0.04}, {-0.27, 0.02}};

    // Direct inverse of M = 1-A*G, followed by U=M^{-1}A.
    const cdouble m00 = 1.0 - A[0] * G[0];
    const cdouble m01 =     - A[1] * G[1];
    const cdouble m10 =     - A[2] * G[0];
    const cdouble m11 = 1.0 - A[3] * G[1];
    const cdouble det = m00 * m11 - m01 * m10;
    cdouble U_dense[D * D];
    U_dense[0] = ( m11 * A[0] - m01 * A[2]) / det;
    U_dense[1] = ( m11 * A[1] - m01 * A[3]) / det;
    U_dense[2] = (-m10 * A[0] + m00 * A[2]) / det;
    U_dense[3] = (-m10 * A[1] + m00 * A[3]) / det;

    // a_0=A; a_(n+1)=a_n*G*A. A 2x2 transfer matrix gives a rational
    // denominator of degree at most two, so [2/2] has enough information.
    cdouble terms[5][D * D]{};
    for (int i = 0; i < D * D; ++i) terms[0][i] = A[i];
    for (int n = 1; n < 5; ++n) {
        for (int i = 0; i < D; ++i)
        for (int j = 0; j < D; ++j)
        for (int k = 0; k < D; ++k)
            terms[n][i * D + j] += terms[n - 1][i * D + k] * G[k] * A[k * D + j];
    }

    double max_diff = 0.0;
    for (int element = 0; element < D * D; ++element) {
        cdouble coefficients[5];
        for (int n = 0; n < 5; ++n) coefficients[n] = terms[n][element];
        const cdouble U_pade = pade_approximant(coefficients, 2, 2, 1.0);
        max_diff = std::max(max_diff, std::abs(U_pade - U_dense[element]));
    }
    check_close("production Pade[2/2] vs dense two-alpha solve", max_diff, 0.0, 1e-11);
}

int main() {
    g_grid.init();
    std::mt19937_64 rng(0xC0FFEEULL);  // fixed seed → reproducible
    g_w1.fill(rng);
    DenseMock m = build_dense_mock(rng);
    std::vector<double> M_math = build_M_math(m);

    std::printf("=== Phase 0 operator-level oracle (DIM=%d) ===\n", DIM);
    test_pointer_table_layout(m);
    test_real_3S1_3D1_CT_storage();
    test_coupled_alpha_W1C(m);
    test_CPVC_col_matches_math(m, M_math);
    test_CPVC_rows_share_algebra(m, M_math);
    test_zero_W1_reproduces_2NF(m);
    test_W1_cache_matches_fallback(m);
    test_operator_ordering(m);
    test_CT_distinction(m, M_math);
    test_dense_solver_identity(m, M_math);
    test_neumann_matches_dense(M_math);
    test_production_pade_matches_dense_toy();

    std::printf("\n=== Phase 0 oracle summary ===\n");
    std::printf("  passes:    %d\n", g_passes);
    std::printf("  failures:  %d\n", g_failures);
    if (g_failures == 0) {
        std::printf("ALL PHASE-0 ORACLE CHECKS PASSED\n");
        return 0;
    }
    std::printf("PHASE-0 ORACLE FAILED\n");
    return 1;
}
