#include "chiral_3nf_recoupling.h"
#include "coupling_coefficients.h"
#include <cmath>

// [EN] Implementation notes:
//   We call the GSL-backed `wigner_*` wrappers (from coupling_coefficients.h)
//   rather than the cached `SixJSymbol` helper — the recoupling helpers run
//   only once per (alpha_r, alpha_c) pair outside the radial hot loop, so
//   the direct GSL path is preferred for simplicity and correctness.
//
// / [CN] 说明：直接调用 GSL 后端而非缓存版本，因为这些重耦合系数只在
// (alpha_r, alpha_c) 通道对外层循环中计算一次，远离径向热循环。

// Helper: (2*J + 1) hat-squared factor.
static inline double two_jp1(int J) { return 2.0 * J + 1.0; }

double recoupling_3nf_scalar(
    int L_2N_r, int S_2N_r, int J_2N_r, int T_2N_r,
    int L_1N_r, int two_J_1N_r, int two_J_3N,
    int L_2N_c, int S_2N_c, int J_2N_c, int T_2N_c,
    int L_1N_c, int two_J_1N_c,
    int two_T_3N)
{
    // [EN] Pair-scalar operator: sigma_2.sigma_3 * tau_2.tau_3 acting on
    // the pair (2,3) with all pair orbital/spin/isospin quantum numbers
    // conserved and the spectator line completely untouched. The full
    // 3N matrix element reduces to a Kronecker product of selection
    // rules times the pair eigenvalues.
    //
    // / [CN] 配对标量算子：sigma_2.sigma_3 * tau_2.tau_3 作用于配对 (2,3)，
    // 保持所有配对轨道/自旋/同位旋量子数，旁观者线完全不变。

    // Selection rules — pair orbital & spin & isospin all conserved,
    // spectator completely unchanged, J_3N and T_3N already assumed equal
    // by the caller but double-check via delta factors below (2T_3N, 2J_3N
    // are shared arguments so nothing to compare against — they only
    // multiply the operator eigenvalues through conservation).
    if (L_2N_r != L_2N_c) return 0.0;
    if (S_2N_r != S_2N_c) return 0.0;
    if (J_2N_r != J_2N_c) return 0.0;
    if (T_2N_r != T_2N_c) return 0.0;
    if (L_1N_r != L_1N_c) return 0.0;
    if (two_J_1N_r != two_J_1N_c) return 0.0;

    // Suppress unused warnings for J_3N/T_3N arguments — they are part of
    // the interface for future generalisation but the current rank-0
    // coefficient is J_3N / T_3N independent (pair eigenvalues only).
    (void)two_J_3N;
    (void)two_T_3N;

    // Pair spin eigenvalue: sigma_2 . sigma_3 = 2 S_pair(S_pair+1) - 3
    //   S_2N = 0 -> -3   (singlet)
    //   S_2N = 1 -> +1   (triplet)
    const double sigma_sigma = 2.0 * S_2N_r * (S_2N_r + 1) - 3.0;

    // Pair isospin eigenvalue: tau_2 . tau_3 = 2 T_pair(T_pair+1) - 3
    //   T_2N = 0 -> -3   (isoscalar pair)
    //   T_2N = 1 -> +1   (isovector pair)
    const double tau_tau = 2.0 * T_2N_r * (T_2N_r + 1) - 3.0;

    return sigma_sigma * tau_tau;
}

double recoupling_3nf_rank2(
    int L_2N_r, int S_2N_r, int J_2N_r, int T_2N_r,
    int L_1N_r, int two_J_1N_r, int two_J_3N,
    int L_2N_c, int S_2N_c, int J_2N_c, int T_2N_c,
    int L_1N_c, int two_J_1N_c,
    int two_T_3N)
{
    // [EN] Rank-2 pair operator [sigma_2 x sigma_3]_2 . [spatial]_2 with
    // an outer tau_2.tau_3 isospin scalar. Formula from
    // formula_reference.md §3.4 / Hebeler PRC 91 (2015) eq. 15.
    //
    // / [CN] 配对 rank-2 算子：[sigma_2 x sigma_3]_2 . [空间]_2 乘以
    // tau_2.tau_3 同位旋标量。

    // --- Selection rules ---
    // Rank-2 spatial Y_2: requires L_2N' in { |L_2N-2|, ..., L_2N+2 }
    // and CG(L_2N 0; 2 0 | L_2N' 0) nonzero -> |L_2N - L_2N'| in {0, 2}
    // with L_2N + L_2N' + 2 even.
    const int dL = std::abs(L_2N_r - L_2N_c);
    if (dL != 0 && dL != 2) return 0.0;
    if ((L_2N_r + L_2N_c + 2) % 2 != 0) return 0.0;

    // Rank-2 pair spin: S=S'=1 only (singlet has no rank-2 on two spins).
    if (S_2N_r != 1 || S_2N_c != 1) return 0.0;

    // Pair isospin is scalar tau_2.tau_3 -> diagonal in T_2N.
    if (T_2N_r != T_2N_c) return 0.0;

    // Spectator line untouched by the pair-only rank-2 operator.
    if (L_1N_r != L_1N_c) return 0.0;
    if (two_J_1N_r != two_J_1N_c) return 0.0;

    // Pair-level rank-2 coupled to scalar: requires J_2N = J_2N'.
    // (The 9j bottom-right 0 enforces this; we short-circuit for speed.)
    if (J_2N_r != J_2N_c) return 0.0;

    // Suppress unused warnings for J_3N — conservation is already assumed
    // by the caller; the rank-2 coefficient is independent of J_3N because
    // the operator is scalar at the 3N level (rank-0 in J_3N space) and
    // the pair/spectator coupling is diagonal.
    (void)two_J_3N;

    // --- Pair isospin eigenvalue: tau_2 . tau_3 ---
    const double tau_tau = 2.0 * T_2N_r * (T_2N_r + 1) - 3.0;
    if (tau_tau == 0.0) return 0.0;

    // --- Also check via 3N-isospin recoupling is trivially satisfied:
    // (tau_2.tau_3) is diagonal in (T_2N, two_T_3N) so the operator value
    // is just the pair eigenvalue times the caller's conservation delta.
    (void)two_T_3N;

    // --- Reduced matrix element < 1 || [sigma_2 x sigma_3]_2 || 1 > = sqrt(30) ---
    // (Standard Racah result; see Edmonds 7.1.8, or formula_reference.md §5.)
    const double rme_sigma_tensor = std::sqrt(30.0);

    // --- Spatial rank-2 reduced matrix element via Clebsch-Gordan ---
    // <L' 0 | 2 0 ; L 0> using twice-integer convention.
    const double cg_spatial = clebsch_gordan(2 * L_2N_c, 4, 2 * L_2N_r,
                                             0, 0, 0);
    if (cg_spatial == 0.0) return 0.0;

    // --- 9j recoupling of (L_r S_r J_r) <-> (L_c S_c J_c) via rank-2x rank-2 coupled to rank-0 ---
    // NineJ(L_r, S_r, J_r; L_c, S_c, J_c; 2, 2, 0)
    // Bottom row (2 2 0) enforces J_r = J_c (scalar at pair level).
    const double w9j = wigner_9j(2 * L_2N_r, 2 * S_2N_r, 2 * J_2N_r,
                                 2 * L_2N_c, 2 * S_2N_c, 2 * J_2N_c,
                                 4,           4,           0);
    if (w9j == 0.0) return 0.0;

    // --- Hat factors from the 9j coupling normalisation ---
    // sqrt((2S+1)(2S'+1)(2J+1)(2J'+1)) per Hebeler eq. 15 style.
    const double hat = std::sqrt(
        two_jp1(S_2N_r) * two_jp1(S_2N_c)
      * two_jp1(J_2N_r) * two_jp1(J_2N_c));

    return tau_tau * rme_sigma_tensor * hat * w9j * cg_spatial;
}
