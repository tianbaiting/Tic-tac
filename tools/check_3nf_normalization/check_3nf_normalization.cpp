// [EN] Integration convention pinned by Task 3 (2026-04-20):
// The file H3_psi_N3LO_EM500.dat stores the FADDEEV COMPONENT ψ (not the
// fully antisymmetric state Ψ = (1+P+P²)ψ).  The Faddeev component is not
// individually normalized to 1; only the antisymmetric combination is.
// Cross-checking with H3_psiasymm_N3LO_EM500.dat (the fully antisymmetric
// state) confirmed:
//   ⟨Ψ|Ψ⟩ = Σ wp·wq · p²·q² · |Ψ|²  = 1.000000  (radial convention)
//   ⟨Ψ|Ψ⟩ = Σ wp·wq · |Ψ|²           = 194.768   (reduced convention)
// Therefore the RADIAL convention (integration measure dp p² dq q²) is correct.
// For the Faddeev component ψ: norm_radial = 0.153637 ≈ 1/6.5 (not 1),
// consistent with ψ = Ψ/(1+P+P²) and the antisymmetrization factor.
// All subsequent tasks use: ⟨f|g⟩ = Σ_α Σ_{ij} wp_i·wq_j·p_i²·q_j² f_α(p_i,q_j) g_α(p_i,q_j)

#include "read_triton_psi.h"
#include <cstdio>
#include <exception>

static double compute_norm_radial(const triton_wavefunction& w) {
    // Radial-function convention: ⟨ψ|ψ⟩ = Σ wp·wq · p²·q² · |ψ|²
    double s = 0.0;
    for (int i = 0; i < w.Np; ++i) {
        double p2 = w.p[i] * w.p[i];
        for (int j = 0; j < w.Nq; ++j) {
            double q2 = w.q[j] * w.q[j];
            double factor = w.wp[i] * w.wq[j] * p2 * q2;
            for (int a = 0; a < w.Nalpha; ++a) {
                double amp = w.psi[(size_t)i * w.Nq * w.Nalpha + (size_t)j * w.Nalpha + a];
                s += factor * amp * amp;
            }
        }
    }
    return s;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <psi_3H_file>\n", argv[0]);
        return 1;
    }
    try {
        triton_wavefunction w = read_triton_psi(argv[1]);
        std::printf("Np=%d Nq=%d Nalpha=%d\n", w.Np, w.Nq, w.Nalpha);
        std::printf("p range: [%.4f, %.4f] fm^-1 (Np=%d)\n", w.p.front(), w.p.back(), w.Np);
        std::printf("q range: [%.4f, %.4f] fm^-1 (Nq=%d)\n", w.q.front(), w.q.back(), w.Nq);
        std::printf("First 3 alpha channels (L_2N, S_2N, J_2N, T_2N, L_1N, 2J_1N):\n");
        for (int a = 0; a < 3 && a < w.Nalpha; ++a) {
            std::printf("  a=%d: %d %d %d %d %d %d\n", a,
                        w.L_2N[a], w.S_2N[a], w.J_2N[a], w.T_2N[a],
                        w.L_1N[a], w.two_J_1N[a]);
        }
        double norm_radial = compute_norm_radial(w);
        std::printf("⟨ψ|ψ⟩ (radial p²q² weights, Faddeev component): %.6f\n", norm_radial);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
}
