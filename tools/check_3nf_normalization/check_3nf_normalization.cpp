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
#include "build_pw_3n_statespace.h"
#include "contract_W1_expectation.h"
#include "three_nucleon_force_model.h"
#include "chiral_N2LO_3NF.h"
#include "constants.h"
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
        pw_3N_statespace pw = build_pw_3n_statespace_from_triton(w);
        std::printf("Built pw_3N_statespace: Nalpha=%d J_2N_max=%d\n", pw.Nalpha, pw.J_2N_max);
        int p0 = pw.P_3N_array[0];
        bool all_same_parity = true;
        for (int a = 0; a < pw.Nalpha; ++a) {
            if (pw.P_3N_array[a] != p0) { all_same_parity = false; break; }
        }
        std::printf("All α share parity P_3N=%+d: %s\n", p0, all_same_parity ? "yes" : "no");
        struct lec_run { const char* name; double c_D, c_E, c_1, c_3; };
        // c_1, c_3 per EM500 (Entem–Machleidt) 2NF: c_1=-0.81, c_3=-3.20 GeV^-1
        // (Witała PRC 77 (2008) 034004)
        const lec_run runs[] = {
            { "c_E_only",    0.0,   -0.205,  0.0,   0.0   },
            { "c_D_only",   -0.20,   0.0,    0.0,   0.0   },
            { "c_1_only",    0.0,    0.0,   -0.81,  0.0   },
            { "c_3_only",    0.0,    0.0,    0.0,  -3.20  },
            { "full_Witala",-0.20,  -0.205, -0.81, -3.20  },
        };
        std::printf("\n=== W^(1) first-order expectation (no (1+P+P²) yet) ===\n");
        std::printf("%-16s  %18s  %18s\n", "channel", "⟨W⟩ [fm^-1]", "⟨W⟩ [MeV]");
        double v_MeV_values[5];
        for (size_t i = 0; i < sizeof(runs)/sizeof(runs[0]); ++i) {
            const auto& r = runs[i];
            chiral_N2LO_3NF tnf(r.c_D, r.c_E, /*Lambda_MeV=*/500.0, r.c_1, r.c_3, /*c4=*/0.0);
            double v_fminv = contract_W1_expectation(w, pw, tnf);
            double v_MeV = v_fminv * hbarc;
            v_MeV_values[i] = v_MeV;
            std::printf("%-16s  %+18.6e  %+18.6e\n", r.name, v_fminv, v_MeV);
        }
        // Additivity check: full_Witala == c_E + c_D + c_1 + c_3
        double sum_individual = v_MeV_values[0] + v_MeV_values[1] + v_MeV_values[2] + v_MeV_values[3];
        double additivity_residual = v_MeV_values[4] - sum_individual;
        std::printf("Additivity check: sum(c_E+c_D+c_1+c_3) = %+e MeV; full = %+e MeV; residual = %+e MeV\n",
                    sum_individual, v_MeV_values[4], additivity_residual);
        free_pw_3n_statespace_triton(pw);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
}
