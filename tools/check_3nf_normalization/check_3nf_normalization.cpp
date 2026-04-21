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
#include <cmath>
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
        std::fprintf(stderr, "usage: %s <H3_psiasymm_*.dat>\n"
                             "  Expects the antisymmetric triton wavefunction Ψ (not the Faddeev ψ).\n"
                             "  Radial-convention norm ⟨Ψ|Ψ⟩ must equal 1.0.\n", argv[0]);
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
        std::printf("⟨Ψ|Ψ⟩ (radial p²q² weights): %.6f\n", norm_radial);
        if (std::fabs(norm_radial - 1.0) > 0.01) {
            std::fprintf(stderr,
                "WARNING: norm=%.6f is not 1.0 ± 0.01. Did you pass H3_psiasymm_*.dat?\n"
                "  H3_psi_*.dat is the Faddeev component (norm≈0.154, NOT usable here).\n",
                norm_radial);
        }
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
        // --- Task-3 single-point verification against Python oracle ---------
        // Diagonal 3S1 channel, p=q=p'=q'=0.5 fm^-1, c_E-only (c_D=c_1=c_3=0).
        // Expected (hand_calc_cE.py, fpi=92.2 MeV, with 1/(2π)³ Fourier
        // normalization mirroring chiral_LO_internal.cpp:59): W1 = -7.264e-03 fm^5.
        {
            int a_3S1 = -1;
            for (int a = 0; a < pw.Nalpha; ++a) {
                if (pw.L_2N_array[a] == 0 && pw.S_2N_array[a] == 1 &&
                    pw.J_2N_array[a] == 1 && pw.T_2N_array[a] == 0 &&
                    pw.L_1N_array[a] == 0 && pw.two_J_1N_array[a] == 1) {
                    a_3S1 = a; break;
                }
            }
            if (a_3S1 >= 0) {
                chiral_N2LO_3NF tnf_cE(/*c_D=*/0.0, /*c_E=*/-0.205, /*Lambda_MeV=*/500.0,
                                        /*c1=*/0.0, /*c3=*/0.0, /*c4=*/0.0);
                double v_fm5 = tnf_cE.W1_element(a_3S1, a_3S1,
                                                 0.5, 0.5, 0.5, 0.5, pw);
                std::printf("\n[Task-3 single-point check]\n");
                std::printf("  W1_cE(3S1, p=q=p'=q'=0.5 fm^-1) = %+.6e fm^5\n", v_fm5);
                std::printf("  oracle (hand_calc_cE.py, fpi=92.2) = -7.264125e-03 fm^5\n");
                double expected = -7.264125e-03;
                double rel = (v_fm5 - expected) / expected * 100.0;
                std::printf("  relative difference = %+.3f%%  (pass within 2%%)\n", rel);
            } else {
                std::printf("\n[Task-3 single-point check] WARNING: 3S1 channel not found.\n");
            }
        }

        // --- Task-4 single-point verification against Python oracle ---------
        // Off-diagonal channel (1S0 bra ↔ 3S1 ket), p=q=p'=q'=0.5 fm^-1.
        // bra: (L_2N=0,S_2N=0,J_2N=0,T_2N=1,L_1N=0,2J_1N=1) — 1S0 + s1=1/2
        // ket: (L_2N=0,S_2N=1,J_2N=1,T_2N=0,L_1N=0,2J_1N=1) — 3S1 + s1=1/2
        // Expected (C++ oracle with corrected phase, python verified):
        //   W1_cD(1S0->3S1, p=q=p'=q'=0.5) = +3.392900e-03 fm^5
        {
            int a_1S0 = -1, a_3S1 = -1;
            for (int a = 0; a < pw.Nalpha; ++a) {
                if (pw.L_2N_array[a] == 0 && pw.S_2N_array[a] == 0 &&
                    pw.J_2N_array[a] == 0 && pw.T_2N_array[a] == 1 &&
                    pw.L_1N_array[a] == 0 && pw.two_J_1N_array[a] == 1 &&
                    pw.two_J_3N_array[a] == 1 && pw.two_T_3N_array[a] == 1) {
                    a_1S0 = a;
                }
                if (pw.L_2N_array[a] == 0 && pw.S_2N_array[a] == 1 &&
                    pw.J_2N_array[a] == 1 && pw.T_2N_array[a] == 0 &&
                    pw.L_1N_array[a] == 0 && pw.two_J_1N_array[a] == 1 &&
                    pw.two_J_3N_array[a] == 1 && pw.two_T_3N_array[a] == 1) {
                    a_3S1 = a;
                }
            }
            std::printf("\n[Task-4 single-point check]\n");
            if (a_1S0 >= 0 && a_3S1 >= 0) {
                chiral_N2LO_3NF tnf_cD(/*c_D=*/-0.20, /*c_E=*/0.0, /*Lambda_MeV=*/500.0,
                                        /*c1=*/0.0, /*c3=*/0.0, /*c4=*/0.0);
                // W1_element(bra=1S0, ket=3S1) — off-diagonal c_D element
                double v_fm5 = tnf_cD.W1_element(a_1S0, a_3S1,
                                                  0.5, 0.5, 0.5, 0.5, pw);
                std::printf("  W1_cD(1S0->3S1, p=q=p'=q'=0.5 fm^-1)  = %+.6e fm^5\n", v_fm5);
                std::printf("  oracle (C++ phase + 24-pt GL)           = +3.392900e-03 fm^5\n");
                double expected_cD = +3.392900e-03;
                double rel_cD = (v_fm5 - expected_cD) / expected_cD * 100.0;
                std::printf("  relative difference = %+.3f%%  (pass within 2%%)\n", rel_cD);
                if (std::fabs(rel_cD) > 2.0) {
                    std::printf("  WARNING: relative difference > 2%% — W1_1pe_contact may have a bug.\n");
                }
            } else {
                std::printf("  WARNING: 1S0 channel (a_1S0=%d) or 3S1 channel (a_3S1=%d) not found.\n",
                            a_1S0, a_3S1);
                std::printf("  The triton wavefunction may use a different 2J_3N sector.\n");
            }
        }

        std::printf("\n=== ΔE_3NF = ⟨Ψ|(1+P+P²)W^(1)|Ψ⟩ = 3·⟨Ψ|W^(1)|Ψ⟩ (antisymmetric Ψ → P|Ψ⟩ = |Ψ⟩) ===\n");
        std::printf("%-16s  %18s  %18s\n", "channel", "⟨W⟩ [MeV]", "3·⟨W⟩ [MeV]");
        double v_MeV_values[5];
        for (size_t i = 0; i < sizeof(runs)/sizeof(runs[0]); ++i) {
            const auto& r = runs[i];
            chiral_N2LO_3NF tnf(r.c_D, r.c_E, /*Lambda_MeV=*/500.0, r.c_1, r.c_3, /*c4=*/0.0);
            double v_fminv = contract_W1_expectation(w, pw, tnf);
            double v_MeV = v_fminv * hbarc;
            v_MeV_values[i] = v_MeV;
            std::printf("%-16s  %+18.6e  %+18.6e\n", r.name, v_MeV, 3.0*v_MeV);
        }
        // Additivity check on the symmetrized ΔE_3NF
        double sum_individual = 3.0 * (v_MeV_values[0] + v_MeV_values[1] + v_MeV_values[2] + v_MeV_values[3]);
        double full_ΔE = 3.0 * v_MeV_values[4];
        std::printf("Additivity check (3·⟨W⟩): sum(c_E+c_D+c_1+c_3) = %+e MeV; full = %+e MeV; residual = %+e MeV\n",
                    sum_individual, full_ΔE, full_ΔE - sum_individual);

        // ---- Comparison vs Epelbaum PRC 66 (2002) reference (see epelbaum_reference.md) ----
        // Reference values are Epelbaum Table 2, ³H column, NNLO Λ=500 MeV, linearly
        // rescaled from Epelbaum's LEC convention (c_D=+3.6, c_E=+0.37) to ours
        // (c_D=-0.20, c_E=-0.205).  c_1/c_3 split is an ESTIMATE (Epelbaum reports
        // only combined c-terms = -0.39 MeV).  Robust combined c_1+c_3 ref = -0.367 MeV.
        struct ref_row { const char* name; double ref_MeV; };
        const ref_row refs[] = {
            { "c_E_only",    +0.410 },   // transcribed + rescaled (-0.74 × -0.205/0.37)
            { "c_D_only",    -0.045 },   // transcribed + rescaled (+0.81 × -0.20/3.6)
            { "c_1_only",    -0.073 },   // ESTIMATED — see epelbaum_reference.md (c_1/c_3 split ~20%)
            { "c_3_only",    -0.294 },   // ESTIMATED — see epelbaum_reference.md (c_1/c_3 split ~80%)
            { "full_Witala", -0.002 },   // sum of rescaled individual contributions at our LECs
        };
        std::printf("\n=== Comparison vs Epelbaum/Witała reference ===\n");
        std::printf("%-16s  %18s  %18s  %10s\n",
                    "channel", "3·<W>_code [MeV]", "<W>_ref [MeV]", "ratio X");
        for (size_t i = 0; i < sizeof(runs)/sizeof(runs[0]); ++i) {
            double code_MeV = 3.0 * v_MeV_values[i];
            double ref_MeV = refs[i].ref_MeV;
            double ratio = (ref_MeV != 0.0) ? (code_MeV / ref_MeV) : 0.0;
            std::printf("%-16s  %+18.6e  %+18.6e  %+10.3f\n",
                        refs[i].name, code_MeV, ref_MeV, ratio);
        }
        // Robust combined c_1+c_3 ratio (does not depend on the c_1/c_3 split estimate).
        {
            double code_combined = 3.0 * (v_MeV_values[2] + v_MeV_values[3]);
            double ref_combined  = -0.367; // Epelbaum Table 2 × (our c_3 / Epelbaum c_3)
            double ratio_comb    = code_combined / ref_combined;
            std::printf("%-16s  %+18.6e  %+18.6e  %+10.3f   (robust: c-terms combined)\n",
                        "c_1+c_3 comb",
                        code_combined, ref_combined, ratio_comb);
        }
        std::printf("\nInterpretation: if |ratio X_i| is the same (within ±5%%) across all four single-LEC rows,\n"
                    "a single multiplicative correction 1/X can be applied in chiral_N2LO_3NF.h.\n"
                    "Different |X_i| values indicate per-term bugs.\n");

        free_pw_3n_statespace_triton(pw);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
}
