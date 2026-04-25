// ===============================================================
// 抽取自仓库 [current]: tools/check_3nf_normalization/check_3nf_normalization.cpp
// 行号区段：41..86
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
