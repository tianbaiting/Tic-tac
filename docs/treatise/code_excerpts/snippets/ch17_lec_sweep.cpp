// ===============================================================
// 抽取自仓库 [current]: tools/check_3nf_normalization/check_3nf_normalization.cpp
// 行号区段：206..260
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
