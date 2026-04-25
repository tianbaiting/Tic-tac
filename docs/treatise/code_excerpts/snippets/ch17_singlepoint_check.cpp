// ===============================================================
// 抽取自仓库 [current]: tools/check_3nf_normalization/check_3nf_normalization.cpp
// 行号区段：88..156
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
                std::printf("  oracle (hand_calc_cE.py, fpi=92.2) = +7.264125e-03 fm^5\n");
                double expected = +7.264125e-03;
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
