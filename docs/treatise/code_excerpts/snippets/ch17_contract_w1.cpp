// ===============================================================
// 抽取自仓库 [current]: tools/check_3nf_normalization/contract_W1_expectation.cpp
// 行号区段：1..43
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
#include "contract_W1_expectation.h"
#include "three_nucleon_force_model.h"
#include <omp.h>

// Radial convention (pinned in Task 3):
// ⟨Ψ|W|Ψ⟩ = Σ wp'·wq'·p'²·q'² · wp·wq·p²·q² · Ψ(α',p',q') · W^(1)(α',α, p',q', p,q) · Ψ(α,p,q)
double contract_W1_expectation(const triton_wavefunction& w,
                               const pw_3N_statespace& pw,
                               const three_nucleon_force_model& tnf) {
    const int Np = w.Np, Nq = w.Nq, Na = w.Nalpha;
    auto psi_at = [&](int i, int j, int a) -> double {
        return w.psi[(size_t)i * Nq * Na + (size_t)j * Na + a];
    };
    double total = 0.0;
    #pragma omp parallel for reduction(+:total) schedule(dynamic) collapse(2)
    for (int ip = 0; ip < Np; ++ip) {
        for (int iq = 0; iq < Nq; ++iq) {
            double wp_r = w.wp[ip], wq_r = w.wq[iq];
            double p_r = w.p[ip], q_r = w.q[iq];
            double row_factor = wp_r * wq_r * p_r * p_r * q_r * q_r;
            for (int jp = 0; jp < Np; ++jp) {
                double wp_c = w.wp[jp], p_c = w.p[jp];
                for (int jq = 0; jq < Nq; ++jq) {
                    double wq_c = w.wq[jq], q_c = w.q[jq];
                    double col_factor = wp_c * wq_c * p_c * p_c * q_c * q_c;
                    double weight = row_factor * col_factor;
                    for (int ar = 0; ar < Na; ++ar) {
                        double psi_r = psi_at(ip, iq, ar);
                        if (psi_r == 0.0) continue;
                        for (int ac = 0; ac < Na; ++ac) {
                            double psi_c = psi_at(jp, jq, ac);
                            if (psi_c == 0.0) continue;
                            double w1 = tnf.W1_element(ar, ac, p_r, q_r, p_c, q_c, pw);
                            if (w1 == 0.0) continue;
                            total += weight * psi_r * w1 * psi_c;
                        }
                    }
                }
            }
        }
    }
    return total;
}
