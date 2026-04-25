// ===============================================================
// 抽取自仓库 [current]: tools/check_3nf_normalization/check_3nf_normalization.cpp
// 行号区段：22..39
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
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
