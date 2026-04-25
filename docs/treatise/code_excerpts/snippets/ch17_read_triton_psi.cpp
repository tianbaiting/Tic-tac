// ===============================================================
// 抽取自仓库 [current]: tools/check_3nf_normalization/read_triton_psi.cpp
// 行号区段：1..48
// 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
// ===============================================================
#include "read_triton_psi.h"
#include <cstdio>
#include <stdexcept>
#include <string>

static void expect_read(int got, int want, const char* ctx) {
    if (got != want) {
        throw std::runtime_error(std::string("read_triton_psi: parse failure at ") + ctx);
    }
}

triton_wavefunction read_triton_psi(const std::string& filename) {
    triton_wavefunction out;
    FILE* f = std::fopen(filename.c_str(), "r");
    if (!f) throw std::runtime_error("read_triton_psi: cannot open " + filename);

    expect_read(std::fscanf(f, "%d %d %d", &out.Np, &out.Nq, &out.Nalpha), 3, "dimensions");

    out.p.resize(out.Np);
    out.wp.resize(out.Np);
    for (int i = 0; i < out.Np; ++i) {
        expect_read(std::fscanf(f, "%lf %lf", &out.p[i], &out.wp[i]), 2, "p row");
    }
    out.q.resize(out.Nq);
    out.wq.resize(out.Nq);
    for (int i = 0; i < out.Nq; ++i) {
        expect_read(std::fscanf(f, "%lf %lf", &out.q[i], &out.wq[i]), 2, "q row");
    }
    out.L_2N.resize(out.Nalpha);
    out.S_2N.resize(out.Nalpha);
    out.J_2N.resize(out.Nalpha);
    out.T_2N.resize(out.Nalpha);
    out.L_1N.resize(out.Nalpha);
    out.two_J_1N.resize(out.Nalpha);
    for (int a = 0; a < out.Nalpha; ++a) {
        expect_read(std::fscanf(f, "%d %d %d %d %d %d",
                                &out.L_2N[a], &out.S_2N[a], &out.J_2N[a],
                                &out.T_2N[a], &out.L_1N[a], &out.two_J_1N[a]),
                    6, "alpha row");
    }
    const size_t N = (size_t)out.Np * out.Nq * out.Nalpha;
    out.psi.assign(N, 0.0);
    for (size_t k = 0; k < N; ++k) {
        expect_read(std::fscanf(f, "%lf", &out.psi[k]), 1, "psi entry");
    }
    std::fclose(f);
    return out;
}
