#include "read_triton_psi.h"
#include <cstdio>
#include <exception>

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
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 2;
    }
}
