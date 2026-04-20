#include "build_pw_3n_statespace.h"

pw_3N_statespace build_pw_3n_statespace_from_triton(const triton_wavefunction& w) {
    pw_3N_statespace s{};
    const int N = w.Nalpha;
    s.Nalpha = N;
    s.J_2N_max = 0;
    for (int a = 0; a < N; ++a) if (w.J_2N[a] > s.J_2N_max) s.J_2N_max = w.J_2N[a];
    s.L_2N_array      = new int[N];
    s.S_2N_array      = new int[N];
    s.J_2N_array      = new int[N];
    s.T_2N_array      = new int[N];
    s.L_1N_array      = new int[N];
    s.two_J_1N_array  = new int[N];
    s.two_J_3N_array  = new int[N];
    s.two_T_3N_array  = new int[N];
    s.P_3N_array      = new int[N];
    for (int a = 0; a < N; ++a) {
        s.L_2N_array[a]     = w.L_2N[a];
        s.S_2N_array[a]     = w.S_2N[a];
        s.J_2N_array[a]     = w.J_2N[a];
        s.T_2N_array[a]     = w.T_2N[a];
        s.L_1N_array[a]     = w.L_1N[a];
        s.two_J_1N_array[a] = w.two_J_1N[a];
        s.two_J_3N_array[a] = 1;   // triton is J=1/2
        s.two_T_3N_array[a] = 1;   // triton is T=1/2
        int parity_exp = w.L_2N[a] + w.L_1N[a];
        s.P_3N_array[a] = (parity_exp % 2 == 0) ? +1 : -1;
    }
    s.chn_3N_idx_array = nullptr;
    s.N_chn_3N = 1;
    return s;
}

void free_pw_3n_statespace_triton(pw_3N_statespace& s) {
    delete[] s.L_2N_array;      s.L_2N_array = nullptr;
    delete[] s.S_2N_array;      s.S_2N_array = nullptr;
    delete[] s.J_2N_array;      s.J_2N_array = nullptr;
    delete[] s.T_2N_array;      s.T_2N_array = nullptr;
    delete[] s.L_1N_array;      s.L_1N_array = nullptr;
    delete[] s.two_J_1N_array;  s.two_J_1N_array = nullptr;
    delete[] s.two_J_3N_array;  s.two_J_3N_array = nullptr;
    delete[] s.two_T_3N_array;  s.two_T_3N_array = nullptr;
    delete[] s.P_3N_array;      s.P_3N_array = nullptr;
    s.Nalpha = 0;
}
