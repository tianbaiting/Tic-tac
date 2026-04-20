#ifndef READ_TRITON_PSI_H
#define READ_TRITON_PSI_H

#include <vector>
#include <string>

// [EN] ψ_3H record on the file's own Gauss-Legendre grid. psi[p*Nq*Nalpha + q*Nalpha + a]
// is the amplitude of partial wave α at momentum (p_i, q_j). All momenta are in fm^-1
// as stored in tests/Cont_Faddeev/Triton_states/H3_psi_N3LO_EM500.dat.
struct triton_wavefunction {
    int Np, Nq, Nalpha;
    std::vector<double> p;   // size Np
    std::vector<double> wp;  // size Np, Gauss-Legendre weights
    std::vector<double> q;   // size Nq
    std::vector<double> wq;  // size Nq
    std::vector<int> L_2N, S_2N, J_2N, T_2N, L_1N, two_J_1N;  // size Nalpha each
    std::vector<double> psi; // size Np*Nq*Nalpha (indexing: p*Nq*Nalpha + q*Nalpha + a)
};

// Throws std::runtime_error on any parse failure.
triton_wavefunction read_triton_psi(const std::string& filename);

#endif
