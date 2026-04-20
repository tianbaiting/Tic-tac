#ifndef CONTRACT_W1_EXPECTATION_H
#define CONTRACT_W1_EXPECTATION_H

#include "read_triton_psi.h"
#include "type_defs.h"

class three_nucleon_force_model;

// [EN] Compute ⟨ψ|W^(1)|ψ⟩ as a 6D sum on ψ's own quadrature grid.
// Radial convention: weight = Π (wp·wq·p²·q²) over (row, col); ψ is the full
// antisymmetric Ψ loaded from H3_psiasymm_*.dat. W1_element returns fm^5.
// Caller multiplies by inv_hbarc5 to convert to MeV.
double contract_W1_expectation(const triton_wavefunction& w,
                               const pw_3N_statespace& pw,
                               const three_nucleon_force_model& tnf);

#endif
