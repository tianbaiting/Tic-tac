#ifndef TICTAC_CORE_THREE_BODY_CHANNEL_H
#define TICTAC_CORE_THREE_BODY_CHANNEL_H

#include <cstddef>

namespace tictac::core {

// [EN] Mathematical description of one three-body partial wave, i.e. one value
// of the composite channel label
//
//     alpha = ( L_pair, S_pair, J_pair, T_pair,        // antisymmetric (2,3) pair
//               lambda,     two_j_spectator,            // spectator (particle 1)
//               two_J, two_T, parity );                 // conserved 3N quantum numbers
//
// (see docs/three_nf_equation_contract.md §4.1). This is a plain value type so
// it can be returned by value from a basis view without owning any storage.
// / [CN] 一个三体分波通道 alpha 的数学描述（纯值类型，不持有存储）。
struct ThreeBodyChannel {
	int L_pair;            // pair relative orbital angular momentum      (L_2N)
	int S_pair;            // pair total spin                              (S_2N)
	int J_pair;            // pair total angular momentum                  (J_2N)
	int T_pair;            // pair total isospin                           (T_2N)
	int lambda;            // spectator orbital angular momentum           (L_1N)
	int two_j_spectator;   // spectator total angular momentum x2          (2j_1N)
	int two_J;             // three-body total angular momentum x2         (2J_3N)  [conserved]
	int two_T;             // three-body total isospin x2                  (2T_3N)  [conserved]
	int parity;            // three-body parity  (+1 / -1)                 (P_3N)   [conserved]
};

// Pair-antisymmetry selection rule: (-1)^(L_pair+S_pair+T_pair) = -1.
// (docs/n2lo_3nf_conventions.md: "Allowed channels obey (-1)^(l+s+t)=-1".)
inline bool pair_is_antisymmetric(const ThreeBodyChannel& ch)
{
	const int sign = ((ch.L_pair + ch.S_pair + ch.T_pair) % 2 + 2) % 2;
	return sign == 1;  // odd -> antisymmetric
}

} // namespace tictac::core

#endif // TICTAC_CORE_THREE_BODY_CHANNEL_H
