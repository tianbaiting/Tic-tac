// Semantic view unit tests (Phase 3 of the structural refactor).
//
// These views are non-owning aliases over the legacy SoA storage
// (pw_3N_statespace / fwp_statespace / the P123 CSC arrays). The tests build
// small synthetic storage, wrap it, and check the accessors return exactly the
// underlying values -- i.e. the views add ZERO numerical impact.
//
// Returns 0 on success.

#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

#include "three_body_channel.h"
#include "partial_wave_basis_view.h"
#include "packet_grid_view.h"
#include "permutation_operator_view.h"

#include "type_defs.h"

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { ++failures; \
	std::printf("FAIL: %s (line %d)\n", #cond, __LINE__); } } while (0)

static void test_three_body_channel_and_pair_antisymmetry()
{
	using namespace tictac::core;
	// Physical 3S1 deuteron-like pair: L=0,S=1,J=1,T=0 -> (-1)^(0+1+0) = -1 ok.
	ThreeBodyChannel deuteron_like{0, 1, 1, 0, 0, 1, 1, 0, +1};
	CHECK(pair_is_antisymmetric(deuteron_like));
	// Forbidden L=0,S=0,T=0 -> (-1)^0 = +1 (symmetric, not a physical NN state).
	ThreeBodyChannel forbidden{0, 0, 0, 0, 0, 1, 1, 0, +1};
	CHECK(!pair_is_antisymmetric(forbidden));
}

static void test_partial_wave_basis_view()
{
	using namespace tictac::core;
	// Two conserved blocks: block 0 has 2 alphas, block 1 has 1 alpha.
	pw_3N_statespace pw{};
	pw.Nalpha = 3;
	pw.J_2N_max = 1;
	int L2[]  = {0, 2, 1};      // block0: 3S1/3D1 pair; block1: 1P1
	int S2[]  = {1, 1, 0};
	int J2[]  = {1, 1, 1};
	int T2[]  = {0, 0, 1};
	int L1[]  = {0, 0, 1};      // spectator lambda
	int j1[]  = {1, 1, 3};      // 2j spectator
	int J3[]  = {1, 1, 3};      // block0 two_J=1 (J=1/2), block1 two_J=3 (J=3/2)
	int T3[]  = {0, 0, 0};
	int P3[]  = {+1, +1, -1};   // block0 positive, block1 negative
	int chn[] = {0, 2, 3};
	pw.L_2N_array = L2; pw.S_2N_array = S2; pw.J_2N_array = J2; pw.T_2N_array = T2;
	pw.L_1N_array = L1; pw.two_J_1N_array = j1;
	pw.two_J_3N_array = J3; pw.two_T_3N_array = T3; pw.P_3N_array = P3;
	pw.chn_3N_idx_array = chn;
	pw.N_chn_3N = 2;

	PartialWaveBasisView basis(pw);
	CHECK(basis.num_channels() == 3);
	CHECK(basis.num_blocks() == 2);
	CHECK(basis.J_2N_max() == 1);

	ThreeBodyChannel a0 = basis.channel(0);
	CHECK(a0.L_pair == 0 && a0.S_pair == 1 && a0.J_pair == 1 && a0.T_pair == 0);
	CHECK(a0.lambda == 0 && a0.two_j_spectator == 1);
	CHECK(a0.two_J == 1 && a0.two_T == 0 && a0.parity == +1);
	CHECK(pair_is_antisymmetric(a0));

	auto r0 = basis.block_range(0);
	CHECK(r0.begin == 0 && r0.end == 2);
	CHECK(basis.block_size(0) == 2);
	CHECK(basis.block_two_J(0) == 1);
	CHECK(basis.block_parity(0) == +1);

	auto r1 = basis.block_range(1);
	CHECK(r1.begin == 2 && r1.end == 3);
	CHECK(basis.block_size(1) == 1);
	CHECK(basis.block_two_J(1) == 3);
	CHECK(basis.block_parity(1) == -1);
}

static void test_packet_grid_view()
{
	using namespace tictac::core;
	// 2 p-packets x 1 node/packet; 3 q-packets x 2 nodes/packet.
	fwp_statespace fwp{};
	fwp.Np_WP = 2; fwp.Np_per_WP = 1;
	fwp.Nq_WP = 3; fwp.Nq_per_WP = 2;
	double pB[] = {0.0, 1.0, 2.0};
	double pA[] = {0.5, 1.5};
	double wp[] = {1.0, 1.0};
	double fp[] = {1.0, 1.0};
	double np[] = {1.0, 1.0};
	double qB[] = {0.0, 1.0, 2.0, 3.0};
	double qA[] = {0.5, 0.5,  1.5, 1.5,  2.5, 2.5};
	double wq[] = {0.5, 0.5,  0.5, 0.5,  0.5, 0.5};
	double fq[] = {1.0, 1.0,  1.0, 1.0,  1.0, 1.0};
	double nq[] = {1.0, 1.0,  1.0, 1.0,  1.0, 1.0};
	fwp.p_WP_array = pB; fwp.p_array = pA; fwp.wp_array = wp; fwp.fp_array = fp; fwp.norm_p_array = np;
	fwp.q_WP_array = qB; fwp.q_array = qA; fwp.wq_array = wq; fwp.fq_array = fq; fwp.norm_q_array = nq;

	PacketGridView grid(fwp);
	CHECK(grid.Np_WP() == 2 && grid.Nq_WP() == 3);
	CHECK(grid.p().num_packets() == 2 && grid.p().nodes_per_packet() == 1);
	CHECK(grid.p().lower(0) == 0.0 && grid.p().upper(0) == 1.0 && grid.p().width(1) == 1.0);
	CHECK(grid.p().node(1, 0) == 1.5);
	CHECK(grid.q().num_packets() == 3 && grid.q().nodes_per_packet() == 2);
	CHECK(grid.q().node(2, 1) == 2.5);
	CHECK(grid.q().upper(2) == 3.0);
	// dense_dim = Nalpha * Nq * Np
	CHECK(grid.dense_dim(/*nalpha=*/3) == static_cast<std::size_t>(3 * 3 * 2));
}

static void test_permutation_operator_view()
{
	using namespace tictac::core;
	// A 3x3 sparse matrix with 4 nonzeros in CSC:
	//   col 0: rows {0,1}; col 1: row {2}; col 2: (empty).
	double val[]  = {2.0, 3.0, 5.0};
	int    row[]  = {0,   1,   2};
	std::size_t cptr[] = {0, 2, 3, 3};
	PermutationOperatorView P(val, row, cptr, /*dim=*/3);
	CHECK(P.dimension() == 3);
	CHECK(P.num_nonzeros() == 3);

	auto c0 = P.column(0);
	CHECK(c0.begin == 0 && c0.end == 2);
	// iterate nonzeros of column 0
	double s0 = 0.0;
	for (std::size_t k = c0.begin; k < c0.end; ++k) s0 += c0.values[k];
	CHECK(std::fabs(s0 - 5.0) < 1e-12);

	auto c1 = P.column(1);
	CHECK(c1.begin == 2 && c1.end == 3);
	CHECK(std::fabs(c1.values[2] - 5.0) < 1e-12);

	auto c2 = P.column(2);
	CHECK(c2.begin == 3 && c2.end == 3);  // empty column
}

int main()
{
	test_three_body_channel_and_pair_antisymmetry();
	test_partial_wave_basis_view();
	test_packet_grid_view();
	test_permutation_operator_view();

	if (failures == 0) {
		std::printf("semantic_views: all checks passed\n");
		return 0;
	}
	std::printf("semantic_views: %d CHECK(s) FAILED\n", failures);
	return 1;
}
