// W1 block-store abstraction test (Phase 7).
//
// Exercises the resumable/distributed exact-construction contract:
//   - contains() is false before save(), true after;
//   - load() round-trips the payload exactly;
//   - a completed block is not "lost" (restart safety);
//   - manifest() is deterministic and signature-scoped;
//   - Hermitian-transpose reuse is recorded explicitly.
// Returns 0 on success.

#include <cstdio>
#include <vector>

#include "w1_block_store.h"

static int failures = 0;
#define CHECK(cond) do { if (!(cond)) { ++failures; \
	std::printf("FAIL: %s (line %d)\n", #cond, __LINE__); } } while (0)

int main()
{
	using namespace tictac::interactions;

	W1BlockSignature sig{};
	sig.tnf_model = "chiral_N2LO_full_factorized";
	sig.two_J_3N = 1; sig.parity = +1; sig.Np_WP = 4; sig.Nq_WP = 3;
	sig.c_D = -0.2; sig.c_E = -0.205; sig.schema_version = 9;

	MemoryW1BlockStore store;
	W1BlockId id0{0, 0};
	W1BlockId id1{0, 1};

	// 1. Resumable loop: miss -> build -> save -> hit.
	CHECK(!store.contains(id0, sig));
	W1Block b0;
	b0.id = id0;
	b0.Np = 4; b0.Nq = 3;
	b0.data.assign(/*Nq*Nq*Np*Np*/ 3 * 3 * 4 * 4, 1.23456789);
	store.save(b0, sig);
	CHECK(store.contains(id0, sig));

	auto loaded = store.load(id0, sig);
	CHECK(loaded.has_value());
	CHECK(loaded->data.size() == b0.data.size());
	bool exact = true;
	for (std::size_t i = 0; i < b0.data.size(); ++i)
		if (loaded->data[i] != b0.data[i]) { exact = false; break; }
	CHECK(exact);  // byte-exact payload round-trip (exact, no compression)

	// 2. Restart safety: a second block built later joins the same store.
	CHECK(!store.contains(id1, sig));
	W1Block b1;
	b1.id = id1; b1.Np = 4; b1.Nq = 3;
	b1.is_transpose_reuse = true;     // explicit Hermitian transpose of (1,0)
	const W1BlockId conj_src{1, 0};
	b1.conjugate_id = conj_src;
	b1.data.assign(3 * 3 * 4 * 4, -0.5);
	store.save(b1, sig);

	auto mf = store.manifest(sig);
	CHECK(mf.size() == 2);
	CHECK(mf[0] == id0 && mf[1] == id1);  // deterministic (alpha_r,alpha_c) order

	auto loaded1 = store.load(id1, sig);
	CHECK(loaded1.has_value());
	CHECK(loaded1->is_transpose_reuse);
	const W1BlockId conj{1, 0};
	CHECK(loaded1->conjugate_id == conj);

	// 3. Signature scoping: a different schema/LEC does not see the blocks.
	W1BlockSignature sig_other = sig;
	sig_other.c_E = 0.0;            // different physics -> different blocks
	CHECK(!store.contains(id0, sig_other));
	CHECK(store.manifest(sig_other).empty());

	// Deterministic serialization: identical signatures -> identical key.
	CHECK(MemoryW1BlockStore::serialize(sig) == MemoryW1BlockStore::serialize(sig));

	if (failures == 0) {
		std::printf("w1_block_store: resumable/transpose-reuse contract OK\n");
		return 0;
	}
	std::printf("w1_block_store: %d CHECK(s) FAILED\n", failures);
	return 1;
}
