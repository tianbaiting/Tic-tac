#ifndef TICTAC_INTERACTIONS_W1_BLOCK_STORE_H
#define TICTAC_INTERACTIONS_W1_BLOCK_STORE_H

#include <cstddef>
#include <map>
#include <optional>
#include <string>
#include <vector>

namespace tictac::interactions {

// [EN] Identity of one independently-addressable, exactly-computed W^(1) block.
//
// Within a fixed basis/grid/coupling context (captured by W1BlockSignature), the
// expensive W^(1) construction factors into independent dense blocks, one per
// allowed (alpha_r, alpha_c) channel pair of a conserved J^pi block.  A block
// built by any worker is reusable by any solver once its identity and signature
// match, which is the foundation of resumable / distributed EXACT construction.
//
// Required properties (docs/complete_n2lo_3nf_status.md, Phase 7):
//   - exact results only (no low-rank, no compression);
//   - deterministic block identity;
//   - safe restart after interruption (a completed block is never recomputed);
//   - Hermitian-transpose reuse remains EXACT and explicit (see W1Block).
// / [CN] 一个可独立寻址、精确计算的 W^(1) 块的身份标识。
struct W1BlockId {
	int alpha_r = -1;   // bra partial-wave channel (row)
	int alpha_c = -1;   // ket partial-wave channel (column)

	bool operator<(const W1BlockId& o) const
	{ return alpha_r != o.alpha_r ? alpha_r < o.alpha_r : alpha_c < o.alpha_c; }
	bool operator==(const W1BlockId& o) const
	{ return alpha_r == o.alpha_r && alpha_c == o.alpha_c; }
};

// [EN] Every physics/numerical input that affects a W^(1) block value. Two
// blocks are interchangeable ONLY if both id and signature match. This mirrors
// the existing hash-keyed HDF5 cache identity (src/io/cache_layer/cache_keys)
// and is what makes cross-worker, cross-run reuse safe.
// / [CN] 所有影响 W^(1) 块取值的物理/数值输入；只有 id 与签名都匹配的块才可互换。
struct W1BlockSignature {
	std::string  tnf_model;          // e.g. "chiral_N2LO_full_factorized"
	int          two_J_3N      = 0;  // conserved 2J of the block
	int          two_T_3N      = 0;
	int          parity        = 0;
	int          J_2N_max      = 0;
	std::size_t  Np_WP         = 0;
	std::size_t  Nq_WP         = 0;
	int          Np_per_WP_W1  = 0;  // radial quadrature order (p)
	int          Nq_per_WP_W1  = 0;  // radial quadrature order (q)
	int          Nangle_3NF    = 0;  // 5D angular projector order
	double       c_1 = 0, c_3 = 0, c_4 = 0, c_D = 0, c_E = 0;
	double       Lambda_3NF = 0, Lambda_chi = 0, g_A = 0, f_pi = 0, m_pi = 0, hbarc = 0;
	std::string  p_grid_hash;        // exact p-boundary array hash
	std::string  q_grid_hash;        // exact q-boundary array hash
	std::string  regulator_kind = "gaussian";
	int          schema_version = 0;
};

// [EN] One exactly-computed dense W^(1) block payload. The dense layout matches
// W1_PW_cache: index = (((blk*Nq + q_r)*Nq + q_c)*Np + p_r)*Np + p_c.
//
// Hermitian-transpose reuse is EXACT and EXPLICIT: if this block was produced
// by transposing its conjugate (alpha_c, alpha_r) rather than direct
// integration, `is_transpose_reuse` is set and `conjugate_id` records the
// source. This keeps the manifest honest about which orientations were
// integrated versus transposed (contract §3.1 / status doc).
// / [CN] 一个精确稠密 W^(1) 块载荷。Hermitian 转置复用是显式且精确的。
struct W1Block {
	W1BlockId       id;
	std::vector<double> data;        // dense payload, size Nq*Nq*Np*Np
	std::size_t     Np = 0;
	std::size_t     Nq = 0;
	bool            is_transpose_reuse = false;
	W1BlockId       conjugate_id;    // populated iff is_transpose_reuse
};

// [EN] Abstract block store: persistence and retrieval of exact W^(1) blocks.
// Concrete stores (in-memory, HDF5, future distributed) implement this so a
// builder and a solver are decoupled from WHERE a block was computed.
// / [CN] 抽象块存储：精确 W^(1) 块的持久化与检索。
class W1BlockStore {
public:
	virtual ~W1BlockStore() = default;

	// Has this exact (id, signature) block already been stored?
	virtual bool contains(const W1BlockId& id,
	                      const W1BlockSignature& sig) const = 0;

	// Load a stored block. Returns std::nullopt if absent (caller may then build
	// and save() it -- the resumable-construction loop).
	virtual std::optional<W1Block> load(const W1BlockId& id,
	                                    const W1BlockSignature& sig) const = 0;

	// Persist an exactly-computed block for future contains()/load().
	virtual void save(const W1Block& block,
	                  const W1BlockSignature& sig) = 0;

	// Completeness manifest: the set of block ids present for a signature.
	virtual std::vector<W1BlockId> manifest(const W1BlockSignature& sig) const = 0;
};

// [EN] In-memory block store. Backs the resumable loop for a single process;
// the existing W1_PW_cache (dense all-blocks-at-once cache) and the hash-keyed
// HDF5 cache (src/io/cache_layer) are the production concrete stores and remain
// the numerically authoritative paths.
// / [CN] 内存块存储；生产仍用 W1_PW_cache 与 HDF5 缓存。
class MemoryW1BlockStore : public W1BlockStore {
public:
	using SigKey = std::string;  // serialized signature

	bool contains(const W1BlockId& id, const W1BlockSignature& sig) const override;
	std::optional<W1Block> load(const W1BlockId& id,
	                            const W1BlockSignature& sig) const override;
	void save(const W1Block& block, const W1BlockSignature& sig) override;
	std::vector<W1BlockId> manifest(const W1BlockSignature& sig) const override;

	// Deterministic signature serialization (stable ordering for cross-run use).
	static SigKey serialize(const W1BlockSignature& s);

private:
	// signature -> (block id -> block)
	std::map<SigKey, std::map<W1BlockId, W1Block>> m_store;
};

} // namespace tictac::interactions

#endif // TICTAC_INTERACTIONS_W1_BLOCK_STORE_H
