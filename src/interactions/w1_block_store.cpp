// W1 block-store abstraction implementation (Phase 7).
//
// Concrete in-memory store + deterministic signature serialization. This is the
// resumable/distributed exact-construction design layer; the production W1 path
// still goes through W1_PW_cache and the hash-keyed HDF5 cache (unchanged).

#include "w1_block_store.h"

#include <sstream>

namespace tictac::interactions {

MemoryW1BlockStore::SigKey MemoryW1BlockStore::serialize(const W1BlockSignature& s)
{
	// Deterministic, field-by-field serialization with explicit delimiters so two
	// equal signatures always produce identical keys (cross-run / cross-worker).
	std::ostringstream os;
	os << s.tnf_model << '|'
	   << s.two_J_3N << '|' << s.two_T_3N << '|' << s.parity << '|' << s.J_2N_max
	   << '|' << s.Np_WP << '|' << s.Nq_WP
	   << '|' << s.Np_per_WP_W1 << '|' << s.Nq_per_WP_W1 << '|' << s.Nangle_3NF
	   << '|' << s.c_1 << '|' << s.c_3 << '|' << s.c_4
	   << '|' << s.c_D << '|' << s.c_E
	   << '|' << s.Lambda_3NF << '|' << s.Lambda_chi
	   << '|' << s.g_A << '|' << s.f_pi << '|' << s.m_pi << '|' << s.hbarc
	   << '|' << s.p_grid_hash << '|' << s.q_grid_hash
	   << '|' << s.regulator_kind << '|' << s.schema_version;
	return os.str();
}

bool MemoryW1BlockStore::contains(const W1BlockId& id, const W1BlockSignature& sig) const
{
	const auto sig_it = m_store.find(serialize(sig));
	if (sig_it == m_store.end()) return false;
	return sig_it->second.find(id) != sig_it->second.end();
}

std::optional<W1Block> MemoryW1BlockStore::load(const W1BlockId& id,
                                                const W1BlockSignature& sig) const
{
	const auto sig_it = m_store.find(serialize(sig));
	if (sig_it == m_store.end()) return std::nullopt;
	const auto blk_it = sig_it->second.find(id);
	if (blk_it == sig_it->second.end()) return std::nullopt;
	return blk_it->second;
}

void MemoryW1BlockStore::save(const W1Block& block, const W1BlockSignature& sig)
{
	m_store[serialize(sig)][block.id] = block;
}

std::vector<W1BlockId> MemoryW1BlockStore::manifest(const W1BlockSignature& sig) const
{
	std::vector<W1BlockId> ids;
	const auto sig_it = m_store.find(serialize(sig));
	if (sig_it != m_store.end()) {
		ids.reserve(sig_it->second.size());
		for (const auto& kv : sig_it->second) ids.push_back(kv.first);
	}
	return ids;  // std::map ordering -> deterministic (alpha_r, then alpha_c)
}

} // namespace tictac::interactions
