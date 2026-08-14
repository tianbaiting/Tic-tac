#ifndef TICTAC_INTERACTIONS_W1_MANIFEST_H
#define TICTAC_INTERACTIONS_W1_MANIFEST_H

#include <cstddef>
#include <string>
#include <vector>

#include "type_defs.h"
#include "io/cache_layer/cache_keys.h"   // W1Key, hash_full, canonical_json
#include "w1_work_plan.h"

class three_nucleon_force_model;

namespace tictac::interactions {

// [EN] Provenance signature for a W^(1) construction campaign.
//
// This is the per-BLOCK cache key (make_w1_key) with the sector-/block-local
// fields (a_r, a_c, two_J_3N, P_3N) pinned to sentinels, so its hash changes iff
// any physics / grid / quadrature / schema input changes -- exactly the
// provenance the task requires.  `two_J_3NF_force_max` is recorded separately
// because it selects the active sector set but is not part of an individual
// block's value.
//
// Two campaigns are interchangeable ONLY if their signature hashes agree.  A
// stale manifest is detected by comparing hashes and never silently reused.
//
// / [CN] W^(1) 构造活动的来源签名。是单块 cache key 去掉扇区/块局部字段后的
// 模板；其哈希当且仅当物理/网格/求积/schema 输入变化时才变。two_J_3NF_force_max
// 单独记录（决定激活扇区集合，但不影响单块取值）。
struct W1Signature {
	tictac::cache::W1Key key_template;   // a_r=a_c=0, two_J_3N=P_3N=0 sentinels
	int two_J_3NF_force_max = -1;
};

// Build the signature from the model + run config + grid hashes.  Pure.
W1Signature make_w1_signature(const three_nucleon_force_model& tnf,
                              const run_params& rp,
                              const pw_3N_statespace& pw,
                              int Np_WP, int Nq_WP,
                              int Np_quad, int Nq_quad,
                              const std::string& p_grid_hash,
                              const std::string& q_grid_hash);

// Stable 64-hex SHA-256 over the canonical signature (reuses the cache layer's
// hash_full so manifest and per-block shard provenance cannot diverge).
std::string signature_hash(const W1Signature& sig);

// Human-readable canonical JSON of the signature (for the manifest file).
std::string signature_json(const W1Signature& sig);

// [EN] Per-sector slice of the work plan + a live completion snapshot.
// `num_evaluator_present` is the count of evaluate-units whose shard is a
// valid cache hit RIGHT NOW (re-derived at status/assemble time, never stored
// stale).  transpose_fill units are never persisted by workers; the solver
// produces them on demand, so completeness = (evaluate present == evaluate total).
// / [CN] 每扇区工作计划切片 + 当前完成快照。完成数实时派生，不存盘。
struct W1ManifestSector {
	int chn = 0;
	int two_J = 0;
	int parity = 0;
	int Nalpha = 0;
	bool active = false;
	std::size_t num_blocks = 0;
	std::size_t num_evaluate = 0;
	std::size_t num_transpose_fill = 0;
	std::size_t num_evaluate_present = 0;   // live cache hit count
	std::vector<W1WorkUnit> units;
};

// [EN] Full manifest report: signature + per-sector plan + live completion.
// This is what `plan`/`status`/`assemble` print and what `plan --manifest` writes.
// / [CN] 完整清单：签名 + 各扇区计划 + 实时完成情况。
struct W1ManifestReport {
	W1Signature signature;
	std::string signature_hash;
	std::vector<W1ManifestSector> sectors;
};

// Write a manifest report as deterministic JSON to `path` (atomic temp+rename,
// so an interrupted write never leaves a half-written manifest).
bool write_manifest_json(const std::string& path, const W1ManifestReport& report);

// [EN] Assembler: the canonical fingerprint of a completed campaign.
// Concatenates the dense double payloads of every evaluate-unit in deterministic
// (sector, a_r, a_c) order and returns the SHA-256 over the raw bytes.  This is
// the exactness comparator: monolithic vs resumable vs multi-worker builds all
// produce the same fingerprint iff they produce bitwise-identical W^(1).
// Missing evaluate-units make the assembler return false (incomplete).
// / [CN] 装配器：完成活动的规范指纹。按确定顺序拼接每块双精度载荷并求 SHA-256，
// 作为精确性比较器。缺块时返回 false。
struct W1Assembler {
	// Load every evaluate-unit shard from the cache and compute the fingerprint.
	// `missing_out` (if non-null) receives the first missing (chn, a_r, a_c).
	static bool fingerprint(const three_nucleon_force_model& tnf,
	                        const run_params& rp,
	                        const pw_3N_statespace& pw_global,
	                        int Np_WP, int Nq_WP,
	                        int Np_quad, int Nq_quad,
	                        const std::string& p_grid_hash,
	                        const std::string& q_grid_hash,
	                        std::string& hash_out,
	                        std::string* missing_out);
};

} // namespace tictac::interactions

#endif // TICTAC_INTERACTIONS_W1_MANIFEST_H
