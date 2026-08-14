// W1Manifest + W1Assembler implementation.
//
// Reuses the per-block cache key (make_w1_key) and the cache layer's SHA-256
// hashing, so the manifest's provenance cannot diverge from the per-block shard
// provenance.  The assembler's fingerprint is a streaming SHA-256 of per-block
// payload hashes -- memory-light and deterministic, so monolithic vs resumable
// vs multi-worker builds compare exactly.

#include "w1_manifest.h"

#include "w1_pw_cache.h"           // make_w1_key
#include "three_nucleon_force_model.h"
#include "io/cache_layer/cache_layer.h"
#include "io/cache_layer/cache_io_w1.h"

#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>   // ::rename
#include <utility>

namespace tictac::interactions {

namespace {

std::string json_escape(const std::string& s)
{
	std::string out; out.reserve(s.size() + 2);
	out.push_back('"');
	for (char c : s) {
		switch (c) {
			case '"':  out += "\\\""; break;
			case '\\': out += "\\\\"; break;
			case '\n': out += "\\n";  break;
			case '\r': out += "\\r";  break;
			case '\t': out += "\\t";  break;
			default:
				if ((unsigned char)c < 0x20) {
					char hex[8]; std::snprintf(hex, sizeof(hex), "\\u%04x", c);
					out += hex;
				} else out.push_back(c);
		}
	}
	out.push_back('"');
	return out;
}

const char* role_str(W1UnitRole r) { return r == W1UnitRole::evaluate ? "evaluate" : "transpose_fill"; }

} // namespace

W1Signature make_w1_signature(const three_nucleon_force_model& tnf,
                              const run_params& rp,
                              const pw_3N_statespace& pw,
                              int Np_WP, int Nq_WP,
                              int Np_quad, int Nq_quad,
                              const std::string& p_grid_hash,
                              const std::string& q_grid_hash)
{
	// Build a representative per-block key, then pin the sector-/block-local
	// fields to sentinels so the hash only varies with physics/grid/quadrature.
	W1Signature sig;
	sig.key_template = make_w1_key(tnf, rp, pw, Np_WP, Nq_WP, Np_quad, Nq_quad,
	                               /*a_r*/ 0, /*a_c*/ 0, p_grid_hash, q_grid_hash);
	sig.key_template.a_r = 0;
	sig.key_template.a_c = 0;
	sig.key_template.two_J_3N = 0;   // sentinel: sector-independent
	sig.key_template.P_3N = 0;       // sentinel: sector-independent
	sig.two_J_3NF_force_max = rp.two_J_3NF_force_max;
	return sig;
}

std::string signature_hash(const W1Signature& sig)
{
	// Reuse the cache layer's canonical-JSON + SHA-256 so the manifest hash and
	// the per-block shard hash share one provenance pipeline.  Append
	// two_J_3NF_force_max (active-sector selector, not in W1Key) into the hash
	// input so a cutoff change is detected too.
	const std::string base = tictac::cache::canonical_json(sig.key_template)
	                       + "|two_J_3NF_force_max=" + std::to_string(sig.two_J_3NF_force_max);
	return tictac::cache::hash_full_raw(
	    reinterpret_cast<const unsigned char*>(base.data()), base.size());
}

std::string signature_json(const W1Signature& sig)
{
	return tictac::cache::canonical_json(sig.key_template)
	     + ",\"two_J_3NF_force_max\":" + std::to_string(sig.two_J_3NF_force_max);
}

bool write_manifest_json(const std::string& path, const W1ManifestReport& report)
{
	// Atomic publication: write to a sibling temp file, then rename.  An
	// interrupted write never leaves a half-written manifest at `path`.
	const std::string tmp = path + ".tmp";
	{
		std::ofstream os(tmp, std::ios::binary | std::ios::trunc);
		if (!os) return false;
		os << "{\n";
		os << "  \"manifest_schema\": 1,\n";
		os << "  \"signature\": " << signature_json(report.signature) << ",\n";
		os << "  \"signature_hash\": " << json_escape(report.signature_hash) << ",\n";
		os << "  \"two_J_3NF_force_max\": " << report.signature.two_J_3NF_force_max << ",\n";
		os << "  \"sectors\": [\n";
		bool first_sector = true;
		for (const auto& s : report.sectors) {
			if (!first_sector) os << ",\n";
			first_sector = false;
			os << "    {\n";
			os << "      \"chn\": " << s.chn << ", \"two_J\": " << s.two_J
			   << ", \"parity\": " << s.parity << ", \"Nalpha\": " << s.Nalpha
			   << ", \"active\": " << (s.active ? "true" : "false") << ",\n";
			os << "      \"num_blocks\": " << s.num_blocks
			   << ", \"num_evaluate\": " << s.num_evaluate
			   << ", \"num_transpose_fill\": " << s.num_transpose_fill
			   << ", \"num_evaluate_present\": " << s.num_evaluate_present << ",\n";
			os << "      \"units\": [\n";
			bool first_unit = true;
			for (const auto& u : s.units) {
				if (!first_unit) os << ",\n";
				first_unit = false;
				os << "        {\"id\": " << u.id
				   << ", \"a_r\": " << u.a_r << ", \"a_c\": " << u.a_c
				   << ", \"role\": " << json_escape(role_str(u.role));
				if (u.role == W1UnitRole::transpose_fill) {
					os << ", \"conj_a_r\": " << u.conj_a_r
					   << ", \"conj_a_c\": " << u.conj_a_c;
				}
				os << "}";
			}
			os << "\n      ]\n";
			os << "    }";
		}
		os << "\n  ]\n";
		os << "}\n";
		if (!os) { std::remove(tmp.c_str()); return false; }
	}
	return std::rename(tmp.c_str(), path.c_str()) == 0;
}

bool W1Assembler::fingerprint(const three_nucleon_force_model& tnf,
                              const run_params& rp,
                              const pw_3N_statespace& pw_global,
                              int Np_WP, int Nq_WP,
                              int Np_quad, int Nq_quad,
                              const std::string& p_grid_hash,
                              const std::string& q_grid_hash,
                              std::string& hash_out,
                              std::string* missing_out)
{
	// Streaming fingerprint: for every active sector, for every evaluate-unit
	// in deterministic (a_r,a_c) order, look up its shard and feed the SHA-256
	// of the raw payload into the running hash.  Memory-light and exactly
	// reproducible: bitwise-identical campaigns yield identical fingerprints.
	std::string accumulator;
	for (int chn = 0; chn < pw_global.N_chn_3N; ++chn) {
		pw_3N_statespace sub = make_channel_view(pw_global, chn);
		const int two_J = sub.two_J_3N_array[0];
		const bool active = (rp.two_J_3NF_force_max < 0 || two_J <= rp.two_J_3NF_force_max);
		if (!active) continue;

		const bool hermitian = tnf.W1_is_exactly_hermitian();
		W1WorkPlan plan(sub, hermitian);
		for (const W1WorkUnit* u : plan.evaluate_units()) {
			const auto key = make_w1_key(tnf, rp, sub, Np_WP, Nq_WP,
			                             Np_quad, Nq_quad, u->a_r, u->a_c,
			                             p_grid_hash, q_grid_hash);
			tictac::cache::W1Block blk{};
			const auto res = tictac::cache::lookup_w1(key, &blk);
			if (!res.hit
			    || blk.Np != Np_WP || blk.Nq != Nq_WP
			    || blk.data.size() != static_cast<std::size_t>(Nq_WP) * Nq_WP * Np_WP * Np_WP)
			{
				if (missing_out) {
					*missing_out = "chn=" + std::to_string(chn)
					             + " a_r=" + std::to_string(u->a_r)
					             + " a_c=" + std::to_string(u->a_c)
					             + " reason=" + (res.hit ? "shape_mismatch" : res.miss_reason);
				}
				return false;
			}
			const std::size_t nbytes = blk.data.size() * sizeof(double);
			accumulator += tictac::cache::hash_full_raw(
			    reinterpret_cast<const unsigned char*>(blk.data.data()), nbytes);
		}
	}
	hash_out = tictac::cache::hash_full_raw(
	    reinterpret_cast<const unsigned char*>(accumulator.data()),
	    accumulator.size());
	return true;
}

} // namespace tictac::interactions
