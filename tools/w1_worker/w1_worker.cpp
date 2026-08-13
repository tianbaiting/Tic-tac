// W1 block-database worker (Phase D/E/F).
//
// Turns exact W1 construction into a resumable, distributable block campaign.
// W1 is energy-independent and block-addressable (per (two_J_3N, P_3N, a_r, a_c));
// this tool enumerates the blocks needed for the active J_3NF sectors, reports
// which are already stored, builds the missing ones (skipping cache hits), and
// verifies the result -- without touching the Faddeev solve.
//
//   w1_worker plan   input.txt [--shard K/N]            list sectors + completeness
//   w1_worker build  input.txt [--sector J P] [--shard K/N]  build missing blocks
//   w1_worker verify input.txt [--shard K/N]            load-check every block
//
// Sectors above two_J_3NF_force_max are pure-2NF and have NO W1 blocks; they are
// never planned/built. Lower-J W1 built under one run is reused verbatim under a
// larger J_3NF_max (the J ladder reuses completed sectors).

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "set_run_parameters.h"
#include "make_pw_symm_states.h"
#include "make_wp_states.h"
#include "w1_pw_cache.h"
#include "three_nucleon_force_model.h"
#include "io/cache_layer/cache_layer.h"
#include "io/cache_layer/cache_keys.h"
#include "io/cache_layer/cache_manifest.h"

// make_channel_view: slice one conserved J^pi block out of the global basis
// (same logic as solver_pipeline; kept local so the worker has no pipeline dep).
static pw_3N_statespace channel_view(const pw_3N_statespace& pw, int chn) {
	const int lo = pw.chn_3N_idx_array[chn];
	const int hi = pw.chn_3N_idx_array[chn + 1];
	pw_3N_statespace v = {};
	v.Nalpha = hi - lo;
	v.J_2N_max = pw.J_2N_max;
	v.L_2N_array     = &pw.L_2N_array[lo];
	v.S_2N_array     = &pw.S_2N_array[lo];
	v.J_2N_array     = &pw.J_2N_array[lo];
	v.T_2N_array     = &pw.T_2N_array[lo];
	v.L_1N_array     = &pw.L_1N_array[lo];
	v.two_J_1N_array = &pw.two_J_1N_array[lo];
	v.two_T_3N_array = &pw.two_T_3N_array[lo];
	v.two_J_3N_array = &pw.two_J_3N_array[lo];
	v.P_3N_array     = &pw.P_3N_array[lo];
	return v;
}

struct SectorInfo {
	int chn;          // sector index in [0, N_chn_3N)
	int two_J;
	int parity;
	int num_blocks;   // Nalpha_s^2
	bool active;      // two_J <= two_J_3NF_force_max (or cutoff<0)
};

static bool parse_shard(const char* s, int& k, int& n) {
	// "K/N"
	const char* slash = std::strchr(s, '/');
	if (!slash) return false;
	k = std::atoi(s);
	n = std::atoi(slash + 1);
	return n > 0 && k >= 0 && k < n;
}

int main(int argc, char** argv) {
	if (argc < 3) {
		std::fprintf(stderr,
		    "usage: w1_worker {plan|build|verify} input.txt [--sector J P] [--shard K/N]\n");
		return 2;
	}
	const std::string cmd = argv[1];

	run_params rp{};
	// Build a synthetic argv that the existing parser consumes: prog input.txt ...
	std::vector<char*> av(argv, argv + argc);
	set_run_parameters(argc, av.data(), rp);

	if (rp.three_nucleon_force == "none" || rp.three_nucleon_force.empty()) {
		std::fprintf(stderr, "w1_worker: three_nucleon_force='%s' -> no W1 to build.\n",
		             rp.three_nucleon_force.c_str());
		return 0;
	}

	// Optional --sector J P / --shard K/N
	int sel_two_J = -1, sel_P = 0; bool have_sector = false;
	int shard_k = 0, shard_n = 1;
	for (int i = 3; i < argc; ++i) {
		std::string a = argv[i];
		if (a == "--sector" && i + 2 < argc) {
			sel_two_J = std::atoi(argv[++i]);
			sel_P = std::atoi(argv[++i]);
			have_sector = true;
		} else if (a == "--shard" && i + 1 < argc) {
			if (!parse_shard(argv[++i], shard_k, shard_n)) {
				std::fprintf(stderr, "bad --shard K/N\n"); return 2;
			}
		}
	}

	pw_3N_statespace pw{};
	construct_symmetric_pw_states(pw, rp);
	fwp_statespace fwp{};
	make_fwp_statespace(fwp, rp);
#if TICTAC_USE_NEW_CACHE_LAYER
	tictac::cache::initialize(rp.cache_root);
#endif
	auto tnf = three_nucleon_force_model::create(rp);

	// Grid hashes (for block-key construction / probing).
	std::string p_hash, q_hash;
#if TICTAC_USE_NEW_CACHE_LAYER
	p_hash = tictac::cache::hash_full_raw(
	    reinterpret_cast<const unsigned char*>(fwp.p_WP_array), (fwp.Np_WP + 1) * sizeof(double));
	q_hash = tictac::cache::hash_full_raw(
	    reinterpret_cast<const unsigned char*>(fwp.q_WP_array), (fwp.Nq_WP + 1) * sizeof(double));
#endif
	const int Np_quad = std::max(1, rp.Np_per_WP_W1);
	const int Nq_quad = std::max(1, rp.Nq_per_WP_W1);

	// Enumerate sectors + active flag.
	std::vector<SectorInfo> sectors;
	for (int c = 0; c < pw.N_chn_3N; ++c) {
		SectorInfo s;
		s.chn = c;
		pw_3N_statespace sub = channel_view(pw, c);
		s.two_J = sub.two_J_3N_array[0];
		s.parity = sub.P_3N_array[0];
		s.num_blocks = sub.Nalpha * sub.Nalpha;
		s.active = (rp.two_J_3NF_force_max < 0 || s.two_J <= rp.two_J_3NF_force_max);
		sectors.push_back(s);
	}

	// Manifest for fast completeness counts.
#if TICTAC_USE_NEW_CACHE_LAYER
	tictac::cache::Manifest manifest = tictac::cache::Manifest::load(rp.cache_root);
#endif

	// Count stored blocks per active sector via the manifest (fast, no HDF5 read).
	auto stored_in_sector = [&](const SectorInfo& s) -> long {
#if TICTAC_USE_NEW_CACHE_LAYER
		long n = 0;
		pw_3N_statespace sub = channel_view(pw, s.chn);
		for (int a_r = 0; a_r < sub.Nalpha; ++a_r)
			for (int a_c = 0; a_c < sub.Nalpha; ++a_c) {
				auto key = make_w1_key(*tnf, rp, sub, fwp.Np_WP, fwp.Nq_WP,
				                       Np_quad, Nq_quad, a_r, a_c, p_hash, q_hash);
				if (tictac::cache::Manifest::load(rp.cache_root).find(tictac::cache::hash_full(key)))
					++n;
			}
		return n;
#else
		return 0;
#endif
	};

	// Shard filter: sector belongs to shard k of n iff (chn % n) == k.
	auto in_shard = [&](const SectorInfo& s) {
		return (s.chn % shard_n) == shard_k;
	};
	// Sector filter.
	auto selected = [&](const SectorInfo& s) {
		if (!s.active) return false;
		if (have_sector && !(s.two_J == sel_two_J && s.parity == sel_P)) return false;
		if (!in_shard(s)) return false;
		return true;
	};

	if (cmd == "plan") {
		std::printf("# W1 block plan: model=%s Np=%d Nq=%d two_J_3N_max=%d two_J_3NF_force_max=%d\n",
		            tnf->name().c_str(), rp.Np_WP, rp.Nq_WP, rp.two_J_3N_max, rp.two_J_3NF_force_max);
		std::printf("# %-6s %-6s %-10s %-12s %-12s %-12s\n",
		            "chn", "two_J", "parity", "num_blocks", "stored", "active");
		long tot_blocks = 0, tot_stored = 0;
		for (const auto& s : sectors) {
			const long stored = s.active ? stored_in_sector(s) : 0;
			tot_blocks += s.active ? s.num_blocks : 0;
			tot_stored += stored;
			std::printf("  %-6d %-6d %-10d %-12d %-12ld %-12s\n",
			            s.chn, s.two_J, s.parity, s.num_blocks, stored,
			            s.active ? "yes" : "no(2NF)");
		}
		std::printf("# active blocks: %ld, stored: %ld, missing: %ld\n",
		            tot_blocks, tot_stored, tot_blocks - tot_stored);
		return 0;
	}

	if (cmd == "build") {
		int built_sectors = 0;
		for (const auto& s : sectors) {
			if (!selected(s)) continue;
			std::printf("[w1_worker] building sector chn=%d J=%d/2 P=%d (Nalpha^2=%d blocks) ...\n",
			            s.chn, s.two_J, s.parity, s.num_blocks);
			pw_3N_statespace sub = channel_view(pw, s.chn);
			W1_PW_cache cache;
			// build() probes the cache per block: hits are skipped, only missing
			// blocks are evaluated and atomically stored. Resumable by construction.
			cache.build(*tnf, fwp.p_WP_array, fwp.Np_WP,
			            fwp.q_WP_array, fwp.Nq_WP, sub, rp);
			++built_sectors;
		}
		std::printf("[w1_worker] built %d sector(s).\n", built_sectors);
		return 0;
	}

	if (cmd == "verify") {
		int bad = 0; long checked = 0;
		for (const auto& s : sectors) {
			if (!selected(s)) continue;
			pw_3N_statespace sub = channel_view(pw, s.chn);
			for (int a_r = 0; a_r < sub.Nalpha; ++a_r)
				for (int a_c = 0; a_c < sub.Nalpha; ++a_c) {
#if TICTAC_USE_NEW_CACHE_LAYER
					auto key = make_w1_key(*tnf, rp, sub, fwp.Np_WP, fwp.Nq_WP,
					                       Np_quad, Nq_quad, a_r, a_c, p_hash, q_hash);
					tictac::cache::W1Block blk;
					auto res = tictac::cache::lookup_w1(key, &blk);
					if (!res.hit) { ++bad; }
					++checked;
#endif
				}
		}
		std::printf("[w1_worker] verify: %ld blocks checked, %d missing/invalid\n", checked, bad);
		return bad == 0 ? 0 : 1;
	}

	std::fprintf(stderr, "unknown command '%s' (use plan|build|verify)\n", cmd.c_str());
	return 2;
}
