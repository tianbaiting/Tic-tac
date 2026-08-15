// W1 block-database worker (Phase D/E/F/G + block-level distribution).
//
// Turns exact W1 construction into a resumable, distributable block campaign.
// W1 is energy-independent and block-addressable (per (two_J_3N, P_3N, a_r, a_c));
// this tool enumerates the blocks needed for the active J_3NF sectors, reports
// which are already stored, builds the missing ones (skipping cache hits), and
// verifies / assembles the result -- without touching the Faddeev solve.
//
// Two complementary distribution modes:
//
//   --shard K/N              sector-level: worker K of N owns sectors with
//                            (chn % N) == K.  Each owned sector is built
//                            monolithically via W1_PW_cache::build (Hermitian
//                            triangle + transpose-fill + store).  Back-compat.
//
//   --worker-index I --worker-count N   block-level: worker I of N owns the
//                            evaluate-units (Hermitian triangle only) whose
//                            global evaluate-index satisfies idx % N == I.
//                            Each owned unit is integrated by W1BlockExecutor
//                            (the same shared per-cell integration as the
//                            monolithic path -> bitwise-identical shards) and
//                            published atomically.  transpose_fill units are
//                            NEVER built by workers; the solver produces them
//                            on demand.  This lets N workers attack one large
//                            sector in parallel.
//
// Commands:
//   w1_worker plan    input.txt [--manifest path] [--sector J P] [--shard K/N]
//   w1_worker build   input.txt [--sector J P] ( --shard K/N | --worker-index I --worker-count N | --blocks LIST )
//   w1_worker status  input.txt [--manifest path]
//   w1_worker assemble input.txt [--manifest path]
//   w1_worker verify  input.txt [--sector J P] [--shard K/N]
//
// --blocks LIST: comma-separated a_r:a_c pairs (e.g. "0:0,1:3,8:9").  Requires
//   --sector J P (exactly one sector).  Builds ONLY the listed evaluate-units in
//   that sector, each via its own W1BlockExecutor::compute_block call (the SAME
//   integrate_w1_channel_blocks path as production — bitwise-identical).  Prints
//   one JSON line per block to stdout for machine-readable timing.  Intended for
//   realistic-grid pilot measurements; NOT a second integration path.
//
// Sectors above two_J_3NF_force_max are pure-2NF and have NO W1 blocks; they are
// never planned/built. Lower-J W1 built under one run is reused verbatim under a
// larger J_3N_max (the J ladder reuses completed sectors).

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

#include "set_run_parameters.h"
#include "make_pw_symm_states.h"
#include "make_wp_states.h"
#include "w1_pw_cache.h"
#include "w1_integrate.h"
#include "w1_work_plan.h"
#include "w1_manifest.h"
#include "three_nucleon_force_model.h"
#include "io/cache_layer/cache_layer.h"
#include "io/cache_layer/cache_keys.h"
#include "io/cache_layer/cache_manifest.h"

using tictac::interactions::W1WorkPlan;
using tictac::interactions::W1WorkUnit;
using tictac::interactions::W1UnitRole;
using tictac::interactions::W1BlockExecutor;
using tictac::interactions::integrate_w1_channel_blocks;
using tictac::interactions::make_channel_view;
using tictac::interactions::make_w1_signature;
using tictac::interactions::signature_hash;
using tictac::interactions::write_manifest_json;
using tictac::interactions::W1ManifestReport;
using tictac::interactions::W1ManifestSector;
using tictac::cache::W1Block;
using tictac::cache::W1Key;

struct SectorInfo {
	int chn;
	int two_J;
	int parity;
	int Nalpha;
	bool active;   // two_J <= two_J_3NF_force_max (or cutoff<0)
};

static bool parse_shard(const char* s, int& k, int& n) {
	const char* slash = std::strchr(s, '/');
	if (!slash) return false;
	k = std::atoi(s);
	n = std::atoi(slash + 1);
	return n > 0 && k >= 0 && k < n;
}

struct WorkerArgs {
	std::string cmd;
	int sel_two_J = -1, sel_P = 0; bool have_sector = false;
	int shard_k = 0, shard_n = 1;                       // sector-level
	int worker_index = 0, worker_count = 1;            // block-level
	bool block_partition = false;
	std::string manifest_path;
	bool have_manifest = false;
	// --blocks "a_r:a_c,a_r:a_c,..." : explicit list of evaluate-units to build.
	// Requires --sector J P (exactly one sector). Empty = not selected.
	std::vector<std::pair<int,int>> selected_blocks;
	bool have_selected_blocks = false;
};

// Is this sector in this worker's sector-level shard?
static bool sector_in_shard(const SectorInfo& s, const WorkerArgs& a) {
	return (s.chn % a.shard_n) == a.shard_k;
}

static bool sector_selected(const SectorInfo& s, const WorkerArgs& a) {
	if (!s.active) return false;
	if (a.have_sector && !(s.two_J == a.sel_two_J && s.parity == a.sel_P)) return false;
	if (!a.block_partition && !sector_in_shard(s, a)) return false;  // sector-level shard
	return true;
}

int main(int argc, char** argv) {
	if (argc < 3) {
		std::fprintf(stderr,
		    "usage: w1_worker {plan|build|status|assemble|verify} input.txt "
		    "[--sector J P] [--shard K/N] [--worker-index I --worker-count N] [--manifest path]\n");
		return 2;
	}
	WorkerArgs args;
	args.cmd = argv[1];

	run_params rp{};
	std::vector<char*> av(argv, argv + argc);
	set_run_parameters(argc, av.data(), rp);

	if (rp.three_nucleon_force == "none" || rp.three_nucleon_force.empty()) {
		std::fprintf(stderr, "w1_worker: three_nucleon_force='%s' -> no W1 to build.\n",
		             rp.three_nucleon_force.c_str());
		return 0;
	}

	for (int i = 3; i < argc; ++i) {
		std::string a = argv[i];
		if (a == "--sector" && i + 2 < argc) {
			args.sel_two_J = std::atoi(argv[++i]);
			args.sel_P = std::atoi(argv[++i]);
			args.have_sector = true;
		} else if (a == "--shard" && i + 1 < argc) {
			if (!parse_shard(argv[++i], args.shard_k, args.shard_n)) {
				std::fprintf(stderr, "bad --shard K/N\n"); return 2;
			}
		} else if (a == "--worker-index" && i + 1 < argc) {
			args.worker_index = std::atoi(argv[++i]);
			args.block_partition = true;
		} else if (a == "--worker-count" && i + 1 < argc) {
			args.worker_count = std::atoi(argv[++i]);
			args.block_partition = true;
		} else if (a == "--manifest" && i + 1 < argc) {
			args.manifest_path = argv[++i];
			args.have_manifest = true;
		} else if (a == "--blocks" && i + 1 < argc) {
			// Parse "a_r:a_c,a_r:a_c,..." into a list of pairs.
			std::string spec = argv[++i];
			std::string token;
			for (std::size_t p = 0, start = 0; p <= spec.size(); ++p) {
				if (p == spec.size() || spec[p] == ',') {
					token = spec.substr(start, p - start);
					start = p + 1;
					if (token.empty()) continue;
					auto colon = token.find(':');
					if (colon == std::string::npos) {
						std::fprintf(stderr, "bad --blocks entry '%s' (need a_r:a_c)\n", token.c_str());
						return 2;
					}
					int ar = std::atoi(token.c_str());
					int ac = std::atoi(token.c_str() + colon + 1);
					args.selected_blocks.emplace_back(ar, ac);
				}
			}
			args.have_selected_blocks = !args.selected_blocks.empty();
		}
	}
	if (args.block_partition) {
		if (args.worker_count <= 0 || args.worker_index < 0 || args.worker_index >= args.worker_count) {
			std::fprintf(stderr, "bad --worker-index/--worker-count: %d/%d\n",
			             args.worker_index, args.worker_count);
			return 2;
		}
		if (args.shard_n != 1) {
			std::fprintf(stderr, "--shard and --worker-index/--worker-count are mutually exclusive\n");
			return 2;
		}
	}
	if (args.have_selected_blocks) {
		if (!args.have_sector) {
			std::fprintf(stderr, "--blocks requires --sector J P (exactly one sector)\n");
			return 2;
		}
		if (args.block_partition) {
			std::fprintf(stderr, "--blocks and --worker-index/--worker-count are mutually exclusive\n");
			return 2;
		}
		if (args.shard_n != 1) {
			std::fprintf(stderr, "--blocks and --shard are mutually exclusive\n");
			return 2;
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

	std::string p_hash, q_hash;
#if TICTAC_USE_NEW_CACHE_LAYER
	p_hash = tictac::cache::hash_full_raw(
	    reinterpret_cast<const unsigned char*>(fwp.p_WP_array), (fwp.Np_WP + 1) * sizeof(double));
	q_hash = tictac::cache::hash_full_raw(
	    reinterpret_cast<const unsigned char*>(fwp.q_WP_array), (fwp.Nq_WP + 1) * sizeof(double));
#endif
	const int Np_quad = std::max(1, rp.Np_per_WP_W1);
	const int Nq_quad = std::max(1, rp.Nq_per_WP_W1);

	// Enumerate sectors.
	std::vector<SectorInfo> sectors;
	for (int c = 0; c < pw.N_chn_3N; ++c) {
		SectorInfo s;
		s.chn = c;
		pw_3N_statespace sub = make_channel_view(pw, c);
		s.two_J = sub.two_J_3N_array[0];
		s.parity = sub.P_3N_array[0];
		s.Nalpha = sub.Nalpha;
		s.active = (rp.two_J_3NF_force_max < 0 || s.two_J <= rp.two_J_3NF_force_max);
		sectors.push_back(s);
	}

	// Signature for the manifest / cross-run identity.
	const auto sig = make_w1_signature(*tnf, rp, pw, fwp.Np_WP, fwp.Nq_WP,
	                                   Np_quad, Nq_quad, p_hash, q_hash);
	const std::string sig_hash = signature_hash(sig);

	// ---- helpers shared by plan/status/assemble ----
	auto build_report = [&]() -> W1ManifestReport {
		W1ManifestReport rep;
		rep.signature = sig;
		rep.signature_hash = sig_hash;
		for (const auto& s : sectors) {
			W1ManifestSector ms;
			ms.chn = s.chn; ms.two_J = s.two_J; ms.parity = s.parity;
			ms.Nalpha = s.Nalpha; ms.active = s.active;
			if (s.active) {
				pw_3N_statespace sub = make_channel_view(pw, s.chn);
				W1WorkPlan plan(sub, tnf->W1_is_exactly_hermitian());
				ms.num_blocks = plan.num_blocks();
				ms.num_evaluate = plan.num_evaluate();
				ms.num_transpose_fill = plan.num_transpose_fill();
				ms.units = plan.units();
#if TICTAC_USE_NEW_CACHE_LAYER
				std::size_t present = 0;
				for (const W1WorkUnit* u : plan.evaluate_units()) {
					auto key = make_w1_key(*tnf, rp, sub, fwp.Np_WP, fwp.Nq_WP,
					                       Np_quad, Nq_quad, u->a_r, u->a_c, p_hash, q_hash);
					W1Block blk;
					if (tictac::cache::lookup_w1(key, &blk).hit) ++present;
				}
				ms.num_evaluate_present = present;
#endif
			}
			rep.sectors.push_back(std::move(ms));
		}
		return rep;
	};

	// ============================ plan ============================
	if (args.cmd == "plan") {
		std::printf("# W1 block plan: model=%s Np=%d Nq=%d two_J_3N_max=%d two_J_3NF_force_max=%d\n",
		            tnf->name().c_str(), rp.Np_WP, rp.Nq_WP, rp.two_J_3N_max, rp.two_J_3NF_force_max);
		std::printf("# signature_hash=%s\n", sig_hash.c_str());
		std::printf("# %-6s %-6s %-8s %-10s %-10s %-12s %-12s %-14s\n",
		            "chn", "two_J", "parity", "Nalpha", "num_blocks", "num_eval", "num_transp", "eval_present");
		std::size_t tot_blocks = 0, tot_eval = 0, tot_transp = 0, tot_present = 0;
		for (const auto& s : sectors) {
			if (!s.active) {
				std::printf("  %-6d %-6d %-8d %-10d %-10s %-12s %-12s %-14s\n",
				            s.chn, s.two_J, s.parity, s.Nalpha, "(2NF)", "-", "-", "-");
				continue;
			}
			pw_3N_statespace sub = make_channel_view(pw, s.chn);
			W1WorkPlan plan(sub, tnf->W1_is_exactly_hermitian());
#if TICTAC_USE_NEW_CACHE_LAYER
			std::size_t present = 0;
			for (const W1WorkUnit* u : plan.evaluate_units()) {
				auto key = make_w1_key(*tnf, rp, sub, fwp.Np_WP, fwp.Nq_WP,
				                       Np_quad, Nq_quad, u->a_r, u->a_c, p_hash, q_hash);
				W1Block blk;
				if (tictac::cache::lookup_w1(key, &blk).hit) ++present;
			}
#else
			std::size_t present = 0;
#endif
			tot_blocks += plan.num_blocks();
			tot_eval += plan.num_evaluate();
			tot_transp += plan.num_transpose_fill();
			tot_present += present;
			std::printf("  %-6d %-6d %-8d %-10d %-10zu %-12zu %-12zu %-14zu\n",
			            s.chn, s.two_J, s.parity, s.Nalpha, plan.num_blocks(),
			            plan.num_evaluate(), plan.num_transpose_fill(), present);
		}
		std::printf("# active: blocks=%zu evaluate=%zu transpose_fill=%zu present=%zu\n",
		            tot_blocks, tot_eval, tot_transp, tot_present);
		if (args.have_manifest) {
			W1ManifestReport rep = build_report();
			if (write_manifest_json(args.manifest_path, rep)) {
				std::printf("# manifest written: %s\n", args.manifest_path.c_str());
			} else {
				std::fprintf(stderr, "w1_worker: failed to write manifest %s\n", args.manifest_path.c_str());
			}
		}
#if TICTAC_USE_NEW_CACHE_LAYER
		tictac::cache::shutdown();
#endif
		return 0;
	}

	// ============================ status ============================
	if (args.cmd == "status") {
		W1ManifestReport rep = build_report();
		std::printf("# W1 status: model=%s signature_hash=%s\n",
		            tnf->name().c_str(), sig_hash.c_str());
		std::printf("# %-6s %-6s %-8s %-10s %-12s %-12s %-14s %s\n",
		            "chn", "two_J", "parity", "Nalpha", "num_eval", "eval_present", "complete", "");
		std::size_t tot_eval = 0, tot_present = 0;
		bool all_complete = true;
		for (const auto& s : rep.sectors) {
			if (!s.active) continue;
			const bool complete = (s.num_evaluate_present == s.num_evaluate);
			if (!complete) all_complete = false;
			tot_eval += s.num_evaluate;
			tot_present += s.num_evaluate_present;
			std::printf("  %-6d %-6d %-8d %-10d %-12zu %-12zu %-14s\n",
			            s.chn, s.two_J, s.parity, s.Nalpha, s.num_evaluate,
			            s.num_evaluate_present, complete ? "yes" : "NO");
		}
		std::printf("# evaluate: %zu/%zu present%s\n", tot_present, tot_eval,
		            all_complete ? " -- COMPLETE" : " -- INCOMPLETE");
		if (args.have_manifest) {
			// Optionally persist a fresh status snapshot next to the manifest.
			const std::string status_path = args.manifest_path + ".status.json";
			write_manifest_json(status_path, rep);
		}
#if TICTAC_USE_NEW_CACHE_LAYER
		tictac::cache::shutdown();
#endif
		return all_complete ? 0 : 1;   // exit code: 0 complete, 1 incomplete
	}

	// ============================ assemble ============================
	if (args.cmd == "assemble") {
		std::string hash, missing;
		bool ok = tictac::interactions::W1Assembler::fingerprint(
		    *tnf, rp, pw, fwp.Np_WP, fwp.Nq_WP, Np_quad, Nq_quad, p_hash, q_hash,
		    hash, &missing);
		if (!ok) {
			std::fprintf(stderr, "[w1_worker] assemble: INCOMPLETE, first missing %s\n", missing.c_str());
#if TICTAC_USE_NEW_CACHE_LAYER
			tictac::cache::shutdown();
#endif
			return 1;
		}
		std::printf("[w1_worker] assemble: COMPLETE signature_hash=%s\n", sig_hash.c_str());
		std::printf("[w1_worker] assemble: W1_fingerprint=%s\n", hash.c_str());
		// Verify Hermitian transpose contract on every off-diagonal pair: the
		// transpose_fill block's payload must equal the exact transpose of its
		// evaluated conjugate.  This is the in-cache reverse-block check.
#if TICTAC_USE_NEW_CACHE_LAYER
		std::size_t checked_pairs = 0, bad_pairs = 0;
		if (tnf->W1_is_exactly_hermitian()) {
			for (const auto& s : sectors) {
				if (!sector_selected(s, args) && args.have_sector) continue;
				if (!s.active) continue;
				pw_3N_statespace sub = make_channel_view(pw, s.chn);
				W1WorkPlan plan(sub, true);
				for (const auto& u : plan.units()) {
					if (u.role != W1UnitRole::transpose_fill) continue;
					auto kf = make_w1_key(*tnf, rp, sub, fwp.Np_WP, fwp.Nq_WP, Np_quad, Nq_quad,
					                     u.conj_a_r, u.conj_a_c, p_hash, q_hash);
					auto kt = make_w1_key(*tnf, rp, sub, fwp.Np_WP, fwp.Nq_WP, Np_quad, Nq_quad,
					                     u.a_r, u.a_c, p_hash, q_hash);
					W1Block bf, bt;
					if (!tictac::cache::lookup_w1(kf, &bf).hit) continue;
					if (!tictac::cache::lookup_w1(kt, &bt).hit) {
						// transpose_fill block absent is fine (solver fills on demand);
						// the Hermitian check only fires when both are present.
						continue;
					}
					++checked_pairs;
					const std::size_t Np = bf.Np, Nq = bf.Nq;
					bool ok_pair = true;
					for (std::size_t iqr = 0; iqr < Nq && ok_pair; ++iqr)
						for (std::size_t iqc = 0; iqc < Nq && ok_pair; ++iqc)
							for (std::size_t ipr = 0; ipr < Np && ok_pair; ++ipr)
								for (std::size_t ipc = 0; ipc < Np && ok_pair; ++ipc) {
									const std::size_t fwd = (((iqr*Nq+iqc)*Np+ipr)*Np+ipc);
									const std::size_t rev = (((iqc*Nq+iqr)*Np+ipc)*Np+ipr);
									if (bf.data[fwd] != bt.data[rev]) ok_pair = false;
								}
					if (!ok_pair) ++bad_pairs;
				}
			}
		}
		std::printf("[w1_worker] assemble: Hermitian transpose pairs checked=%zu bad=%zu\n",
		            checked_pairs, bad_pairs);
		if (bad_pairs > 0) {
			std::fprintf(stderr, "[w1_worker] assemble: FAILED Hermitian contract\n");
			tictac::cache::shutdown();
			return 1;
		}
#endif
		if (args.have_manifest) {
			W1ManifestReport rep = build_report();
			write_manifest_json(args.manifest_path, rep);
		}
#if TICTAC_USE_NEW_CACHE_LAYER
		tictac::cache::shutdown();
#endif
		return 0;
	}

	// ============================ build: --blocks (selected pilot) ============
	if (args.cmd == "build" && args.have_selected_blocks) {
		// Build ONLY the explicitly-listed evaluate-units in the single selected
		// sector. Each block goes through W1BlockExecutor::compute_block, which
		// calls the SAME integrate_w1_channel_blocks as the monolithic and
		// block-level paths -- bitwise-identical by construction.  Per-block
		// wall time is printed as a JSON line for machine-readable pilot timing.
		std::size_t built = 0, hits = 0, skipped_not_eval = 0, skipped_not_found = 0;
		const int omp_threads = omp_get_max_threads();
		for (const auto& s : sectors) {
			if (!s.active) continue;
			if (!(s.two_J == args.sel_two_J && s.parity == args.sel_P)) continue;
			pw_3N_statespace sub = make_channel_view(pw, s.chn);
			W1WorkPlan plan(sub, tnf->W1_is_exactly_hermitian());
			const std::size_t per_block = static_cast<std::size_t>(fwp.Nq_WP) * fwp.Nq_WP
			                            * fwp.Np_WP * fwp.Np_WP;
			for (const auto& [ar, ac] : args.selected_blocks) {
				const W1WorkUnit* u = plan.unit_for(ar, ac);
				if (!u) { ++skipped_not_found; continue; }
				if (u->role != W1UnitRole::evaluate) {
					++skipped_not_eval;
					std::printf("{\"block\":\"%d:%d\",\"status\":\"transpose_fill\",\"a_r\":%d,\"a_c\":%d}\n",
					            ar, ac, ar, ac);
					continue;
				}
				auto key = make_w1_key(*tnf, rp, sub, fwp.Np_WP, fwp.Nq_WP,
				                       Np_quad, Nq_quad, ar, ac, p_hash, q_hash);
				#if TICTAC_USE_NEW_CACHE_LAYER
				W1Block blk;
				if (tictac::cache::lookup_w1(key, &blk).hit) {
					++hits;
					std::printf("{\"block\":\"%d:%d\",\"status\":\"cache_hit\",\"a_r\":%d,\"a_c\":%d,"
					            "\"payload_bytes\":%zu,\"Np\":%d,\"Nq\":%d}\n",
					            ar, ac, ar, ac, blk.data.size()*sizeof(double),
					            blk.Np, blk.Nq);
					continue;
				}
				#endif
				std::vector<double> buf(per_block, 0.0);
				const double t0 = omp_get_wtime();
				W1BlockExecutor::compute_block(*tnf, fwp.p_WP_array, fwp.Np_WP,
				                               fwp.q_WP_array, fwp.Nq_WP, sub, rp,
				                               ar, ac, buf.data());
				const double t1 = omp_get_wtime();
				#if TICTAC_USE_NEW_CACHE_LAYER
				W1Block out{};
				out.Nq = fwp.Nq_WP; out.Np = fwp.Np_WP;
				out.a_r = ar; out.a_c = ac;
				out.data = std::move(buf);
				tictac::cache::store_w1(key, out);
				#endif
				++built;
				std::printf("{\"block\":\"%d:%d\",\"status\":\"built\",\"a_r\":%d,\"a_c\":%d,"
				            "\"wall_seconds\":%.6f,\"omp_threads\":%d,"
				            "\"payload_bytes\":%zu,\"Np\":%d,\"Nq\":%d,"
				            "\"two_J\":%d,\"parity\":%d,\"Nalpha\":%d,"
				            "\"Np_WP\":%d,\"Nq_WP\":%d,\"Np_per_WP_W1\":%d,\"Nq_per_WP_W1\":%d,"
				            "\"Nangle_3NF\":%d,\"signature_hash\":\"%s\"}\n",
				            ar, ac, ar, ac, t1 - t0, omp_threads,
				            per_block * sizeof(double), fwp.Np_WP, fwp.Nq_WP,
				            s.two_J, s.parity, s.Nalpha,
				            fwp.Np_WP, fwp.Nq_WP, Np_quad, Nq_quad,
				            rp.Nangle_3NF, sig_hash.c_str());
			}
		}
		std::fprintf(stderr, "[w1_worker] build --blocks: built=%zu hits=%zu "
		             "skipped_not_eval=%zu skipped_not_found=%zu\n",
		             built, hits, skipped_not_eval, skipped_not_found);
		#if TICTAC_USE_NEW_CACHE_LAYER
		tictac::cache::shutdown();
		#endif
		return 0;
	}

	// ============================ build: --worker-index (block-level) ========
	if (args.cmd == "build" && args.block_partition) {
		// Block-level: each worker integrates its evaluate-units.  To recover
			// the cross-channel orbital-cache reuse that the monolithic
			// W1_PW_cache::build enjoys, a worker batches ALL its owned missing
			// channels in a sector into ONE integrate_w1_channel_blocks call (the
			// same shared per-cell integration).  This is bitwise-identical to
			// per-block calls (channels are independent in the cell loop) but
			// avoids re-populating the momentum-dependent orbital cache per block.
			std::size_t total_owned = 0, total_built = 0, total_hits = 0;
			std::size_t global_eval = 0;   // running index across selected sectors
			for (const auto& s : sectors) {
				if (!s.active) continue;
				if (args.have_sector && !(s.two_J == args.sel_two_J && s.parity == args.sel_P)) continue;
				pw_3N_statespace sub = make_channel_view(pw, s.chn);
				W1WorkPlan plan(sub, tnf->W1_is_exactly_hermitian());
				const auto evals = plan.evaluate_units();
				// 1) partition + cache probe: collect this worker's missing channels.
				std::vector<std::pair<int, int>> miss_channels;
				std::vector<W1Key> miss_keys;
				std::size_t sector_owned = 0, sector_hits = 0;
				for (std::size_t li = 0; li < evals.size(); ++li, ++global_eval) {
					if ((global_eval % args.worker_count) != args.worker_index) continue;
					++sector_owned;
					const W1WorkUnit* u = evals[li];
#if TICTAC_USE_NEW_CACHE_LAYER
					auto key = make_w1_key(*tnf, rp, sub, fwp.Np_WP, fwp.Nq_WP,
					                       Np_quad, Nq_quad, u->a_r, u->a_c, p_hash, q_hash);
					W1Block blk;
					if (tictac::cache::lookup_w1(key, &blk).hit) { ++sector_hits; continue; }
					miss_channels.emplace_back(u->a_r, u->a_c);
					miss_keys.push_back(key);
#else
					miss_channels.emplace_back(u->a_r, u->a_c);
#endif
				}
				// 2) batched integration of all missing owned channels in one call.
				std::size_t sector_built = miss_channels.size();
				if (!miss_channels.empty()) {
					const std::size_t per_block = static_cast<std::size_t>(fwp.Nq_WP) * fwp.Nq_WP
					                            * fwp.Np_WP * fwp.Np_WP;
					std::vector<std::vector<double>> bufs(miss_channels.size(),
					                                     std::vector<double>(per_block));
					std::vector<double*> ptrs;
					ptrs.reserve(bufs.size());
					for (auto& b : bufs) ptrs.push_back(b.data());
					integrate_w1_channel_blocks(*tnf, fwp.p_WP_array, fwp.Np_WP,
					                            fwp.q_WP_array, fwp.Nq_WP, sub, rp,
					                            miss_channels, ptrs);
#if TICTAC_USE_NEW_CACHE_LAYER
					for (std::size_t k = 0; k < miss_channels.size(); ++k) {
						W1Block out{};
						out.Nq = fwp.Nq_WP; out.Np = fwp.Np_WP;
						out.a_r = miss_channels[k].first; out.a_c = miss_channels[k].second;
						out.data = std::move(bufs[k]);
						tictac::cache::store_w1(miss_keys[k], out);
					}
#endif
				}
				if (sector_owned > 0) {
					std::fprintf(stderr,
					    "[3NF W1] evaluating %zu of %zu missing channel blocks: "
					    "sector J=%d/2 P=%d worker=%d/%d built=%zu hits=%zu\n",
					    sector_built, sector_owned, s.two_J, s.parity,
					    args.worker_index, args.worker_count, sector_built, sector_hits);
				}
				total_owned += sector_owned; total_built += sector_built; total_hits += sector_hits;
			}
			std::printf("[w1_worker] build: worker=%d/%d owned=%zu built=%zu cache_hits=%zu\n",
			            args.worker_index, args.worker_count, total_owned, total_built, total_hits);
#if TICTAC_USE_NEW_CACHE_LAYER
			tictac::cache::shutdown();
#endif
			return 0;
		}

	// ============================ build: --shard (sector-level) ============
	if (args.cmd == "build") {
		// Sector-level (legacy --shard): build each owned sector monolithically
		// via W1_PW_cache::build, which probes the cache per block (hits skipped),
		// evaluates the Hermitian triangle, transpose-fills reverse blocks, and
		// stores all of them atomically.
		int built_sectors = 0;
		for (const auto& s : sectors) {
			if (!sector_selected(s, args)) continue;
			std::printf("[w1_worker] building sector chn=%d J=%d/2 P=%d (Nalpha^2=%d blocks) ...\n",
			            s.chn, s.two_J, s.parity, s.Nalpha * s.Nalpha);
			pw_3N_statespace sub = make_channel_view(pw, s.chn);
			W1_PW_cache cache;
			cache.build(*tnf, fwp.p_WP_array, fwp.Np_WP,
			            fwp.q_WP_array, fwp.Nq_WP, sub, rp);
			++built_sectors;
		}
		std::printf("[w1_worker] built %d sector(s).\n", built_sectors);
#if TICTAC_USE_NEW_CACHE_LAYER
		tictac::cache::shutdown();
#endif
		return 0;
	}

	// ============================ verify ============================
	if (args.cmd == "verify") {
		int bad = 0; long checked = 0;
		for (const auto& s : sectors) {
			if (!sector_selected(s, args)) continue;
			pw_3N_statespace sub = make_channel_view(pw, s.chn);
			W1WorkPlan plan(sub, tnf->W1_is_exactly_hermitian());
			// Verify the EVALUATE units (transpose_fill units are produced by the
			// solver on demand, so their absence is not a defect). A complete
			// campaign is one where every evaluate-unit shard is a valid hit.
			for (const W1WorkUnit* u : plan.evaluate_units()) {
#if TICTAC_USE_NEW_CACHE_LAYER
				auto key = make_w1_key(*tnf, rp, sub, fwp.Np_WP, fwp.Nq_WP,
				                       Np_quad, Nq_quad, u->a_r, u->a_c, p_hash, q_hash);
				W1Block blk;
				auto res = tictac::cache::lookup_w1(key, &blk);
				if (!res.hit || blk.Np != fwp.Np_WP || blk.Nq != fwp.Nq_WP
				    || blk.data.size() != static_cast<std::size_t>(fwp.Nq_WP)*fwp.Nq_WP*fwp.Np_WP*fwp.Np_WP) {
					++bad;
				}
				++checked;
#endif
			}
		}
		std::printf("[w1_worker] verify: %ld evaluate-blocks checked, %d missing/invalid\n", checked, bad);
#if TICTAC_USE_NEW_CACHE_LAYER
		tictac::cache::shutdown();
#endif
		return bad == 0 ? 0 : 1;
	}

	std::fprintf(stderr, "unknown command '%s' (use plan|build|status|assemble|verify)\n", args.cmd.c_str());
	return 2;
}
