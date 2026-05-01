#include "cache_layer.h"
#include "cache_manifest.h"

#include <sys/stat.h>
#include <sys/types.h>
#include <mutex>
#include <iostream>
#include <cerrno>
#include <cstring>
#include <fstream>

namespace tictac::cache {

static std::string         g_root;
static Manifest            g_manifest;
static std::mutex          g_mu;
static std::string         g_writer_version = "tictac-unversioned";

static void mkdir_p(const std::string& path) {
    std::string cur;
    for (size_t i = 0; i < path.size(); ++i) {
        if (path[i] == '/') {
            if (!cur.empty()) {
                ::mkdir(cur.c_str(), 0755);
            }
        }
        cur.push_back(path[i]);
    }
    if (!cur.empty()) ::mkdir(cur.c_str(), 0755);
}

static int64_t file_size(const std::string& path) {
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return 0;
    return (int64_t)st.st_size;
}

void initialize(const std::string& cache_root_in) {
    std::lock_guard<std::mutex> g(g_mu);
    g_root = cache_root_in;
    mkdir_p(g_root);
    mkdir_p(g_root + "/p123");
    mkdir_p(g_root + "/w1");
    g_manifest = Manifest::load(g_root);
}

void shutdown() {
    std::lock_guard<std::mutex> g(g_mu);
    if (!g_root.empty()) g_manifest.save(g_root);
}

const std::string& cache_root() { return g_root; }

static std::string p123_filename(const P123Key& k) {
    return "p123/" + filename_prefix(k) + "__" + hash_short(hash_full(k)) + ".h5";
}
static std::string w1_filename(const W1Key& k) {
    return "w1/" + filename_prefix(k) + "__" + hash_short(hash_full(k)) + ".h5";
}

LookupResult lookup_p123(const P123Key& key, P123Arrays* out) {
    LookupResult r{}; r.hit = false;
    if (g_root.empty()) { r.miss_reason = "not_initialized"; return r; }
    auto rel = p123_filename(key);
    auto abs = g_root + "/" + rel;
    r.source_path = abs;
    bool ok = read_p123_h5(abs, key, out, &r.miss_reason);
    if (ok) {
        r.hit = true; r.miss_reason.clear();
        std::lock_guard<std::mutex> g(g_mu);
        g_manifest.touch(hash_full(key));
    }
    return r;
}

void store_p123(const P123Key& key, const P123Arrays& in) {
    if (g_root.empty()) return;
    auto rel = p123_filename(key);
    auto abs = g_root + "/" + rel;
    write_p123_h5(abs, key, in, g_writer_version);

    ManifestEntry e{};
    e.kind = "p123"; e.filename = rel;
    e.key_hash_full = hash_full(key);
    e.key_json = canonical_json(key);
    e.size_bytes = file_size(abs);
    e.schema_version = key.schema_version;
    e.created_utc = now_utc_iso8601();
    e.last_used_utc = e.created_utc;
    e.legacy_imported = false;
    e.writer_version = g_writer_version;
    {
        std::lock_guard<std::mutex> g(g_mu);
        g_manifest.upsert(e);
        g_manifest.save(g_root);   // flush eagerly so concurrent readers see it
    }
}

LookupResult lookup_w1(const W1Key& key, W1Block* out) {
    LookupResult r{}; r.hit = false;
    if (g_root.empty()) { r.miss_reason = "not_initialized"; return r; }
    auto rel = w1_filename(key);
    auto abs = g_root + "/" + rel;
    r.source_path = abs;
    bool ok = read_w1_h5(abs, key, out, &r.miss_reason);
    if (ok) {
        r.hit = true; r.miss_reason.clear();
        std::lock_guard<std::mutex> g(g_mu);
        g_manifest.touch(hash_full(key));
    }
    return r;
}

void store_w1(const W1Key& key, const W1Block& in) {
    if (g_root.empty()) return;
    auto rel = w1_filename(key);
    auto abs = g_root + "/" + rel;
    write_w1_h5(abs, key, in, g_writer_version);

    ManifestEntry e{};
    e.kind = "w1"; e.filename = rel;
    e.key_hash_full = hash_full(key);
    e.key_json = canonical_json(key);
    e.size_bytes = file_size(abs);
    e.schema_version = key.schema_version;
    e.created_utc = now_utc_iso8601();
    e.last_used_utc = e.created_utc;
    e.legacy_imported = false;
    e.writer_version = g_writer_version;
    {
        std::lock_guard<std::mutex> g(g_mu);
        g_manifest.upsert(e);
        g_manifest.save(g_root);
    }
}

CacheStats summary() {
    CacheStats s{};
    std::lock_guard<std::mutex> g(g_mu);
    for (const auto& e : g_manifest.entries()) {
        if (e.kind == "p123") ++s.n_p123;
        if (e.kind == "w1")   ++s.n_w1;
        s.total_bytes += e.size_bytes;
    }
    return s;
}

}  // namespace tictac::cache
