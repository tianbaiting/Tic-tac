#pragma once
#include <string>
#include <vector>
#include <cstdint>

namespace tictac::cache {

struct ManifestEntry {
    std::string kind;             // "p123" | "w1"
    std::string filename;         // relative to cache root, e.g. "p123/Np30_..._a3f9c2.h5"
    std::string key_hash_full;    // 64-char hex
    std::string key_json;         // canonical_json output
    int64_t     size_bytes;
    int         schema_version;
    std::string created_utc;      // ISO8601, e.g. "2026-05-01T12:34:56Z"
    std::string last_used_utc;    // ISO8601
    bool        legacy_imported;
    std::string writer_version;   // git describe or "legacy_imported_v1"
};

class Manifest {
public:
    static constexpr int MANIFEST_VERSION = 1;

    static Manifest load(const std::string& cache_root);
    void save(const std::string& cache_root) const;
    const ManifestEntry* find(const std::string& key_hash_full) const;
    bool upsert(const ManifestEntry& entry);  // true if inserted, false if updated existing
    void touch(const std::string& key_hash_full);
    bool erase(const std::string& key_hash_full);

    const std::vector<ManifestEntry>& entries() const { return m_entries; }

private:
    std::vector<ManifestEntry> m_entries;
};

std::string now_utc_iso8601();

}  // namespace tictac::cache
