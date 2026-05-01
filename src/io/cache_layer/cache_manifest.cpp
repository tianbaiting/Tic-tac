#include "cache_manifest.h"
#include "third_party/json.hpp"

#include <chrono>
#include <ctime>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>

namespace tictac::cache {

using json = nlohmann::json;

std::string now_utc_iso8601() {
    auto now = std::chrono::system_clock::now();
    std::time_t t = std::chrono::system_clock::to_time_t(now);
    std::tm tm_utc;
    gmtime_r(&t, &tm_utc);
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", &tm_utc);
    return std::string(buf);
}

static std::string manifest_path(const std::string& cache_root) {
    return cache_root + "/manifest.json";
}

Manifest Manifest::load(const std::string& cache_root) {
    Manifest m;
    std::ifstream f(manifest_path(cache_root));
    if (!f.good()) return m;
    json j;
    try { f >> j; } catch (...) { return m; }
    if (!j.contains("entries") || !j["entries"].is_array()) return m;
    for (const auto& e : j["entries"]) {
        ManifestEntry me{};
        me.kind             = e.value("kind", "");
        me.filename         = e.value("filename", "");
        me.key_hash_full    = e.value("key_hash_full", "");
        me.key_json         = e.value("key_json", "");
        me.size_bytes       = e.value("size_bytes", (int64_t)0);
        me.schema_version   = e.value("schema_version", 0);
        me.created_utc      = e.value("created_utc", "");
        me.last_used_utc    = e.value("last_used_utc", "");
        me.legacy_imported  = e.value("legacy_imported", false);
        me.writer_version   = e.value("writer_version", "");
        if (!me.key_hash_full.empty()) m.m_entries.push_back(std::move(me));
    }
    return m;
}

void Manifest::save(const std::string& cache_root) const {
    json j;
    j["manifest_version"] = MANIFEST_VERSION;
    j["entries"] = json::array();
    for (const auto& e : m_entries) {
        json je;
        je["kind"]            = e.kind;
        je["filename"]        = e.filename;
        je["key_hash_full"]   = e.key_hash_full;
        try { je["key"] = json::parse(e.key_json); }
        catch (...) { je["key"] = e.key_json; }
        je["key_json"]        = e.key_json;
        je["size_bytes"]      = e.size_bytes;
        je["schema_version"]  = e.schema_version;
        je["created_utc"]     = e.created_utc;
        je["last_used_utc"]   = e.last_used_utc;
        je["legacy_imported"] = e.legacy_imported;
        je["writer_version"]  = e.writer_version;
        j["entries"].push_back(std::move(je));
    }
    std::string tmp = manifest_path(cache_root) + ".tmp";
    {
        std::ofstream out(tmp);
        out << j.dump(2);
    }
    if (std::rename(tmp.c_str(), manifest_path(cache_root).c_str()) != 0) {
        std::perror("manifest rename");
    }
}

const ManifestEntry* Manifest::find(const std::string& key_hash_full) const {
    for (const auto& e : m_entries)
        if (e.key_hash_full == key_hash_full) return &e;
    return nullptr;
}

bool Manifest::upsert(const ManifestEntry& entry) {
    for (auto& e : m_entries) {
        if (e.key_hash_full == entry.key_hash_full) {
            e = entry;
            return false;
        }
    }
    m_entries.push_back(entry);
    return true;
}

void Manifest::touch(const std::string& key_hash_full) {
    for (auto& e : m_entries) {
        if (e.key_hash_full == key_hash_full) {
            e.last_used_utc = now_utc_iso8601();
            return;
        }
    }
}

bool Manifest::erase(const std::string& key_hash_full) {
    for (auto it = m_entries.begin(); it != m_entries.end(); ++it) {
        if (it->key_hash_full == key_hash_full) {
            m_entries.erase(it);
            return true;
        }
    }
    return false;
}

}  // namespace tictac::cache
