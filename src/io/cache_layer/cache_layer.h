#pragma once
#include "cache_keys.h"
#include "cache_io_p123.h"
#include "cache_io_w1.h"

#include <functional>
#include <string>

namespace tictac::cache {

struct LookupResult {
    bool        hit;
    std::string source_path;
    std::string miss_reason;   // "" on hit; otherwise tag from cache_io_*
};

// Lifecycle: must be called once at program start.
void initialize(const std::string& cache_root);
void shutdown();
const std::string& cache_root();

// P123
LookupResult lookup_p123(const P123Key& key, P123Arrays* out);
void         store_p123 (const P123Key& key, const P123Arrays& in);

// W1 (per JP block)
LookupResult lookup_w1  (const W1Key& key, W1Block* out);
void         store_w1   (const W1Key& key, const W1Block& in);

// Stats
struct CacheStats { size_t n_p123; size_t n_w1; int64_t total_bytes; };
CacheStats summary();

}  // namespace tictac::cache
