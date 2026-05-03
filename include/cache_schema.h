#pragma once

namespace tictac::cache {

// Bump when on-disk format or physics definition changes.
// Version is part of the hash key, so old caches are preserved (not overwritten)
// when this is bumped.
constexpr int P123_SCHEMA_VERSION = 1;
constexpr int W1_SCHEMA_VERSION   = 2;

}  // namespace tictac::cache
