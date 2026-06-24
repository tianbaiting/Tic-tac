#pragma once

namespace tictac::cache {

// Bump when on-disk format or physics definition changes.
// Version is part of the hash key, so old caches are preserved (not overwritten)
// when this is bumped.
constexpr int P123_SCHEMA_VERSION = 1;
// Bumped 3 → 4 on 2026-06-21 (3NF audit B5): added c_1, c_3, c_4 and grid
// hashes to W1Key so caches cannot be wrongly reused across different
// Hamiltonians / momentum grids. Old v3 caches are read as misses.
constexpr int W1_SCHEMA_VERSION   = 4;

}  // namespace tictac::cache
