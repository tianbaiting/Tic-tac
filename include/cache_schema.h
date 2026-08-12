#pragma once

namespace tictac::cache {

// Bump when on-disk format or physics definition changes.
// Version is part of the hash key, so old caches are preserved (not overwritten)
// when this is bumped.
constexpr int P123_SCHEMA_VERSION = 1;
// Bumped 3 → 4 on 2026-06-21 (3NF audit B5): added c_1, c_3, c_4 and grid
// hashes to W1Key so caches cannot be wrongly reused across different
// Hamiltonians / momentum grids. Old v3 caches are read as misses.
// Bumped 4 -> 5 on 2026-08-13: cE now uses the exact Epelbaum A-4 angular
// factor with the two-Jacobi-coordinate Fourier normalization.  Old v4 W1
// blocks have a cE coefficient larger by pi/2 and must remain cache misses.
// Bumped 5 -> 6 on 2026-08-13: complete-reference angular order and all
// dimensionful/dimensionless chiral kernel constants are now part of W1Key.
constexpr int W1_SCHEMA_VERSION   = 6;

}  // namespace tictac::cache
