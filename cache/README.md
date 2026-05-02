# Tic-tac cache directory

This directory holds hash-keyed HDF5 files produced by the Tic-tac cache layer
(Task 8/9).  It is **gitignored** — only this README is tracked.

## Layout

```
cache/
  p123/          # P123 partial-wave arrays, keyed by P123Key hash
  w1/            # W1 per-block arrays, keyed by W1BlockKey hash
  manifest.json  # Human-readable index of all cached entries
  README.md      # This file (force-committed despite gitignore)
```

## How the cache works

On every solver run, `cache_layer_initialize(root)` opens (or creates) the
manifest.  For each P123 or W1 computation the solver would otherwise do from
scratch, it first calls `p123_lookup` / `w1_lookup`:

1. Serialize the key struct → canonical JSON → SHA-256 (first 8 hex chars form
   the filename suffix).
2. Look up `cache/p123/<label>__<hash>.h5` (or `cache/w1/…`).
3. **Hit**: deserialize from HDF5, skip computation.
4. **Miss**: run computation, write HDF5, update manifest.

On shutdown, `cache_layer_shutdown()` flushes the manifest.

## Cold vs warm runs

| Run | P123 files present | Typical speedup |
|-----|-------------------|-----------------|
| Cold | none | baseline (1×) |
| Warm | all hits | 3–10× (disk I/O only) |

Speedup depends on Np, Nq, J2max, and disk speed.  Small test fixtures
(Np = Nq = 5) show a smaller ratio because Python/process overhead dominates.

## Disabling the cache

```cmake
cmake -DTICTAC_USE_NEW_CACHE_LAYER=OFF ..
```

When `OFF`, the cache calls compile to no-ops and no HDF5 files are touched.

## Migrating legacy P123 files

If you have pre-existing `P123_sparse_JP_*.h5` files in `CPP/Output/`, the
migration script can copy them into the new layout without recomputing:

```bash
# Dry-run first (default, no files written):
python3 scripts/migrate_legacy_cache.py \
    --legacy-roots CPP/Output \
    --target-root cache

# Copy with verification (originals preserved):
python3 scripts/migrate_legacy_cache.py \
    --legacy-roots CPP/Output \
    --target-root cache \
    --apply --copy
```

Use `--apply --copy` — never `--apply` alone — to preserve original files.
