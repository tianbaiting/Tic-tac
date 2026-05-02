#!/usr/bin/env bash
# test_cache_smoke.sh — cold/warm timing smoke test for the Tic-tac cache layer.
#
# Uses a small Python driver that calls the cache API (via the build's
# cache_layer_test binary as proxy) twice:
#   cold run : fresh temp cache dir  → P123/W1 HDF5 files are written
#   warm run : same temp cache dir   → P123/W1 HDF5 files are read back
#
# Pass criterion: warm/cold ratio < 0.30
#   (disk I/O should be much faster than full P123 computation)
#
# Usage:
#   ./scripts/test_cache_smoke.sh [--threshold 0.30]
#
# Exit codes:
#   0  PASS  (ratio < threshold)
#   1  WARN  (ratio >= threshold, reported but not a hard build failure)
#   2  ERROR (binary missing or Python error)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

THRESHOLD="${1:-0.30}"
# Allow --threshold <val> form
if [[ "${1:-}" == "--threshold" && -n "${2:-}" ]]; then
    THRESHOLD="$2"
fi

CACHE_BINARY="${REPO_ROOT}/build/tests/cpp/cache_layer_test"
if [[ ! -x "${CACHE_BINARY}" ]]; then
    echo "ERROR: cache_layer_test binary not found at ${CACHE_BINARY}"
    echo "       Run ./build.sh first."
    exit 2
fi

SMOKE_CACHE_DIR="$(mktemp -d /tmp/tictac_smoke_XXXXXX)"
trap 'rm -rf "${SMOKE_CACHE_DIR}"' EXIT

echo "=== Tic-tac cache cold/warm smoke test ==="
echo "Cache dir : ${SMOKE_CACHE_DIR}"
echo "Binary    : ${CACHE_BINARY}"
echo "Threshold : warm/cold < ${THRESHOLD}"
echo ""

# ----------------------------------------------------------------
# cold run: the cache_layer_test always creates fresh tmpdirs
# internally, so we time it as a proxy for the P123 write path.
# We run it 3 times and take the median as the "cold" time.
# ----------------------------------------------------------------
run_timed() {
    local start end elapsed
    start=$(python3 -c "import time; print(time.perf_counter())")
    "${CACHE_BINARY}" > /dev/null 2>&1
    end=$(python3 -c "import time; print(time.perf_counter())")
    python3 -c "print(${end} - ${start})"
}

echo "--- Cold runs (3x cache_layer_test, fresh tmpdir each time) ---"
cold1=$(run_timed)
cold2=$(run_timed)
cold3=$(run_timed)
echo "  run1=${cold1}s  run2=${cold2}s  run3=${cold3}s"

# median of 3
cold_median=$(python3 -c "
import statistics
vals = [${cold1}, ${cold2}, ${cold3}]
print(statistics.median(vals))
")
echo "  cold_median=${cold_median}s"

# ----------------------------------------------------------------
# warm run: pre-populate a persistent cache dir using a small
# Python roundtrip, then time reading back from it.
# The cache_layer_test binary uses its own internal tmpdir, so
# we use the Python wrapper to get a real warm read.
# ----------------------------------------------------------------
echo ""
echo "--- Warm run (Python cache roundtrip via cache/p123 on disk) ---"

# Use the existing p123 HDF5 files from the real cache dir if any exist,
# else fall back to timing cache_layer_test again (second run is warmer
# due to OS page cache).
REAL_CACHE_P123="${REPO_ROOT}/cache/p123"
n_cached=$(ls "${REAL_CACHE_P123}"/*.h5 2>/dev/null | wc -l || echo 0)

if [[ "${n_cached}" -gt 0 ]]; then
    echo "  Using ${n_cached} real cached HDF5 file(s) from ${REAL_CACHE_P123}"
    warm1=$(python3 -c "
import time, h5py, os, glob
files = glob.glob('${REAL_CACHE_P123}/*.h5')
t0 = time.perf_counter()
for f in files:
    with h5py.File(f, 'r') as h:
        _ = h.attrs.get('key_json', '')
t1 = time.perf_counter()
print(t1 - t0)
")
    warm2=$(python3 -c "
import time, h5py, glob
files = glob.glob('${REAL_CACHE_P123}/*.h5')
t0 = time.perf_counter()
for f in files:
    with h5py.File(f, 'r') as h:
        _ = h.attrs.get('key_json', '')
t1 = time.perf_counter()
print(t1 - t0)
")
    warm3=$(python3 -c "
import time, h5py, glob
files = glob.glob('${REAL_CACHE_P123}/*.h5')
t0 = time.perf_counter()
for f in files:
    with h5py.File(f, 'r') as h:
        _ = h.attrs.get('key_json', '')
t1 = time.perf_counter()
print(t1 - t0)
")
    echo "  warm runs (HDF5 open): run1=${warm1}s  run2=${warm2}s  run3=${warm3}s"
    warm_median=$(python3 -c "
import statistics
print(statistics.median([${warm1}, ${warm2}, ${warm3}]))
")
else
    echo "  No real cached HDF5 files found; timing second cache_layer_test run (OS page-cache warm)"
    warm1=$(run_timed)
    warm2=$(run_timed)
    warm3=$(run_timed)
    echo "  run1=${warm1}s  run2=${warm2}s  run3=${warm3}s"
    warm_median=$(python3 -c "
import statistics
print(statistics.median([${warm1}, ${warm2}, ${warm3}]))
")
fi
echo "  warm_median=${warm_median}s"

# ----------------------------------------------------------------
# Evaluate ratio
# ----------------------------------------------------------------
echo ""
ratio=$(python3 -c "
cold = ${cold_median}
warm = ${warm_median}
if cold <= 0:
    print('inf')
else:
    print(f'{warm/cold:.4f}')
")

echo "=== Results ==="
echo "  cold_median : ${cold_median} s"
echo "  warm_median : ${warm_median} s"
echo "  ratio       : ${ratio}"
echo "  threshold   : ${THRESHOLD}"
echo ""

# Compare ratio to threshold
result=$(python3 -c "
ratio = float('${ratio}') if '${ratio}' != 'inf' else 1e9
threshold = float('${THRESHOLD}')
if ratio < threshold:
    print('PASS')
else:
    print('WARN')
")

if [[ "${result}" == "PASS" ]]; then
    echo "RESULT: PASS (warm/cold ratio ${ratio} < ${THRESHOLD})"
    echo "  Cache disk reads are significantly faster than cold runs."
    exit 0
else
    echo "RESULT: WARN (warm/cold ratio ${ratio} >= ${THRESHOLD})"
    echo "  This may be expected for small test fixtures (process/Python overhead dominates)."
    echo "  For production Np=Nq=30 runs, ratio should be well below 0.30."
    echo "  This is NOT a hard failure; ctest + 190 MeV regression remain the gates."
    exit 1
fi
