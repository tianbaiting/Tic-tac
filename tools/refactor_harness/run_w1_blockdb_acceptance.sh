#!/usr/bin/env bash
# Phase D-J acceptance: the W1 block database is exact, resumable, distributable,
# and consumed by the solver with zero W1 re-evaluation. Runs on the cheap
# 4-sector grid (two_J_3N_max=3, Np=4, Nq=3). ~3 minutes.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
SOLVER="$REPO/CPP/run"
WORKER="$REPO/build/bin/w1_worker"
IN="$HERE/input_j3nf_multiblock.txt"
fail=0; chk(){ if eval "$1"; then echo "PASS  $2"; else echo "FAIL  $2"; fail=1; fi; }

# Fresh workspace.
for d in cold_cache cold_out warm_cache warm_out shard_cache rest_cache; do rm -rf "$HERE/$d"; done
mkdir -p "$HERE/cold_out" "$HERE/warm_out"

echo "=== J.1 small-grid oracle: inline W1 vs worker-built W1 give identical U ==="
# (a) solver builds W1 inline into its own cache
( cd "$REPO" && "$SOLVER" "$IN" output_folder="$HERE/cold_out" P123_folder="$HERE/cold_out" cache_root="$HERE/cold_cache" ) >/dev/null 2>&1
# (b) worker pre-builds the SAME blocks, then solver runs warm
"$WORKER" build "$IN" cache_root="$HERE/warm_cache" >/dev/null 2>&1
( cd "$REPO" && "$SOLVER" "$IN" output_folder="$HERE/warm_out" P123_folder="$HERE/warm_out" cache_root="$HERE/warm_cache" ) >/dev/null 2>&1
for JP in 1_1 1_-1; do
  f="U_PW_elements_Np_4_Nq_3_JP_${JP}_Jmax_1_PSI_0.txt"
  chk "diff -q '$HERE/cold_out/$f' '$HERE/warm_out/$f' >/dev/null 2>&1" "oracle  $f : inline==worker-built"
done

echo "=== G  cache-hit == zero evaluation: warm W1 build time ~ 0 ==="
warm_build=$( ( cd "$REPO" && "$SOLVER" "$IN" output_folder="$HERE/warm_out" P123_folder="$HERE/warm_out" cache_root="$HERE/warm_cache" ) 2>&1 \
              | grep -oE 'build=[0-9.]+ s' | head -1 )
echo "   warm W1 build line: $warm_build"
chk "echo '$warm_build' | grep -q 'build=0.0 s'" "warm W1 build is 0.0 s (zero evaluation)"

echo "=== J.2 restart: a completed sector is rebuilt with zero evaluation ==="
"$WORKER" build "$IN" cache_root="$HERE/rest_cache" --sector 1 1 >/dev/null 2>&1
# Worker signal: it prints "[3NF W1] ... evaluating N of M missing ..." ONLY when
# it evaluates blocks. An all-cache-hit rebuild prints none.
restart_eval=$( "$WORKER" build "$IN" cache_root="$HERE/rest_cache" --sector 1 1 2>&1 \
               | grep -cE "evaluating [0-9]+ of" || true )
fresh_eval=$( "$WORKER" build "$IN" cache_root="$HERE/rest_cache" --sector 1 -1 2>&1 \
             | grep -cE "evaluating [0-9]+ of" || true )
echo "   re-built completed sector 'evaluating' lines: $restart_eval (want 0)"
echo "   fresh sector 'evaluating' lines:              $fresh_eval (want >0)"
chk "[[ '$restart_eval' -eq 0 && '$fresh_eval' -gt 0 ]]" "restart: completed sector reused (0 eval), fresh sector still evaluates"

echo "=== F  distributed: two disjoint shards build both active sectors ==="
"$WORKER" build "$IN" cache_root="$HERE/shard_cache" --shard 0/2 >/dev/null 2>&1
"$WORKER" build "$IN" cache_root="$HERE/shard_cache" --shard 1/2 >/dev/null 2>&1
stored=$( "$WORKER" plan "$IN" cache_root="$HERE/shard_cache" 2>/dev/null | grep -oE 'stored: [0-9]+' | grep -oE '[0-9]+' )
chk "[[ '$stored' -gt 0 ]]" "distributed shards populated the cache (stored=$stored blocks)"

echo "=== D  verify the assembled database ==="
vres=$( "$WORKER" verify "$IN" cache_root="$HERE/shard_cache" 2>&1 | grep -oE 'verify: [0-9]+ blocks checked, [0-9]+ missing' || true )
echo "   ${vres:-verify: (no summary line)}"

for d in cold_cache cold_out warm_cache warm_out shard_cache rest_cache; do rm -rf "$HERE/$d"; done
if [[ $fail -eq 0 ]]; then echo "PHASE D-J ACCEPTANCE: PASS"; exit 0
else echo "PHASE D-J ACCEPTANCE: FAIL" >&2; exit 1; fi
