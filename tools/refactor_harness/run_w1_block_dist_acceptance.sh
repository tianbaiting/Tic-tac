#!/usr/bin/env bash
# Block-level distributed exact W^(1) construction acceptance.
#
# Demonstrates the production deliverable end-to-end on the cheap 4-sector grid:
#   plan -> 16 block-level workers in parallel -> status -> assemble (fingerprint)
#   interrupted run -> resume (zero recomputation of completed blocks)
#   exactness: block-level fingerprint == monolithic sector-level fingerprint
#
# The fingerprint is a SHA-256 over every evaluate-block payload in
# deterministic order, so identical fingerprints <=> bitwise-identical W^(1).
# ~3 minutes.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
WORKER="$REPO/build/bin/w1_worker"
SOLVER="$REPO/CPP/run"
IN="$HERE/input_j3nf_multiblock.txt"
fail=0; chk(){ if eval "$1"; then echo "PASS  $2"; else echo "FAIL  $2"; fail=1; fi; }

for d in bdist mono resume; do rm -rf "$HERE/$d"; done
mkdir -p "$HERE/bdist/out" "$HERE/mono/out" "$HERE/resume/out"

echo "=== 1. plan (work decomposition + signature) ==="
$WORKER plan "$IN" cache_root="$HERE/bdist/cache" --manifest "$HERE/bdist/manifest.json" 2>/dev/null | tail -4

echo "=== 2. build with 16 block-level workers in parallel ==="
for i in $(seq 0 15); do
  OMP_NUM_THREADS=6 $WORKER build "$IN" cache_root="$HERE/bdist/cache" \
    --worker-index $i --worker-count 16 \
    output_folder="$HERE/bdist/out" P123_folder="$HERE/bdist/out" >"$HERE/bdist/w$i.log" 2>&1 &
done
wait
echo "   built=$(grep -hoE 'built=[0-9]+' "$HERE/bdist"/w*.log | grep -oE '[0-9]+' | awk '{s+=$1} END{print s}') blocks across 16 workers"

echo "=== 3. status + assemble (fingerprint) ==="
$WORKER status "$IN" cache_root="$HERE/bdist/cache" 2>/dev/null | grep evaluate:
BDIST_FP=$($WORKER assemble "$IN" cache_root="$HERE/bdist/cache" 2>/dev/null | grep -oE 'W1_fingerprint=[0-9a-f]+')
echo "   block-level 16-worker fingerprint: $BDIST_FP"

echo "=== 4. exactness: block-level vs monolithic sector-level fingerprint ==="
$WORKER build "$IN" cache_root="$HERE/mono/cache" --shard 0/1 \
  output_folder="$HERE/mono/out" P123_folder="$HERE/mono/out" >/dev/null 2>&1
MONO_FP=$($WORKER assemble "$IN" cache_root="$HERE/mono/cache" 2>/dev/null | grep -oE 'W1_fingerprint=[0-9a-f]+')
echo "   monolithic sector-level fingerprint: $MONO_FP"
chk "[[ '$BDIST_FP' == '$MONO_FP' ]]" "block-level == monolithic W1 (bitwise fingerprint)"

echo "=== 5. downstream U equality: block-level cache vs monolithic cache -> solver ==="
( cd "$REPO" && "$SOLVER" "$IN" output_folder="$HERE/bdist/out" P123_folder="$HERE/bdist/out" cache_root="$HERE/bdist/cache" ) >/dev/null 2>&1
( cd "$REPO" && "$SOLVER" "$IN" output_folder="$HERE/mono/out" P123_folder="$HERE/mono/out" cache_root="$HERE/mono/cache" ) >/dev/null 2>&1
for JP in 1_1 1_-1; do
  f="U_PW_elements_Np_4_Nq_3_JP_${JP}_Jmax_1_PSI_0.txt"
  chk "diff -q '$HERE/bdist/out/$f' '$HERE/mono/out/$f' >/dev/null 2>&1" "U bitwise-identical  $f"
done

echo "=== 6. resume: Run A (worker 0 only) -> STOP -> Run B (all 16) -> no recomputation ==="
OMP_NUM_THREADS=6 $WORKER build "$IN" cache_root="$HERE/resume/cache" \
  --worker-index 0 --worker-count 16 \
  output_folder="$HERE/resume/out" P123_folder="$HERE/resume/out" >/dev/null 2>&1
A=$($WORKER status "$IN" cache_root="$HERE/resume/cache" 2>/dev/null | grep -oE 'evaluate: [0-9]+/[0-9]+' || true)
echo "   Run A (worker 0 only): $A"
for i in $(seq 0 15); do
  OMP_NUM_THREADS=6 $WORKER build "$IN" cache_root="$HERE/resume/cache" \
    --worker-index $i --worker-count 16 \
    output_folder="$HERE/resume/out" P123_folder="$HERE/resume/out" >"$HERE/resume/w$i.log" 2>&1 &
done
wait
B=$($WORKER status "$IN" cache_root="$HERE/resume/cache" 2>/dev/null | grep -oE 'evaluate: [0-9]+/[0-9]+' || true)
echo "   Run B (all 16):         $B"
W0=$(grep 'build:' "$HERE/resume/w0.log")
echo "   worker-0 resume line: $W0"
RESUME_FP=$($WORKER assemble "$IN" cache_root="$HERE/resume/cache" 2>/dev/null | grep -oE 'W1_fingerprint=[0-9a-f]+')
chk "echo '$W0' | grep -q 'built=0'" "resume: worker-0 built 0 (completed blocks skipped)"
chk "[[ '$B' == *'112/112'* ]]" "resume: campaign complete after restart"
chk "[[ '$RESUME_FP' == '$MONO_FP' ]]" "resumed W1 == monolithic W1 (bitwise)"

echo "=== 7. Hermitian contract (assemble reverse-block check) ==="
chk "$WORKER assemble "$IN" cache_root="$HERE/mono/cache" 2>/dev/null | grep -q 'Hermitian transpose pairs checked=90 bad=0'" \
     "monolithic: 90 reverse blocks == exact transpose of conjugate (bad=0)"

for d in bdist mono resume; do rm -rf "$HERE/$d"; done
if [[ $fail -eq 0 ]]; then echo "BLOCK-LEVEL DISTRIBUTED ACCEPTANCE: PASS"; exit 0
else echo "BLOCK-LEVEL DISTRIBUTED ACCEPTANCE: FAIL" >&2; exit 1; fi
