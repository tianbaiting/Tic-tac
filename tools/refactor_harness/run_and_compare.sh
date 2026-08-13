#!/usr/bin/env bash
# Bitwise regression harness for the Tic-tac structural refactor.
#
# Runs two deterministic reduced-grid solves (2NF-only and complete-factorized-3NF)
# and compares every output artifact byte-for-byte against the recorded baseline.
#
# Usage:
#   tools/refactor_harness/run_and_compare.sh          # run + compare
#   tools/refactor_harness/run_and_compare.sh --record # run + (re)write baseline
#
# The physics contract (docs/three_nf_equation_contract.md) is NOT exercised here
# for correctness; this harness only guards BITWISE stability of deterministic
# solver artifacts across the structural refactor. Operator-algebra, W1 cache
# parity, dense-vs-Pade, etc. remain covered by the CTest/Python suites.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
EXE="$REPO/CPP/run"
BASELINE="$HERE/baseline_hashes.txt"
RECORD=0
[[ "${1:-}" == "--record" ]] && RECORD=1

run_case() {
    local input="$1"; local out="$2"; local cache="$3"
    rm -rf "$HERE/$out" "$HERE/$cache"
    mkdir -p "$HERE/$out"
    # Override the IO paths inline so each case lands in an isolated dir;
    # the solver accepts key=value overrides on the command line.
    ( cd "$REPO" && "$EXE" "$HERE/$input" \
        "output_folder=$HERE/$out" \
        "P123_folder=$HERE/$out" \
        "cache_root=$HERE/$cache" ) >/dev/null 2>&1
}

hash_dir() {
    local dir="$1"
    ( cd "$dir" && find . -type f ! -name '*.h5' -printf '%P\n' | sort \
        | xargs -r sha256sum )
}

run_case input_golden_2nf.txt  out_run_2nf  cache_run_2nf
run_case input_golden_3nf.txt  out_run_3nf  cache_run_3nf

CURRENT="$HERE/current_hashes.txt"
: > "$CURRENT"
{ echo "# 2NF-only reduced-grid solve"; hash_dir "$HERE/out_run_2nf"; } >> "$CURRENT"
{ echo "# complete-factorized-3NF reduced-grid solve"; hash_dir "$HERE/out_run_3nf"; } >> "$CURRENT"

if [[ $RECORD -eq 1 ]]; then
    cp "$CURRENT" "$BASELINE"
    echo "RECORDED baseline -> $BASELINE"
    exit 0
fi

if [[ ! -f "$BASELINE" ]]; then
    echo "ERROR: no baseline at $BASELINE (run with --record first)" >&2
    exit 2
fi

if diff -u "$BASELINE" "$CURRENT"; then
    echo "BITWISE REGRESSION: PASS (all deterministic artifacts identical)"
    exit 0
else
    echo "BITWISE REGRESSION: FAIL (deterministic artifacts changed)" >&2
    exit 1
fi
