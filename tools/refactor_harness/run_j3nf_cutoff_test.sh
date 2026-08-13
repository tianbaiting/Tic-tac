#!/usr/bin/env bash
# Phase B acceptance: prove the J_3NF cutoff selects the correct kernel per block.
#
# Runs a 4-sector grid (J=1/2,3/2 x parity) three ways on IDENTICAL numerical
# settings, then verifies byte-for-byte:
#   (a) full-3NF        (two_J_3NF_force_max=-1): every sector gets 3NF
#   (b) cutoff=1        (two_J_3NF_force_max=1):  J=1/2 3NF, J=3/2 pure-2NF
#   (c) pure-2NF        (three_nucleon_force=none): every sector 2NF
#
# Required equalities (the physics contract of the cutoff):
#   J=3/2 U in (b) == J=3/2 U in (c)   [high-J is exactly the 2NF kernel]
#   J=1/2 U in (b) == J=1/2 U in (a)   [low-J keeps the full 3NF]
# And: (b)'s W1 cache contains ONLY J=1/2 sectors (no high-J W1 was built).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
EXE="$REPO/CPP/run"
IN="$HERE/input_j3nf_multiblock.txt"

run() { # tag  override...
  local tag="$1"; shift
  rm -rf "$HERE/out_$tag" "$HERE/cache_$tag"
  mkdir -p "$HERE/out_$tag"
  ( cd "$REPO" && "$EXE" "$IN" \
      "output_folder=$HERE/out_$tag" \
      "P123_folder=$HERE/out_$tag" \
      "cache_root=$HERE/cache_$tag" "$@" ) >/dev/null 2>&1
}

echo "[Phase B] running full-3NF (a) ...";     run a  two_J_3NF_force_max=-1
echo "[Phase B] running cutoff=1 (b) ...";     run b  two_J_3NF_force_max=1
echo "[Phase B] running pure-2NF (c) ...";     run c  three_nucleon_force=none

# Equality helpers.
same() { diff -q "$HERE/out_$1/$3" "$HERE/out_$2/$4" >/dev/null 2>&1; }

U() { echo "U_PW_elements_Np_4_Nq_3_JP_${1}_${2}_Jmax_1_PSI_0.txt"; }

fail=0
# (1) high-J (J=3/2, both parities): cutoff == pure-2NF
for P in 1 -1; do
  if same c b "$(U 3 $P)" "$(U 3 $P)"; then
    echo "PASS  J=3/2,P=$P: cutoff(b) == pure-2NF(c)"
  else echo "FAIL  J=3/2,P=$P: cutoff(b) != pure-2NF(c)"; fail=1; fi
done
# (2) low-J (J=1/2, both parities): cutoff == full-3NF
for P in 1 -1; do
  if same a b "$(U 1 $P)" "$(U 1 $P)"; then
    echo "PASS  J=1/2,P=$P: cutoff(b) == full-3NF(a)"
  else echo "FAIL  J=1/2,P=$P: cutoff(b) != full-3NF(a)"; fail=1; fi
done
# (3) cutoff run must NOT equal full-3NF on J=3/2 (proves the cutoff actually
#     changed the high-J kernel, i.e. 3NF was contributing there in (a))
for P in 1 -1; do
  if ! same a b "$(U 3 $P)" "$(U 3 $P)"; then
    echo "PASS  J=3/2,P=$P: cutoff(b) != full-3NF(a)  [cutoff changed the high-J result]"
  else echo "FAIL  J=3/2,P=$P: cutoff(b) == full-3NF(a) unexpectedly"; fail=1; fi
done
# (4) W1 cache in (b) has ONLY J=1/2 sectors (no high-J W1 built)
n_j32_w1=$(ls "$HERE/cache_b/w1/" 2>/dev/null | grep -cE "JP3[-+][0-9]+" || true)
n_j12_w1=$(ls "$HERE/cache_b/w1/" 2>/dev/null | grep -cE "JP1[-+][0-9]+" || true)
if [[ "$n_j32_w1" -eq 0 && "$n_j12_w1" -gt 0 ]]; then
  echo "PASS  W1 cache: J=1/2 sectors present ($n_j12_w1 files), J=3/2 sectors absent (0 files)"
else echo "FAIL  W1 cache: J=3/2 files=$n_j32_w1, J=1/2 files=$n_j12_w1"; fail=1; fi

if [[ $fail -eq 0 ]]; then
  echo "PHASE B CUTOFF TEST: PASS"; exit 0
else
  echo "PHASE B CUTOFF TEST: FAIL" >&2; exit 1
fi
