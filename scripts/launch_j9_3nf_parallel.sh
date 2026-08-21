#!/usr/bin/env bash
# Launch the 7 remaining j9+3NF channels IN PARALLEL (one process per channel).
# Each gets its own output_folder to avoid shared-file collisions.
set -u
cd /data/tian/workspace/dpol/Tic-tac
BASE=output/ay_diagnosis_current_binary/j9_3nf/par
INPUT=CPP/Input/input_j9_3nf.txt
# chn_3N (0-indexed) values still missing after the serial run did 0,1,2
for chn in 3 4 5 6 7 8 9; do
  outdir="$BASE/ch$chn/solver_out"
  mkdir -p "$outdir"
  logdir="$BASE/ch$chn"
  setsid bash -c "./CPP/run $INPUT parallel_run=true channel_idx=$chn output_folder=$outdir > $logdir/run.log 2>&1" < /dev/null > /dev/null 2>&1 &
done
echo "launched 7 parallel channel solves"
