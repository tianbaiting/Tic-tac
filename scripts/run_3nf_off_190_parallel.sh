#!/bin/bash
# Launch all 10 channels (JP=1±1, 3±1, 5±1, 7±1, 9±1) in parallel for 3NF-off at 190 MeV.
# Reuses the existing P123 cache (nijmegen 2NF, Np=20, Nq=20, J_2N_max=3).
set -u
cd /data/tian/workspace/dpol/Tic-tac

PIDS=()
for chn in 0 1 2 3 4 5 6 7 8 9; do
    OUTDIR=CPP/Output/dpol190_np20_3nf_off_chn${chn}
    mkdir -p $OUTDIR
    # Per-channel input file: same as the parent, but with channel_idx=$chn and output_folder=$OUTDIR
    NEW=$OUTDIR/input.txt
    sed -e "s|^channel_idx=.*|channel_idx=${chn}|" \
        -e "s|^output_folder=.*|output_folder=${OUTDIR}|" \
        CPP/Input/input_dpol190_np20_3nf_off.txt > $NEW
    echo "Launching chn=${chn} -> $OUTDIR"
    env OMP_NUM_THREADS=8 ./CPP/run $NEW > $OUTDIR/run.log 2>&1 &
    PIDS+=($!)
done

echo "---waiting for ${#PIDS[@]} channels---"
FAIL=0
for i in "${!PIDS[@]}"; do
    if ! wait "${PIDS[$i]}"; then
        echo "  chn=${i} (pid ${PIDS[$i]}) FAILED"
        FAIL=$((FAIL + 1))
    fi
done
echo "---done. failures=${FAIL}---"
ls CPP/Output/dpol190_np20_3nf_off_chn*/U_PW_elements_*.txt 2>&1 | wc -l
