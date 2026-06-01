#!/bin/bash
# Re-launch 6 parallel 3NF runs with P123 caches reused (calculate_and_store_P123=false).
# Stagger launch by 60s to spread memory pressure during initial loading.
set -e
cd /data/tian/workspace/dpol/Tic-tac

for chn in 4 5 6 7 8 9; do
    OUTDIR=CPP/Output/labenpg_3NF_J3N9_par_chn${chn}
    test -f $OUTDIR/input.txt || { echo "missing input $OUTDIR/input.txt"; exit 1; }
    test -f $OUTDIR/P123_sparse_JP_*_*_Np_20_Nq_20_J2max_3.h5 || { echo "missing P123 cache $OUTDIR"; exit 1; }

    # Update input to read existing P123 (no re-compute, no over-write)
    NEW=$OUTDIR/input_resume.txt
    sed -e 's|^calculate_and_store_P123=.*|calculate_and_store_P123=false|' \
        -e 's|^P123_recovery=.*|P123_recovery=false|' \
        $OUTDIR/input.txt > $NEW

    # Truncate previous log so we can tell new from old
    > $OUTDIR/run.log

    nohup env OMP_NUM_THREADS=14 ./CPP/run $NEW \
        > $OUTDIR/run.log 2>&1 &
    PID=$!
    renice +10 -p $PID >/dev/null 2>&1
    echo "Re-launched chn=${chn} pid=$PID (renice +10)"
    sleep 60   # stagger launches
done

sleep 5
echo "---all CPP/run processes---"
ps -eo pid,ni,etime,pcpu,pmem,cmd | grep "CPP/run" | grep -v grep
echo "---mem snapshot---"
free -h | head -3
