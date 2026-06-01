#!/bin/bash
# =====================================================================
# labenpg_run.sh — 复现 之前在 labenpg 上跑的 Tic-tac 3NF

#
# 用法：
#   bash labenpg_admin_setup.sh
#
# 跑起来后会在后台启 6 个 ./CPP/run 进程，nohup + renice +10 + 60s
# stagger。脚本本身退出，进程继续跑。后续监视方式跟 tian 之前一样：
#
#   ps -eo pid,ni,etime,pcpu,pmem,cmd | grep CPP/run | grep -v grep
#   for c in 4 5 6 7 8 9; do
#       d=/data/tian/workspace/dpol/Tic-tac/CPP/Output/labenpg_3NF_J3N9_par_chn$c
#       echo "chn=$c: $(ls $d/U_PW_elements_*.txt 2>/dev/null | wc -l) blocks"
#   done
# =====================================================================
set -e

CONDA_PREFIX_PATH=/data/tian/conda/envs/anaroot-env
TICTAC_DIR=/data/tian/workspace/dpol/Tic-tac

# ---------------------------------------------------------------------
# 1. conda env
# ---------------------------------------------------------------------
export CONDA_PREFIX="$CONDA_PREFIX_PATH"
export PATH="$CONDA_PREFIX/bin:$PATH"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"

# ---------------------------------------------------------------------
# 2. Build CPP/run（Makefile 的 LDLIBS := 硬赋值，必须命令行整体覆盖

# ---------------------------------------------------------------------
cd "$TICTAC_DIR/CPP"
make -j32 LDLIBS="-Wl,--no-as-needed -lgomp -lgsl -lpthread -lm -ldl \
    -lgfortran -lhdf5_hl_cpp -lhdf5_cpp -lhdf5_hl -lhdf5 -lstdc++fs \
    -llapacke -llapack -lcblas -lblas -lcurl"

# ---------------------------------------------------------------------
# 3. Launch 6 个并行 chn 4..9（J=5/2±, 7/2±, 9/2±）
#    P123 cache 已在盘上，用 calculate_and_store_P123=false 直接读
# ---------------------------------------------------------------------
cd "$TICTAC_DIR"

for chn in 4 5 6 7 8 9; do
    OUTDIR=CPP/Output/labenpg_3NF_J3N9_par_chn${chn}

    sed -e 's|^calculate_and_store_P123=.*|calculate_and_store_P123=false|' \
        -e 's|^P123_recovery=.*|P123_recovery=false|' \
        "$OUTDIR/input.txt" > "$OUTDIR/input_resume.txt"

    > "$OUTDIR/run.log"
    nohup env OMP_NUM_THREADS=14 ./CPP/run "$OUTDIR/input_resume.txt" \
        > "$OUTDIR/run.log" 2>&1 &
    PID=$!
    renice +10 -p $PID >/dev/null 2>&1 || true
    echo "Launched chn=${chn} pid=${PID} (renice +10)"
    sleep 60
done

# ---------------------------------------------------------------------
# 4. 启动确认
# ---------------------------------------------------------------------
sleep 5
echo
echo "== currently running =="
ps -eo pid,ni,etime,pcpu,pmem,cmd | grep "CPP/run" | grep -v grep
