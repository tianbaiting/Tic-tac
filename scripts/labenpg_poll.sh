#!/bin/bash
# Poll labenpg 3NF run progress. Usage: bash labenpg_poll.sh
set -e
ssh labenpg-hk 'cd /data/tian/workspace/dpol/Tic-tac
echo "=== uptime / load ==="
uptime
echo
echo "=== per-chn progress ==="
DONE=0
TOTAL=0
for c in 4 5 6 7 8 9; do
  d=CPP/Output/labenpg_3NF_J3N9_par_chn$c
  TOTAL=$((TOTAL+1))
  u=$(ls $d/U_PW_elements_*.txt 2>/dev/null | wc -l)
  if [ "$u" -ge 1 ]; then DONE=$((DONE+1)); fi
  pade_steps=$(grep -c "Working on Pade approximant" $d/run.log 2>/dev/null || echo 0)
  kn_max=$(grep "A\*K^n for n=" $d/run.log 2>/dev/null | sed "s/.*n=//" | sed "s/[^0-9].*//" | sort -n | tail -1)
  [ -z "$kn_max" ] && kn_max=0
  echo "chn=$c: U_PW=$u  Padé_steps=$pade_steps  max_K_n=$kn_max  | last: $(tail -2 $d/run.log | tr "\r" "\n" | tail -1 | head -c 80)"
done
echo
echo "=== chn 0-3 (already on disk) ==="
ls CPP/Output/labenpg_3NF_J3N9/U_PW_elements_*.txt 2>/dev/null | wc -l
echo
echo "=== procs running ==="
N=$(ps -eo cmd | grep "CPP/run" | grep -v grep | wc -l)
echo "CPP/run procs: $N"
ps -eo pid,ni,etime,pcpu,pmem,rss,cmd --sort=-rss | grep "CPP/run" | grep -v grep | head -10
echo
echo "=== U_PW completion: $DONE/$TOTAL chn 4-9 done ==="
'
