#!/bin/bash
# =============================================================================
# Sean-Miller-compatible pure-2NF convergence campaign
#
# Runs three independent convergence ladders with three_nucleon_force=none:
#   Phase 3: J_2N_max = 1, 2, 3, 4  (at Np=82, Nq=12, two_J_3N_max=1)
#   Phase 4: two_J_3N_max = 1, 3, 5, 7, 9, 11, 13, 15, 17
#            (at converged J_2N_max, Np=82, Nq=12)
#   Phase 5: Nq = 8, 12, 16, 24, 32  (at converged J_2N_max, two_J_3N_max, Np=82)
#
# RESTARTABLE: each rung checks for existing U files before running.
# All work is pure 2NF (three_nucleon_force=none).
#
# Usage:  bash tools/convergence/run_sean_2nf_campaign.sh
# =============================================================================
set -euo pipefail

REPO="/data/tian/workspace/dpol/Tic-tac"
SOLVER="$REPO/CPP/run"
INPUT="$REPO/CPP/Input/input_realistic_82_12.txt"
ENERGIES="$REPO/CPP/Input/lab_energies.txt"
CAMPAIGN="$REPO/output/sean_2nf_convergence"
STATUS_FILE="$CAMPAIGN/status/status.json"
LOG_FILE="$CAMPAIGN/logs/campaign.log"
GIT_SHA=$(cd "$REPO" && git rev-parse --short HEAD)
HOSTNAME=$(hostname)
OMP_THREADS=${OMP_NUM_THREADS:-48}

export OMP_NUM_THREADS=$OMP_THREADS

cd "$REPO"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

# ---- status tracking ---------------------------------------------------------
init_status() {
    if [ ! -f "$STATUS_FILE" ]; then
        cat > "$STATUS_FILE" << 'EOF'
{
  "tmux_session": "tictac-2nf-conv",
  "hostname": "",
  "git_sha": "",
  "omp_threads": 0,
  "j2n_ladder": {},
  "j3n_ladder": {},
  "nq_ladder": {},
  "last_updated": ""
}
EOF
    fi
}

update_status() {
    # $1 = ladder name, $2 = rung, $3 = status, $4 = extra_json
    python3 -c "
import json, sys, datetime
from pathlib import Path
sf = Path('$STATUS_FILE')
d = json.loads(sf.read_text())
d['hostname'] = '$HOSTNAME'
d['git_sha'] = '$GIT_SHA'
d['omp_threads'] = $OMP_THREADS
ladder = sys.argv[1]
rung = sys.argv[2]
status = sys.argv[3]
extra = json.loads(sys.argv[4]) if len(sys.argv) > 4 and sys.argv[4] else {}
if ladder not in d:
    d[ladder] = {}
d[ladder][rung] = {'status': status, 'updated': datetime.datetime.now().isoformat(), **extra}
d['last_updated'] = datetime.datetime.now().isoformat()
sf.write_text(json.dumps(d, indent=2))
" "$1" "$2" "$3" "${4:-}"
}

# ---- run one solver rung -----------------------------------------------------
# Args: $1=ladder $2=rung_label $3=work_dir $4=p123_folder $5..=overrides
# Checks for existing U files; skips if both parities present.
run_rung() {
    local ladder="$1"
    local rung="$2"
    local work_dir="$3"
    local p123_folder="$4"
    shift 4
    local overrides=("$@")

    local out_dir="$work_dir/out"
    local log="$work_dir/run.log"
    mkdir -p "$out_dir"

    # Check if both parities' U files exist
    local u_pos=$(ls "$out_dir"/U_PW_elements_*JP_1_1*_PSI_0.txt 2>/dev/null | head -1)
    local u_neg=$(ls "$out_dir"/U_PW_elements_*JP_1_-1*_PSI_0.txt 2>/dev/null | head -1)
    if [ -n "$u_pos" ] && [ -n "$u_neg" ]; then
        log "$rung: SKIP (both parities already have U files)"
        update_status "$ladder" "$rung" "DONE" '{"skipped": true}'
        return 0
    fi

    # Check if P123 already exists (for P123_recovery)
    local p123_exists=""
    if ls "$p123_folder"/P123_sparse_*J2max*.h5 >/dev/null 2>&1; then
        p123_exists="true"
    fi

    log "$rung: START (P123 exists: ${p123_exists:-no})"

    # Build command
    local cmd_args=(
        "$SOLVER" "$INPUT"
        "output_folder=$out_dir"
        "P123_folder=$p123_folder"
        "energy_input_file=$ENERGIES"
        "three_nucleon_force=none"
        "P123_omp_num_threads=$OMP_THREADS"
        "pade_max_order=24"
    )
    # Add P123 reuse if P123 files exist (read HDF5 directly)
    if [ -n "$p123_exists" ]; then
        cmd_args+=("calculate_and_store_P123=false")
    fi
    cmd_args+=("${overrides[@]}")

    local t0=$(date +%s)
    /usr/bin/time -v "${cmd_args[@]}" > "$log" 2>&1 || {
        local rc=$?
        local t1=$(date +%s)
        log "$rung: FAILED (rc=$rc, $((t1-t0))s)"
        update_status "$ladder" "$rung" "FAILED" "{\"wall_seconds\": $((t1-t0)), \"rc\": $rc}"
        return 1
    }
    local t1=$(date +%s)
    local wall=$((t1-t0))

    # Extract peak RSS from /usr/bin/time output
    local peak_rss=$(grep "Maximum resident" "$log" | grep -o '[0-9]*' || echo "0")
    local u_count=$(ls "$out_dir"/U_PW_elements_*_PSI_0.txt 2>/dev/null | wc -l)

    log "$rung: DONE (${wall}s, ${u_count} U files, RSS ${peak_rss} kB)"
    update_status "$ladder" "$rung" "DONE" \
        "{\"wall_seconds\": $wall, \"u_files\": $u_count, \"peak_rss_kb\": $peak_rss, \"p123_reused\": ${p123_exists:-false}}"
    return 0
}

# =============================================================================
# Phase 3: J_2N_max ladder (Np=82, Nq=12, two_J_3N_max=1, pure 2NF)
# =============================================================================
run_j2n_ladder() {
    log "========== Phase 3: J_2N_max ladder =========="
    local base="$CAMPAIGN/j2n_ladder"

    # Reuse existing J2N=1 and J2N=2 outputs from the previous campaign
    local prev_base="$REPO/output/realistic_2nf_convergence/j2n_ladder"

    for j2n in 1 2 3 4; do
        local work="$base/J2N_${j2n}"
        mkdir -p "$work/out"

        # Copy existing results from previous campaign if available
        local prev_rung="$prev_base/rung_$((j2n-1))_J2N=$j2n/out"
        if [ -d "$prev_rung" ]; then
            cp -n "$prev_rung/"* "$work/out/" 2>/dev/null || true
        fi

        # Run with P123 in the output dir (so P123_recovery can reuse)
        run_rung "j2n_ladder" "J2N_${j2n}" "$work" "$work/out" \
            "J_2N_max=$j2n" "two_J_3N_max=1"
    done
}

# =============================================================================
# Phase 4: two_J_3N_max ladder (at converged J2N, Np=82, Nq=12, pure 2NF)
# =============================================================================
run_j3n_ladder() {
    local j2n="${1:-3}"  # default to J2N=3, but can be overridden
    log "========== Phase 4: two_J_3N_max ladder (J2N=$j2n) =========="
    local base="$CAMPAIGN/j3n_ladder"

    for tj3 in 1 3 5 7 9 11 13 15 17; do
        local work="$base/J3N_${tj3}"
        mkdir -p "$work/out"

        # Copy P123 files from previous rung (P123 is per-sector, reusable).
        # IMPORTANT: do NOT copy U files — they must be regenerated for each
        # two_J_3N_max because the solver produces U for ALL sectors up to two_J_3N_max.
        local prev_tj3=$((tj3 - 2))
        if [ "$prev_tj3" -ge 1 ] && [ -d "$base/J3N_${prev_tj3}/out" ]; then
            cp -n "$base/J3N_${prev_tj3}/out/"P123_*.h5 "$work/out/" 2>/dev/null || true
        elif [ "$tj3" -eq 1 ] && [ -d "$CAMPAIGN/j2n_ladder/J2N_${j2n}/out" ]; then
            cp -n "$CAMPAIGN/j2n_ladder/J2N_${j2n}/out/"P123_*.h5 "$work/out/" 2>/dev/null || true
        fi

        # For the skip check: count expected sectors vs existing U files.
        # At two_J_3N_max=tj3, the number of J^pi sectors is roughly tj3+1
        # (both parities for each J=1/2, 3/2, ..., tj3/2).
        # We check for U files with the correct Jmax in the filename.
        local expected_jp_count=$(( (tj3 + 1) ))  # approximate
        local actual_u_count=$(ls "$work/out"/U_PW_elements_*_Jmax_${j2n}_PSI_0.txt 2>/dev/null | wc -l)

        if [ "$actual_u_count" -ge "$expected_jp_count" ]; then
            log "J3N_${tj3}: SKIP ($actual_u_count U files already present)"
            update_status "j3n_ladder" "J3N_${tj3}" "DONE" "{\"skipped\": true, \"u_files\": $actual_u_count}"
            continue
        fi

        run_rung "j3n_ladder" "J3N_${tj3}" "$work" "$work/out" \
            "J_2N_max=$j2n" "two_J_3N_max=$tj3"
    done
}

# =============================================================================
# Phase 5: Nq ladder (at converged J2N, two_J_3N_max, Np=82, pure 2NF)
# WARNING: changing Nq moves the on-shell midpoint — compare observables, not U
# =============================================================================
run_nq_ladder() {
    local j2n="${1:-3}"
    local tj3="${2:-1}"
    log "========== Phase 5: Nq ladder (J2N=$j2n, 2J3N=$tj3) =========="
    local base="$CAMPAIGN/nq_ladder"

    for nq in 8 12 16 24 32; do
        local work="$base/Nq_${nq}"
        mkdir -p "$work"
        run_rung "nq_ladder" "Nq_${nq}" "$work" "$work/out" \
            "J_2N_max=$j2n" "two_J_3N_max=$tj3" "Nq_WP=$nq"
    done
}

# =============================================================================
# Main
# =============================================================================
main() {
    init_status
    log "Sean-compatible 2NF convergence campaign starting"
    log "  git_sha=$GIT_SHA  hostname=$HOSTNAME  OMP=$OMP_THREADS"
    log "  input=$INPUT  energies=$ENERGIES"

    # Phase 3: J_2N_max ladder (HIGHEST PRIORITY)
    run_j2n_ladder

    # Phase 4: two_J_3N_max ladder (at J2N=3 or 4 depending on convergence)
    # Use J2N=3 as the starting point (may upgrade to 4 if J2N=3->4 shows convergence)
    run_j3n_ladder 3

    # Phase 5: Nq ladder (deferred until angular axes are settled)
    # run_nq_ladder 3 1

    log "Campaign complete."
    update_status "" "campaign" "COMPLETE"
}

main "$@"
