#!/usr/bin/env bash

# Purpose:
#   Convenience entrypoint for common example workflows.
#
# Data flow (high level):
#   command token
#     -> call one or more python scripts under examples/
#     -> produce outputs under output/*
#     -> optional plotting through micromamba env
#
# Calls:
#   - quick_Ay_test.py
#   - deuteron_proton_Ay.py
#   - compare_Ay_experiment.py
#   - run_dpol_p_observables.py
#   - plot_validation_curves.py / plot_dpol_p_observables.py
#
# Usage:
#   ./examples/run_examples.sh smoke
#   ./examples/run_examples.sh validate190
#   ./examples/run_examples.sh multi
#   ./examples/run_examples.sh clean

set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PYTHON_BIN="${PYTHON_BIN:-python3}"

print_info() {
    echo "[INFO] $1"
}

print_warn() {
    echo "[WARN] $1"
}

run_smoke() {
    # Fast pipeline health check: solver + comparison (and optional summary parsing).
    print_info "Run quick Tlab=190 MeV smoke validation"
    "$PYTHON_BIN" "$ROOT_DIR/examples/quick_Ay_test.py" --work-dir output/quick_Ay_test
}

run_validate_190() {
    # Single-energy Tlab=190 MeV run and validation artifact generation.
    print_info "Run Tlab=190 MeV solver and experiment comparison"
    "$PYTHON_BIN" "$ROOT_DIR/examples/deuteron_proton_Ay.py" \
        --work-dir output/deuteron_proton_Ay \
        --target-tlab-mev 190

    "$PYTHON_BIN" "$ROOT_DIR/examples/compare_Ay_experiment.py" \
        --work-dir output/deuteron_proton_Ay \
        --solver-out-dir output/deuteron_proton_Ay/solver_out \
        --target-tlab-mev 190

    if command -v micromamba >/dev/null 2>&1; then
        print_info "Plot with micromamba env anaroot-env"
        micromamba run -n anaroot-env python "$ROOT_DIR/examples/plot_validation_curves.py" \
            --work-dir output/deuteron_proton_Ay || print_warn "Plot command failed"
    else
        print_warn "micromamba not found, skip plotting"
    fi
}

run_multi_energy() {
    # Multi-energy production path for 70/135/190 MeV target Tlab values.
    print_info "Run multi-energy observables (target Tlab = 70/135/190 MeV)"
    "$PYTHON_BIN" "$ROOT_DIR/examples/run_dpol_p_observables.py" \
        --work-dir output/dpol_p_observables \
        --target-tlabs-mev 70,135,190

    if command -v micromamba >/dev/null 2>&1; then
        print_info "Plot multi-energy observables"
        micromamba run -n anaroot-env python "$ROOT_DIR/examples/plot_dpol_p_observables.py" \
            --work-dir output/dpol_p_observables || print_warn "Plot command failed"
    else
        print_warn "micromamba not found, skip plotting"
    fi
}

clean_outputs() {
    # Remove generated outputs only; source code remains untouched.
    print_info "Clean output folder"
    rm -rf "$ROOT_DIR/output"/*
}

show_help() {
    cat <<'EOF'
Tic-tac example runner

Usage:
  ./examples/run_examples.sh <command>

Commands:
  smoke        quick Tlab=190 MeV smoke validation
  validate190  run Tlab=190 MeV validation workflow
  multi        run 70/135/190 MeV target-Tlab observables workflow
  clean        remove files under output/
  help         show this help
EOF
}

case "${1:-help}" in
    smoke)
        run_smoke
        ;;
    validate190)
        run_validate_190
        ;;
    multi)
        run_multi_energy
        ;;
    clean)
        clean_outputs
        ;;
    help|--help|-h)
        show_help
        ;;
    *)
        print_warn "Unknown command: $1"
        show_help
        exit 1
        ;;
esac
