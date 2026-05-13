#!/usr/bin/env python3
"""
Purpose:
  Drive the Tic-tac solver in three 3NF configurations at a single Tlab and
  collect U-matrix outputs side-by-side so the 3NF effect on observables can
  be visualised with plot_3nf_effect.py.

The three configurations:
  1) no_3nf     — three_nucleon_force=none (reference, must match 2NF baseline)
  2) zero_lec   — three_nucleon_force=chiral_N2LO with c_D=c_E=0 (sanity: W^(1)
                  from c_i still nonzero but small; reduces to pure 2NF only when
                  the 2NF model also has c1=c3=c4=0, e.g. LO_internal)
  3) witala     — three_nucleon_force=chiral_N2LO with Witala typical LECs
                  (c_D=-0.20, c_E=-0.205, Λ_3NF=500 MeV; PRC 77 (2008) 034004)

Shared P123 cache:
  All configs share the same state space, so P123 is computed once in the first
  run and reused via calculate_and_store_P123=false in subsequent runs.

Usage:
  python3 examples/run_3nf_sweep.py \
    --target-tlab-mev 64.5 \
    --potential-model N2LOopt \
    --work-dir output/3nf_sweep_64p5
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import List, Optional

from tlab_utils import format_tlab_label


@dataclass
class SweepConfig:
    name: str                   # directory label under work-dir
    three_nucleon_force: str    # "none" or "chiral_N2LO"
    c_D: float = 0.0
    c_E: float = 0.0
    Lambda_3NF: float = 500.0
    description: str = ""


DEFAULT_SWEEP: List[SweepConfig] = [
    SweepConfig(
        name="no_3nf",
        three_nucleon_force="none",
        description="baseline: 2NF only",
    ),
    SweepConfig(
        name="zero_lec",
        three_nucleon_force="chiral_N2LO",
        c_D=0.0,
        c_E=0.0,
        description="3NF path active with c_D=c_E=0; isolates c_i 2PE contribution",
    ),
    SweepConfig(
        name="witala",
        three_nucleon_force="chiral_N2LO",
        c_D=-0.20,
        c_E=-0.205,
        description="Witala et al. PRC 77 (2008) 034004 typical LECs",
    ),
]


@dataclass
class RunResult:
    name: str
    returncode: int
    solver_out_dir: str
    log_file: str
    config_file: str
    u_files: List[str] = field(default_factory=list)
    has_u_data: bool = False


def _rel_to_cpp(path: Path, root: Path) -> str:
    cpp_dir = root / "CPP"
    return os.path.relpath(path.resolve(), cpp_dir.resolve())


def _ensure_solver(root: Path, solver_path: Path) -> None:
    if solver_path.exists():
        return
    if solver_path != root / "CPP" / "run":
        raise FileNotFoundError(f"Solver binary not found: {solver_path}")
    result = subprocess.run(
        ["make", "-C", str(root / "CPP"), "-j"],
        cwd=root,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "Failed to build CPP/run\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )


def _write_config(
    *,
    root: Path,
    target_tlab_mev: float,
    two_j_3n_max: int,
    j_2n_max: int,
    np_wp: int,
    nq_wp: int,
    nphi: int,
    nx: int,
    np_per_wp: int,
    nq_per_wp: int,
    chebyshev_s: float,
    chebyshev_t: float,
    threads: int,
    potential_model: str,
    sweep: SweepConfig,
    shared_p123_dir: Path,
    solver_out_dir: Path,
    calculate_p123: bool,
    w1_scale: float,
) -> tuple[Path, Path]:
    cpp_input_dir = root / "CPP" / "Input"
    cpp_input_dir.mkdir(parents=True, exist_ok=True)

    tlab_label = format_tlab_label(target_tlab_mev)
    slug = f"sweep_{sweep.name}_tlab_{tlab_label}"

    energy_file = cpp_input_dir / f"tlab_{tlab_label}_3nf_sweep_autogen.txt"
    if not energy_file.exists():
        energy_file.write_text(f"{target_tlab_mev:.6f}\n", encoding="utf-8")

    config_file = cpp_input_dir / f"input_{slug}_autogen.txt"
    solver_out_dir.mkdir(parents=True, exist_ok=True)
    shared_p123_dir.mkdir(parents=True, exist_ok=True)

    lines = [
        f"# 3NF sweep config: {sweep.name}",
        f"# target Tlab {target_tlab_mev:.6f} MeV",
        f"two_J_3N_max={two_j_3n_max}",
        f"J_2N_max={j_2n_max}",
        f"Np_WP={np_wp}",
        f"Nq_WP={nq_wp}",
        f"Nphi={nphi}",
        f"Nx={nx}",
        f"Np_per_WP={np_per_wp}",
        f"Nq_per_WP={nq_per_wp}",
        f"chebyshev_s={chebyshev_s}",
        f"chebyshev_t={chebyshev_t}",
        "p_grid_type=chebyshev",
        "q_grid_type=chebyshev",
        f"P123_omp_num_threads={threads}",
        "parallel_run=false",
        "P123_recovery=false",
        "tensor_force=true",
        "isospin_breaking_1S0=true",
        "midpoint_approx=false",
        f"calculate_and_store_P123={'true' if calculate_p123 else 'false'}",
        "include_breakup_channels=false",
        "solve_faddeev=true",
        "solve_dense=false",
        "production_run=true",
        f"potential_model={potential_model}",
        "parameter_walk=false",
        "parameter_file=",
        "PSI_start=-1",
        "PSI_end=-1",
        f"three_nucleon_force={sweep.three_nucleon_force}",
        f"c_D={sweep.c_D}",
        f"c_E={sweep.c_E}",
        f"Lambda_3NF={sweep.Lambda_3NF}",
        f"w1_scale={w1_scale}",
        f"energy_input_file={_rel_to_cpp(energy_file, root)}",
        f"output_folder={_rel_to_cpp(solver_out_dir, root)}",
        f"P123_folder={_rel_to_cpp(shared_p123_dir, root)}",
        "",
    ]
    config_file.write_text("\n".join(lines), encoding="utf-8")
    return energy_file, config_file


def _scan_u_files(solver_out_dir: Path) -> tuple[List[Path], bool]:
    u_files = sorted(solver_out_dir.glob("U_PW_elements_*.txt"))
    has_data = False
    for u_file in u_files:
        for raw_line in u_file.read_text(encoding="utf-8").splitlines():
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 7:
                continue
            try:
                float(parts[0])
                float(parts[1])
                int(parts[2])
                complex(parts[3])
                has_data = True
                break
            except Exception:
                continue
        if has_data:
            break
    return u_files, has_data


def run_single_config(
    *,
    root: Path,
    solver_path: Path,
    work_dir: Path,
    target_tlab_mev: float,
    sweep: SweepConfig,
    shared_p123_dir: Path,
    calculate_p123: bool,
    timeout_s: int,
    w1_scale: float,
    **solver_kwargs,
) -> RunResult:
    solver_out_dir = work_dir / sweep.name / "solver_out"
    energy_file, config_file = _write_config(
        root=root,
        target_tlab_mev=target_tlab_mev,
        sweep=sweep,
        shared_p123_dir=shared_p123_dir,
        solver_out_dir=solver_out_dir,
        calculate_p123=calculate_p123,
        w1_scale=w1_scale,
        **solver_kwargs,
    )

    log_file = work_dir / sweep.name / "solver_run.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)

    cpp_cwd = root / "CPP"
    solver_rel = os.path.relpath(solver_path.resolve(), cpp_cwd.resolve())
    if not solver_rel.startswith("."):
        solver_rel = f"./{solver_rel}"
    config_rel = os.path.relpath(config_file.resolve(), cpp_cwd.resolve())

    env = os.environ.copy()
    env["HDF5_DISABLE_VERSION_CHECK"] = "2"

    with log_file.open("w", encoding="utf-8") as log_handle:
        result = subprocess.run(
            [solver_rel, config_rel],
            cwd=cpp_cwd,
            env=env,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=timeout_s,
        )

    u_files, has_u_data = _scan_u_files(solver_out_dir)
    return RunResult(
        name=sweep.name,
        returncode=result.returncode,
        solver_out_dir=str(solver_out_dir),
        log_file=str(log_file),
        config_file=str(config_file),
        u_files=[str(p) for p in u_files],
        has_u_data=has_u_data,
    )


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run Tic-tac across three 3NF configurations at one Tlab")
    p.add_argument("--work-dir", default="output/3nf_sweep")
    p.add_argument("--solver", default="CPP/run")
    p.add_argument("--target-tlab-mev", type=float, default=64.5)
    p.add_argument("--two-j-3n-max", type=int, default=1)
    p.add_argument("--j-2n-max", type=int, default=2)
    p.add_argument("--np", type=int, default=20)
    p.add_argument("--nq", type=int, default=20)
    p.add_argument("--nphi", type=int, default=24)
    p.add_argument("--nx", type=int, default=24)
    p.add_argument("--np-per-wp", type=int, default=6)
    p.add_argument("--nq-per-wp", type=int, default=6)
    p.add_argument("--chebyshev-s", type=float, default=100.0)
    p.add_argument("--chebyshev-t", type=float, default=1.0)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--potential-model", default="N2LOopt",
                   help="2NF model (N2LOopt / Idaho_N3LO / LO_internal); chiral_N2LO 3NF inherits c_i from here")
    p.add_argument("--timeout", type=int, default=7200)
    p.add_argument("--c-d", type=float, default=-0.20, help="Witala c_D override")
    p.add_argument("--c-e", type=float, default=-0.205, help="Witala c_E override")
    p.add_argument("--lambda-3nf-mev", type=float, default=500.0)
    p.add_argument("--w1-scale", type=float, default=1.0,
                   help="Diagnostic scaling on W^(1) (1.0 = physical; smaller values tame Neumann divergence)")
    p.add_argument(
        "--reuse-p123",
        action="store_true",
        help="Skip P123 calculation (requires existing cache in work_dir/p123)",
    )
    p.add_argument(
        "--only",
        default="",
        help="Comma-separated subset of run names to execute (default: all)",
    )
    return p


def main() -> int:
    args = build_parser().parse_args()
    root = Path(__file__).resolve().parents[1]
    work_dir = (root / args.work_dir).resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    solver = (root / args.solver).resolve() if not Path(args.solver).is_absolute() else Path(args.solver)
    _ensure_solver(root, solver)

    shared_p123 = work_dir / "p123"

    # Tune Witala config from CLI values
    sweeps: List[SweepConfig] = []
    for s in DEFAULT_SWEEP:
        if s.name == "witala":
            sweeps.append(SweepConfig(
                name=s.name,
                three_nucleon_force=s.three_nucleon_force,
                c_D=args.c_d,
                c_E=args.c_e,
                Lambda_3NF=args.lambda_3nf_mev,
                description=s.description,
            ))
        else:
            sweeps.append(s)

    if args.only.strip():
        keep = {x.strip() for x in args.only.split(",") if x.strip()}
        sweeps = [s for s in sweeps if s.name in keep]
        if not sweeps:
            raise SystemExit(f"No sweeps matched --only='{args.only}'")

    solver_kwargs = dict(
        two_j_3n_max=args.two_j_3n_max,
        j_2n_max=args.j_2n_max,
        np_wp=args.np,
        nq_wp=args.nq,
        nphi=args.nphi,
        nx=args.nx,
        np_per_wp=args.np_per_wp,
        nq_per_wp=args.nq_per_wp,
        chebyshev_s=args.chebyshev_s,
        chebyshev_t=args.chebyshev_t,
        threads=args.threads,
        potential_model=args.potential_model,
    )

    results: List[RunResult] = []
    failed: List[str] = []
    for i, sweep in enumerate(sweeps):
        # Compute P123 only on first config if user didn't pass --reuse-p123.
        calc_p123 = (i == 0) and (not args.reuse_p123)
        print(f"[3nf_sweep] ({i+1}/{len(sweeps)}) running {sweep.name}: "
              f"tnf={sweep.three_nucleon_force} c_D={sweep.c_D} c_E={sweep.c_E} "
              f"Λ={sweep.Lambda_3NF} (calc_P123={calc_p123})",
              flush=True)
        try:
            res = run_single_config(
                root=root,
                solver_path=solver,
                work_dir=work_dir,
                target_tlab_mev=args.target_tlab_mev,
                sweep=sweep,
                shared_p123_dir=shared_p123,
                calculate_p123=calc_p123,
                timeout_s=args.timeout,
                w1_scale=args.w1_scale,
                **solver_kwargs,
            )
        except subprocess.TimeoutExpired:
            print(f"[3nf_sweep]   !! {sweep.name} TIMED OUT after {args.timeout}s; skipping",
                  flush=True)
            failed.append(sweep.name)
            continue

        results.append(res)
        print(f"[3nf_sweep]   -> rc={res.returncode} u_files={len(res.u_files)} has_data={res.has_u_data}",
              flush=True)
        if res.returncode != 0:
            print(f"[3nf_sweep]   !! {sweep.name} exited with rc={res.returncode}; see {res.log_file}",
                  flush=True)
            failed.append(sweep.name)

    summary = {
        "target_tlab_mev": args.target_tlab_mev,
        "potential_model": args.potential_model,
        "two_j_3n_max": args.two_j_3n_max,
        "j_2n_max": args.j_2n_max,
        "np_wp": args.np,
        "nq_wp": args.nq,
        "runs": [asdict(r) for r in results],
        "sweeps": [asdict(s) for s in sweeps],
        "shared_p123_dir": str(shared_p123),
    }
    summary["failed_runs"] = failed
    summary["w1_scale"] = args.w1_scale
    summary_file = work_dir / f"sweep_summary_tlab_{format_tlab_label(args.target_tlab_mev)}.json"
    summary_file.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"[3nf_sweep] summary -> {summary_file}")
    if failed:
        print(f"[3nf_sweep] NOTE: {len(failed)} run(s) failed: {failed}")
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
