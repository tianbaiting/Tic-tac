# ===============================================================
# 抽取自仓库 [current]: examples/deuteron_proton_Ay.py
# 行号区段：148..203
# 由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑
# ===============================================================
def run_cpp_solver(cfg: SolverRunConfig) -> dict:
    ensure_cpp_solver_binary(cfg.root, cfg.solver)
    energy_file, config_file, solver_out_dir, log_file = write_solver_inputs(cfg)

    cpp_cwd = cfg.root / "CPP"
    solver_rel = os.path.relpath(cfg.solver.resolve(), cpp_cwd.resolve())
    if not solver_rel.startswith("."):
        solver_rel = f"./{solver_rel}"
    config_rel = os.path.relpath(config_file.resolve(), cpp_cwd.resolve())

    cmd = [solver_rel, config_rel]
    env = os.environ.copy()
    env["HDF5_DISABLE_VERSION_CHECK"] = "2"

    with log_file.open("w", encoding="utf-8") as log_handle:
        result = subprocess.run(
            cmd,
            cwd=cpp_cwd,
            env=env,
            stdout=log_handle,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=cfg.timeout_s,
        )

    u_files = sorted(solver_out_dir.glob("U_PW_elements_*.txt"))
    has_u_data = False
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
                has_u_data = True
                break
            except Exception:
                continue
        if has_u_data:
            break

    return {
        "returncode": result.returncode,
        "log_file": log_file,
        "energy_file": energy_file,
        "config_file": config_file,
        "solver_out_dir": solver_out_dir,
        "u_files": u_files,
        "has_u_data": has_u_data,
    }
