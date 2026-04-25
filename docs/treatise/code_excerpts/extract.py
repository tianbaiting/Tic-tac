#!/usr/bin/env python3
"""
extract.py
==========

工具：从仓库源文件按行号区段抽取代码片段，写入 treatise/code_excerpts/snippets/
为 LaTeX listings 提供 \\lstinputlisting 输入。

设计目标：
    1. 任何 §x.5 代码落地小节必须通过本工具抽取，不允许手写复制
    2. 输出片段携带原始文件路径与行号注释（首行注释）
    3. 支持 origin / current 两个仓库

用法
----
单条抽取：
    ./extract.py --repo current --src src/core/state_space/make_permutation_matrix.cpp \\
                 --start 120 --end 180 --out p123_core.cpp

通过 manifest 批量抽取：
    ./extract.py --manifest manifest.yaml

manifest.yaml 形如：
    snippets:
      - out: ch09_p123_origin.cpp
        repo: origin
        src: CPP/make_permutation_matrix.cpp
        start: 120
        end: 180
      - out: ch09_p123_current.cpp
        repo: current
        src: src/core/state_space/make_permutation_matrix.cpp
        start: 120
        end: 180
"""

from __future__ import annotations
import argparse
import os
import sys
from pathlib import Path

# ----------------------------------------------------------------------
# 路径常量
# ----------------------------------------------------------------------
TREATISE_DIR = Path(__file__).resolve().parent.parent
DPOL_ROOT = TREATISE_DIR.parent.parent.parent          # /home/tian/workspace/dpol
REPO_ROOTS = {
    "current": DPOL_ROOT / "Tic-tac",
    "origin":  DPOL_ROOT / "tictac-origin" / "Tic-tac",
}
SNIPPETS_DIR = TREATISE_DIR / "code_excerpts" / "snippets"


def extract(repo: str, src: str, start: int, end: int, out: str) -> Path:
    """抽取代码片段，返回输出文件路径。"""
    if repo not in REPO_ROOTS:
        raise SystemExit(f"unknown repo '{repo}'; expected one of {list(REPO_ROOTS)}")

    src_path = REPO_ROOTS[repo] / src
    if not src_path.is_file():
        raise SystemExit(f"source not found: {src_path}")

    with src_path.open("r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    if start < 1 or end > len(lines) or start > end:
        raise SystemExit(
            f"invalid range [{start}..{end}] for {src_path} "
            f"(file has {len(lines)} lines)"
        )

    SNIPPETS_DIR.mkdir(parents=True, exist_ok=True)
    out_path = SNIPPETS_DIR / out

    ext = Path(src).suffix.lstrip(".")
    comment = "// " if ext in {"cpp", "h", "hpp", "cc"} else ("# " if ext in {"py", "f", "f90"} else "// ")
    if ext in {"f", "f90"}:
        comment = "! "

    header = (
        f"{comment}===============================================================\n"
        f"{comment}抽取自仓库 [{repo}]: {src}\n"
        f"{comment}行号区段：{start}..{end}\n"
        f"{comment}由 docs/treatise/code_excerpts/extract.py 自动生成；勿手动编辑\n"
        f"{comment}===============================================================\n"
    )

    body = "".join(lines[start - 1 : end])
    out_path.write_text(header + body, encoding="utf-8")
    return out_path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--repo", choices=list(REPO_ROOTS), help="repo: current or origin")
    ap.add_argument("--src", help="path inside the repo")
    ap.add_argument("--start", type=int, help="start line (1-indexed, inclusive)")
    ap.add_argument("--end", type=int, help="end line (inclusive)")
    ap.add_argument("--out", help="output filename inside snippets/")
    ap.add_argument("--manifest", help="YAML manifest for batch extraction")
    args = ap.parse_args()

    if args.manifest:
        try:
            import yaml  # type: ignore
        except ImportError:
            print("[extract] PyYAML not installed; install via 'pip install pyyaml'", file=sys.stderr)
            return 2
        with open(args.manifest, "r", encoding="utf-8") as f:
            manifest = yaml.safe_load(f) or {}
        items = manifest.get("snippets", [])
        for item in items:
            out = extract(item["repo"], item["src"], int(item["start"]), int(item["end"]), item["out"])
            print(f"[extract] {item['repo']}:{item['src']} [{item['start']}..{item['end']}] -> {out}")
        return 0

    needed = [args.repo, args.src, args.start, args.end, args.out]
    if not all(x is not None for x in needed):
        ap.error("either --manifest, or all of --repo --src --start --end --out are required")

    out = extract(args.repo, args.src, args.start, args.end, args.out)
    print(f"[extract] -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
