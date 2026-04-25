# Tic-tac 长文工程

中文 LaTeX 长文，覆盖 `Tic-tac` 与 `tictac-origin` 两仓库的物理、算法、数值化与上手实战。

## 编译

```bash
cd Tic-tac/docs/treatise
make build       # 全量编译（xelatex + biber + xelatex x2）
xdg-open _build/main.pdf
```

依赖：

- `xelatex`（TeX Live 2023+）
- `biber`
- `latexmk`
- 中文字体：`Noto Serif CJK SC`、`Noto Sans CJK SC`、`Noto Sans Mono CJK SC`
- 西文字体：`TeX Gyre Termes`、`TeX Gyre Heros`、`DejaVu Sans Mono`

## 工程结构

| 路径 | 内容 |
|---|---|
| `main.tex` | 主文档；`\input` 各章 |
| `preamble.tex` | 宏包加载、字体、定理环境（**禁止 subagent 修改**） |
| `physics_macros.tex` | 物理符号宏（**只能添加新宏**） |
| `code_macros.tex` | listings 风格定义 |
| `chapters/00..20_*.tex` | 21 个章节 |
| `appendices/A..F_*.tex` | 6 个附录 |
| `figures/tikz/` | TikZ 调用图、流程图源 |
| `figures/plots/` | matplotlib 输出 PDF |
| `figures/snapshots/` | 数字闭环 CSV → 自动生成的 .tex |
| `code_excerpts/extract.py` | 代码片段抽取工具 |
| `code_excerpts/snippets/` | 抽取产物（供 `\lstinputlisting`） |
| `references.bib` | BibTeX 参考文献 |
| `latexmkrc` | latexmk 配置 |
| `Makefile` | 构建入口 |

## 写作约定

- 章节按 8 段式模板：物理动机 / 数学推导 / 离散化 / 算法骨架 / 代码落地 / 数字闭环 / 接口与依赖 / 常见陷阱
- 详见 `Tic-tac/docs/superpowers/specs/2026-04-25-tictac-treatise-design.md`
- 实施波次：`Tic-tac/docs/superpowers/plans/2026-04-25-tictac-treatise-plan.md`

## 抽取代码

```bash
./code_excerpts/extract.py \
    --repo current \
    --src src/core/state_space/make_permutation_matrix.cpp \
    --start 120 --end 180 \
    --out p123_core.cpp
```

或通过 manifest 批量抽取（YAML 格式见 `extract.py` 头注释）。

## 状态

- Wave 0 已完成：工程骨架可编译
- Wave 1+ 由 subagents 按计划填充章节
