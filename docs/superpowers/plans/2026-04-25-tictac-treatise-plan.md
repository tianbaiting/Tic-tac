# Tic-tac & tictac-origin 长文 —— 实施计划

**关联设计文档**：`docs/superpowers/specs/2026-04-25-tictac-treatise-design.md`
**日期**：2026-04-25
**总目标**：在 `Tic-tac/docs/treatise/` 产出一份可编译的中文 LaTeX 长文，覆盖 Tic-tac 与 tictac-origin 两仓库的物理、算法、数值化与上手实战。

---

## 总览（依赖图）

```
Wave 0 (scaffold) ──> Wave 1 (Ch 1-4 物理基础)
                       │
                       └─> Wave 2 (Ch 5-8 离散化)
                            │
                            └─> Wave 3 (Ch 9-11 三体核心)
                                 │
                                 ├─> Wave 4 (Ch 12-14 可观测量)
                                 │
                                 └─> Wave 5 (Ch 15-17 3NF)
                                      │
                                      └─> Wave 6 (Ch 18-20 仓库演化)
                                           │
                                           ├─> Wave 7 (附录 A-F)
                                           │
                                           └─> Wave 8 (终极编译/润色)
```

每一波内部并行；波与波之间串行。

---

## Wave 0 — 工程骨架（Claude 主线）

**输出**：可编译的 LaTeX 工程骨架（main.tex 编译通过，但内容为占位符）

**任务**：
1. `mkdir -p Tic-tac/docs/treatise/{chapters,appendices,figures/{tikz,plots,snapshots},code_excerpts/snippets}`
2. 写 `preamble.tex`：加载所有宏包，配置 ctex 中文字体（Noto/SourceHan fallback）
3. 写 `physics_macros.tex`：常用物理符号定义（详见下文 "符号约定"）
4. 写 `code_macros.tex`：listings 风格定义（cpp / fortran / python，行号、彩色）
5. 写 `main.tex`：标题页 + 目录 + `\input{}` 占位每一章 + 附录 + 参考文献
6. 写 21 个章节占位 .tex（每个一行 `\chapter{标题}` + `% TODO Wave N`）
7. 写 6 个附录占位 .tex
8. 写 `code_excerpts/extract.py`：CLI 工具，输入 `(file, start_line, end_line)`，输出片段到 `snippets/`
9. 写 `references.bib`：Sean 博士论文 + 4 篇论文（占位条目，DOI 待填）+ 经典文献（Glöckle, Witała, Epelbaum, Navrátil）
10. 写 `latexmkrc`（xelatex 流水线）
11. 写 `Makefile`：`build / clean / snippets / figures / watch`
12. 写 `README.md`：构建说明
13. 验证：`make build` 产出 `treatise.pdf`（占位内容），xelatex 退出码 0

**符号约定（physics_macros.tex）**：

| 宏 | 渲染 | 含义 |
|---|---|---|
| `\Tlab` | $T_{\text{lab}}$ | 实验室能量 |
| `\Np`, `\Nq` | $N_p$, $N_q$ | WP 网格点数 |
| `\Pop` | $\hat{\mathcal{P}}$ | 置换算符 P123 |
| `\Top` | $\hat{T}$ | T-算符 |
| `\Uop` | $\hat{U}$ | AGS U-算符 |
| `\Vop` | $\hat{V}$ | 相互作用 |
| `\Gzero` | $G_0$ | 自由 Green 函数 |
| `\hbar c` | $\hbar c$ | 197.327 MeV·fm |
| `\bra{·}`, `\ket{·}` | $\langle·|$, $|·\rangle$ | 来自 physics 包 |
| `\jacobi{p}{q}` | $|p,q\rangle$ | Jacobi 动量本征态 |
| `\pwstate{α}` | $|\alpha\rangle$ | 部分波通道索引 |
| `\WPbin{i}` | $|x_i\rangle$ | 第 i 个 WP bin |

**完成判据**：`Tic-tac/docs/treatise/treatise.pdf` 存在且可打开；占位章节标题正确。

---

## Wave 1 — 第 1–4 章（物理与算法基础）

**前置**：Wave 0 完成
**并行度**：2 个 subagents

### Wave 1 / Agent A — 第 1–2 章

**负责**：
- 第 1 章「量子散射的算符结构与 Lippmann–Schwinger 方程」
- 第 2 章「动量表象中的两体散射数值化」

**对应代码**（§x.5 必须真实抽取）：
- 第 2 章：`src/core/potential/make_potential_matrix.cpp` 的两体势矩阵入口
- 第 2 章：`tictac-origin/Tic-tac/CPP/make_potential_matrix.cpp` 同段

**关键内容指引**：
- 第 1 章侧重物理基础：Møller 算符、T-算符、LS 方程的存在唯一性、入态/出态、boost 不变性；§1.5 简短引用代码（potential_model.h 的接口签名）
- 第 2 章重点：动量空间求积、Cutkosky 主值积分、相移提取、Numerov 在动量空间的等价物；§2.5 真实代码片段对照

**篇幅**：第 1 章 25–35 页，第 2 章 30–40 页。

### Wave 1 / Agent B — 第 3–4 章

**负责**：
- 第 3 章「三体问题、Jacobi 坐标与部分波展开」
- 第 4 章「Faddeev 分解、AGS 方程与 U-算符」

**对应代码**：
- 第 3 章：`src/core/state_space/make_pw_symm_states.cpp` 的通道索引构造
- 第 4 章：`src/core/faddeev_solver/solve_faddeev.h` 的接口签名

**关键内容指引**：
- 第 3 章：三种 Jacobi 树（12-3, 23-1, 31-2）、置换关系、$|αpq\rangle$ 态、自旋-同位旋耦合、反对称化
- 第 4 章：Faddeev 三分量、AGS 方程的 U_αβ 形式、与 LS 的对比、连接矩阵 1+P；§4.5 引代码接口签名

**篇幅**：第 3 章 30–40 页，第 4 章 25–35 页。

### Wave 1 共同约束

- 必须 `\input{preamble}` 不允许改它
- 必须用 `\cref` 系列；定义新公式标签时用 `\label{eq:lsequ}` 命名空间约定
- 引代码用 `\lstinputlisting{code_excerpts/snippets/...}`，不允许直接 `\begin{lstlisting}` 内联粘贴
- 物理符号必须用 `physics_macros.tex` 定义的宏；新宏先加进 `physics_macros.tex` 并在 commit message 里说明
- 提交时一并 commit `code_excerpts/snippets/` 内的抽取产物

---

## Wave 2 — 第 5–8 章（连续谱离散化）

**前置**：Wave 1 完成（章节交叉引用第 1–4 章公式）
**并行度**：2 个 subagents

### Wave 2 / Agent A — 第 5–6 章

- 第 5 章「WPCD 数学基础」：wave-packet 的定义、完备性、WP/SWP 的差别、Sturmian-Wave-Packet 推导
- 第 6 章「WP / SWP 基矢的代码构造」：`make_wp_states.cpp`、`make_swp_states.cpp` 双仓库对照

### Wave 2 / Agent B — 第 7–8 章

- 第 7 章「部分波对称态空间」：`make_pw_symm_states.cpp` 完整解读（J^π 块对角化、反对称化判据）
- 第 8 章「两体势矩阵在 WP 基下的构造」：`make_potential_matrix.cpp` + 5 种势模型的 §x.5 对照（chiral_LO_internal / N2LOopt / Idaho_N3LO / nijmegen / malfliet_tjon）

**篇幅**：每章 25–35 页。

---

## Wave 3 — 第 9–11 章（三体核心机器）

**前置**：Wave 2 完成
**并行度**：3 个 subagents（每人 1 章）

### Agent A — 第 9 章「置换算符 P123」
- `make_permutation_matrix.cpp` 解读
- 几何函数 G 的推导
- 稀疏存储与 HDF5 cache 机制
- §9.6 数字闭环：给一个 J^π=1/2^+ 的 P123 范数与稀疏度

### Agent B — 第 10 章「通道 Resolvent」
- `make_resolvent.cpp` 解读
- R + Q 分解、迭代 Möller 公式
- 复数 ε 处方与主值
- §10.6 数字闭环

### Agent C — 第 11 章「Faddeev 求解器」
- `solve_faddeev.cpp` 解读（Neumann 级数 + Padé 加速）
- 与 dense LAPACK 求解的对比模式
- 收敛诊断输出
- §11.6 数字闭环：190 MeV 在固定 J^π 通道上的 Padé 收敛序列

**篇幅**：每章 25–40 页。

---

## Wave 4 — 第 12–14 章（可观测量）

**前置**：Wave 3 完成
**并行度**：3 个 subagents

### Agent A — 第 12 章「U-算符到弹性振幅与截面」
### Agent B — 第 13 章「极化观测量」（iT11, T20, T21, T22, Ay）
### Agent C — 第 14 章「Miller Gate 1 相移提取」（用户最新 WIP 工作）

每章 §x.6 均要嵌入 PDF 矢量图（从 `examples/plot_*.py` 重新生成）+ CSV 数据表。

**篇幅**：每章 25–35 页。

---

## Wave 5 — 第 15–17 章（3NF 专章）

**前置**：Wave 4 完成（部分需要 Ch 11 的 solver 上下文）
**并行度**：3 个 subagents

### Agent A — 第 15 章「3NF 物理」
- 1π-exchange contact、2π-exchange (Fujita-Miyazawa)、contact term
- 低能常数 $c_D$, $c_1$, $c_3$, $c_E$ 的物理意义与 Navrátil/Witała 约定
- 与 Epelbaum N2LO 的对照

### Agent B — 第 16 章「3NF 部分波投影代码实现」
- rank-0 + rank-2 投影的完整推导
- `chiral_3nf_pw_kernels.h` / `three_nucleon_force_model.cpp` 解读
- 用户最近 5 个 commits 逐一映射到代码段：
  - bdb239d: 1/(2π)³ Fourier 归一化
  - 8db06bf: W1_1pe_contact rank-0 + rank-2 重写 (c_D)
  - 1f05031: W1_2pe rank-0 + rank-2 重写 (c_1, c_3)
  - f4df358: c_E 符号修正（Navrátil/Witała +½ 约定）

### Agent C — 第 17 章「3NF 数值正确性」
- `tools/check_3nf_normalization/` 工具完整解读
- 三体束缚态 ψ_t 的归一化与期望值检查
- 数值 sanity test：W1 期望值与 Witała benchmarks 对照
- §17.6 数字闭环：给出某 J=1/2 通道 W1 矩阵元的数值

**篇幅**：第 15 章 30–40 页，第 16 章 35–50 页（最长），第 17 章 25–35 页。

---

## Wave 6 — 第 18–20 章（仓库演化与上手实战）

**前置**：Wave 5 完成
**并行度**：3 个 subagents

### Agent A — 第 18 章「tictac-origin 代码骨架」
- 单层 CPP/ 树的模块映射表
- main.cpp 入口的执行流（步骤 1–10）
- Sean 论文与代码段的一一对应表

### Agent B — 第 19 章「Tic-tac 重构与新增功能」
- 从 origin 到 current 的逐文件 diff 总结
- src/ 三层架构（core / config / io / interactions / utils）的 motivation
- CMake 构建路径 vs Makefile 路径的 trade-off
- Python workflows 的设计（examples/ + tests/）
- 多能量观测量 + 3NF + Miller Gate 三大新增模块的高层概览

### Agent C — 第 20 章「上手实战路径」
- 4 个完整实验场景：
  1. 190 MeV 验证（`python examples/deuteron_proton_Ay.py`）
  2. 多能量观测量（`run_dpol_p_observables.py`）
  3. 3NF sweep（`examples/run_3nf_sweep.py`）
  4. Miller Gate 1 相移提取
- 每个场景：命令 → 期望输出 → 故障排查 → 性能预期

**篇幅**：第 18 章 30–40 页，第 19 章 30–45 页，第 20 章 25–35 页。

---

## Wave 7 — 附录 A–F

**前置**：Wave 6 完成
**并行度**：3 个 subagents

- Agent A：附录 A（常数）+ 附录 B（Wigner）
- Agent B：附录 C（求积）+ 附录 D（GSL/LAPACK/HDF5）
- Agent C：附录 E（参数字典，从 `set_run_parameters.cpp` 反向提取）+ 附录 F（已知问题）

**篇幅**：每个附录 10–25 页。

---

## Wave 8 — 终极编译、交叉引用修复、润色

**前置**：Wave 7 完成
**并行度**：1 个 subagent（也可 Claude 主线）

**任务**：
1. `make build` 全量编译，修所有 `??` 与未定义引用
2. 修所有 underfull/overfull hbox（最严重的 20 个）
3. 统一术语：抽取一个术语索引，扫所有章节，纠正不一致
4. 写"各章导读"放在每部分首页（Part I/II/III/.../VI）
5. 整理参考文献：补全 DOI、ISBN、URL
6. 最终页数核验、目录页码、索引页生成
7. 生成最终 `treatise.pdf` 的 git tag `treatise-v1.0`

---

## 共享约束（适用于所有 Wave 1–7 subagents）

每个 subagent prompt 必须包含：

1. **本设计文档完整路径**（让 agent 自己 Read）：
   - `docs/superpowers/specs/2026-04-25-tictac-treatise-design.md`
2. **本计划文档完整路径**（让 agent 看到自己的位置）：
   - `docs/superpowers/plans/2026-04-25-tictac-treatise-plan.md`
3. **依赖章节列表**：哪些 .tex 已经写完，可以 `\cref` 它的标签
4. **代码文件清单**：自己负责章节对应的真实代码路径
5. **不允许修改的文件清单**：preamble.tex / main.tex / Makefile / latexmkrc
6. **允许修改但需要 commit 时说明的文件**：physics_macros.tex（新加宏时）
7. **写作硬性约束**（与设计 §4 末尾一致）
8. **完成判据**：
   - 自己负责的 .tex 文件全部写完
   - `code_excerpts/snippets/` 内对应抽取产物存在
   - `cd Tic-tac/docs/treatise && xelatex -interaction=nonstopmode main.tex` 不报致命错误
   - commit 自己写的所有内容

## 时间预估（参考）

| Wave | 预估实际时间 |
|---|---|
| 0 | 30–60 分钟（Claude 主线） |
| 1 | 1–2 小时（2 个 agents 并行） |
| 2 | 1–2 小时 |
| 3 | 1.5–2.5 小时 |
| 4 | 1.5–2 小时 |
| 5 | 2–3 小时（最重） |
| 6 | 1.5–2.5 小时 |
| 7 | 1–1.5 小时 |
| 8 | 1–2 小时 |
| **合计** | **10–18 小时跨多次会话** |

不强制单次会话完成。每波结束 Claude 向用户汇报进度，用户决定是否触发下一波。

## 第一步行动

立即开始：
1. Wave 0 由 Claude 主线执行（不拆 agent）
2. Wave 0 完成后 commit
3. 触发 Wave 1（2 个并行 subagent）
4. 完成后向用户汇报，等待用户许可触发 Wave 2
