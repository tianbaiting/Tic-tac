# Tic-tac & tictac-origin 内部技术长文 —— 设计规范

**日期**：2026-04-25
**作者**：Baiting Tian（与 Claude 协作）
**状态**：设计已经经过用户分节审阅，进入实现阶段

---

## 1. 目标

为实验室内部读者沉淀一份覆盖 `Tic-tac`（当前 fork：`tianbaiting/Tic-tac`）与 `tictac-origin`（原版：`seanbsm/Tic-tac`）两个仓库的长文技术参考，要求同时讲清：

1. 物理（量子三体散射 + Faddeev/AGS + Wave-Packet Continuum Discretization 方法）
2. 算法（连续谱离散化、置换算符 P123、resolvent、Padé 重求和）
3. 代码数值化（C++ / Fortran 90 / Fortran 77 / Python 在两个仓库中的真实落地与差异）
4. 上手实战（build → 190 MeV 验证 → 多能量观测量 → 3NF sweep）

读者前置假设：研究生量子力学 + 一点散射理论；能看懂 C++ / Fortran 源码。
读者目标姿态：读完后能独立修改 / 调试求解器、并知道每一段公式对应哪一行代码。

## 2. 关键设计抉择（已与用户确认）

| 选项 | 决定 |
|---|---|
| 与现有 `docs/faddeev/md_blocks/` 9 章中文长文的关系 | **完全重写**，不复用。新工程独立目录。 |
| 文档语言 | **纯中文** |
| 排版 | **LaTeX**，xelatex + ctex 编译，必须能产出 PDF |
| 结构轴 | **物理驱动 + 代码并排**：每章末尾"代码落地"小节并排展示 origin / current 代码差异 |
| 具体化等级 | **完全具体化**：真实代码片段 + 数字闭环 + 真实验证图；不写伪代码冒充真代码 |
| 仓库演化处理 | 每章末尾局部对照 + 全局两章总览（第 18 / 19 章） |
| 3NF 处理 | **专设第 V 部分**（第 15–17 章），覆盖物理推导、PW 投影代码、数值正确性踩坑 |
| 工作量分摊 | 大量使用 subagent，按依赖分波次执行 |

## 3. 全书目录架构

**书名（暂拟）**：《三体核散射的 Wave-Packet Faddeev/AGS 数值实现 —— Tic-tac & tictac-origin 内部技术参考》

**预估规模**：350–500 页（A4 + 中文正文 11pt + 代码 listing 9pt）

```
卷首
  扉页 / 版权 / 致谢 / 序：本文档定位与阅读路径
  符号与单位约定一览表
  读者前置准备

第 I 部分  物理与算法基础（Foundation）
  第 1 章  量子散射的算符结构与 Lippmann–Schwinger 方程
  第 2 章  动量表象中的两体散射数值化
  第 3 章  三体问题、Jacobi 坐标与部分波展开
  第 4 章  Faddeev 分解、AGS 方程与 U-算符

第 II 部分  连续谱离散化（Discretization）
  第 5 章  Wave-Packet Continuum Discretization (WPCD) 数学基础
  第 6 章  WP / SWP 基矢的代码构造（make_wp_states / make_swp_states）
  第 7 章  部分波对称态空间（make_pw_symm_states）
  第 8 章  两体势矩阵在 WP 基下的构造（make_potential_matrix）

第 III 部分  三体核心机器（Three-Body Machinery）
  第 9 章   置换算符 P123 的解析与代码实现（make_permutation_matrix）
  第 10 章  通道 Resolvent 的 R+Q 分解（make_resolvent）
  第 11 章  Faddeev 求解器：Neumann 级数 + Padé 重求和（solve_faddeev）

第 IV 部分  从 U-amplitude 到可观测量（Observables）
  第 12 章  U-算符到弹性散射振幅与散射截面
  第 13 章  极化观测量（iT11, T20, T21, T22, Ay）的部分波合成
  第 14 章  Miller Gate 1 相移提取（用户新工作）

第 V 部分  三核子力专章（Three-Nucleon Force, 3NF）
  第 15 章  3NF 物理：1π-exchange contact、2π-exchange、contact term
  第 16 章  3NF 部分波投影的代码实现（rank-0 + rank-2 分解、c_D / c_1 / c_3 / c_E）
  第 17 章  3NF 数值正确性：Fourier 归一化 1/(2π)³、c_E Navrátil/Witała 符号、check_3nf_normalization 工具

第 VI 部分  仓库演化与上手实战（Repo Evolution & Hands-on）
  第 18 章  tictac-origin 代码骨架解读（Sean 原版 baseline）
  第 19 章  Tic-tac 重构与新增功能（src/ 分层、CMake 路径、Python workflows）
  第 20 章  上手实战路径（build → 190 MeV 验证 → 多能量观测量 → 3NF sweep）

附录
  A. 物理常数与单位换算（与 include/constants.h 对齐）
  B. 角动量耦合系数与 Wigner 符号（CG / 6j / 9j）
  C. Gauss–Legendre 求积与 Chebyshev 非均匀网格
  D. GSL / LAPACK / HDF5 接口最简手册
  E. 完整输入参数字典（与 src/config/set_run_parameters.cpp 同步）
  F. 已知问题、踩坑记录、未来工作

参考文献（BibTeX）
  - Sean B.S. Miller 博士论文 + 4 篇 PRC/JPG
  - Glöckle / Witała / Epelbaum 等核心文献
  - 相关数学物理教科书
```

## 4. 每章统一模板（8 段式）

每一章按以下骨架写，保证读者建立稳定阅读节奏：

| 段 | 段名 | 功能 | 典型篇幅 |
|---|---|---|---|
| §x.1 | 物理动机 | 这一步为什么必须存在；上下游边界 | 1–2 页 |
| §x.2 | 数学推导 | 从已建立的物理出发，一步步推到可数值化的形式（含中间步骤、不跳跃） | 3–10 页 |
| §x.3 | 离散化方案 | 连续量 → 离散数组：网格、求积、截断；说清近似阶 | 2–5 页 |
| §x.4 | 算法骨架 | 算法级伪码（输入 / 输出 / 复杂度 / 内存）；用 `algorithm2e` | 1–2 页 |
| §x.5 | 代码落地 | 真实代码片段 + 行号 + 文件路径；origin / current 双栏并排（diff 章节）；逐行注释关键行 | 3–8 页 |
| §x.6 | 数字闭环 | 一个具体输入 → 给出关键中间数组的真实数值（表）+ 验证图（PDF） | 1–3 页 |
| §x.7 | 接口与依赖 | 这一模块对上游/下游模块的 I/O 契约；TikZ 调用图 | 1 页 |
| §x.8 | 常见陷阱 | 容易写错的地方，含真实 commit 教训 | 1–2 页 |

**特殊章节偏离模板**：
- 第 1–4 章（物理基础）：以 §x.1–§x.4 为主；§x.5 代码落地缩短
- 第 18, 19 章（仓库演化）：完全不走模板；专做仓库地图 + 模块映射表 + 重构 motivation
- 第 20 章（上手实战）：折叠为"实验场景 → 命令 → 期望输出 → 故障排查"四段

**写作硬性约束**：
1. 不写伪代码冒充真代码——所有 §x.5 的代码必须来自真实文件，行号准确（通过 `extract.py` 抽取）
2. 不写"读者自证"——所有数学推导给完整中间步骤
3. 不写"通过验证"——§x.6 必须有具体可复现命令 + 真实输出数值
4. 中文段落 + 英文专业术语首次出现时双语标注（如：散射算符 (T-operator)）
5. 公式在文中编号，跨章引用用 `\cref{eq:LS-form}` 而非"上式"

## 5. LaTeX 工程结构

位置：`Tic-tac/docs/treatise/`（独立 LaTeX 工程，不与现有 `docs/faddeev/` 冲突）。

```
Tic-tac/docs/treatise/
├── main.tex                      # 主文档：\input chapters & appendices
├── preamble.tex                  # 全部宏包加载与全局设置
├── physics_macros.tex            # 自定义符号 (\Tlab, \Np, \Nq, \Pop, ...)
├── code_macros.tex               # listings 风格定义 (cpp/fortran/python)
├── chapters/                     # 21 个 .tex（00 preface + 01–20）
├── appendices/                   # 6 个 .tex（A–F）
├── figures/
│   ├── tikz/                     # 调用图、流程图源
│   ├── plots/                    # matplotlib 输出 PDF
│   └── snapshots/                # 数字闭环表（CSV → .tex 自动转）
├── code_excerpts/
│   ├── extract.py                # 给定 (file, line_start, line_end) 抽取活代码
│   └── snippets/                 # 输出供 \lstinputlisting 使用
├── references.bib                # Sean 论文 + 经典文献
├── latexmkrc                     # xelatex pipeline
├── Makefile                      # make / make clean / make snippets / make figures
└── README.md                     # 构建说明
```

**宏包栈**：`ctex` + `xeCJK` + `listings` + `xcolor` + `algorithm2e` + `tikz` + `pgfplots` + `biblatex` + `biber` + `cleveref` + `physics` + `bm` + `mathtools` + `siunitx`

**复现性硬约束（写进 Makefile / CI）**：
1. 所有代码片段通过 `extract.py` 从真实仓库源文件抽出（带 git SHA pin），不允许手写复制
2. 所有 §x.6 数字闭环对应一个 `figures/snapshots/<chapter>_<topic>.csv` + 一个生成命令注释
3. `make build` 必须一次性产出 `treatise.pdf`，xelatex 退出码 0
4. 提供 `make watch`（latexmk -pvc）便于增量预览

## 6. 实施波次

22 章 + 6 附录 ≈ 28 个独立 .tex 文件。按依赖严格分波：

| 波 | 内容 | 章数 | 并行度 |
|---|---|---|---|
| 0 | 工程骨架（preamble、Makefile、extract.py、main.tex 占位、references.bib 种子） | — | 1（Claude 主线） |
| 1 | 第 I 部分：第 1–4 章 | 4 | 2 agents |
| 2 | 第 II 部分：第 5–8 章 | 4 | 2 agents |
| 3 | 第 III 部分：第 9–11 章 | 3 | 3 agents |
| 4 | 第 IV 部分：第 12–14 章 | 3 | 3 agents |
| 5 | 第 V 部分：第 15–17 章 | 3 | 3 agents |
| 6 | 第 VI 部分：第 18–20 章 | 3 | 3 agents |
| 7 | 附录 A–F | 6 | 3 agents |
| 8 | 整体编译、交叉引用修复、润色 | — | 1 agent |

每个 subagent 拿到：
- 完整的本设计文档（章节清单 + 8 段模板）
- preamble.tex / physics_macros.tex 当前版本（保证符号一致）
- 依赖章节的关键标签清单（"第 N 章已经写了 `\cref{eq:foo}`，你可以引用"）
- 自己负责章节对应的真实代码文件路径列表
- 写作硬约束清单（同 §4 末尾）
- 禁止改 preamble；用 `extract.py` 抽代码；§x.6 必须给可复现命令

## 7. 验收标准

1. `make build` 成功产出 `treatise.pdf`，xelatex 退出码 0、warning ≤ 50（主要来自 underfull hbox）
2. 每章具备 8 段（或明确说明的偏离模板）
3. 每章 §x.5 至少 1 段真实代码（origin 或 current），通过 `extract.py` 抽取
4. 每章 §x.6 至少 1 个可复现命令 + 1 个真实数值表/图
5. 全文交叉引用无 `??`
6. 参考文献至少包含 Sean 博士论文 + 4 篇论文 + Glöckle 等 5 部经典著作
7. 全书页数在 300–600 之间（具体值不做硬指标，但不应低于 250 页或超过 700 页）

## 8. 风险与缓解

| 风险 | 缓解 |
|---|---|
| Subagent 写出的章节符号不一致 | preamble.tex / physics_macros.tex 提前定义 + 给每个 agent 强制注入 |
| 真实代码漂移导致 listing 错位 | extract.py 用 git SHA pin + 行号校验 |
| 引用未来章节产生悬空 \cref | 依赖严格按波次执行；Wave 8 集中修 |
| 数字闭环命令需要可执行环境 | Wave 6/7 标注哪些命令在当前 CI/本地可重跑 |
| 篇幅失控 | 每章篇幅上限写进 prompt；Wave 8 整理时削减冗余 |

## 9. 后续

本规范完成后立即进入实施：

1. **Wave 0**（Claude 主线）：建立 LaTeX 工程骨架
2. **Wave 1**（并行 subagents）：第 1–4 章
3. 执行完毕后向用户汇报，等待触发 Wave 2

后续波次按相同模式推进。
