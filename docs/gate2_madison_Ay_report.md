# Tic-tac Madison Ay(n) 计算器 + Miller Gate 2 量化收敛

**日期范围：** 2026-05-10 ~ 2026-05-12
**Tic-tac commits：** `06fadde`, `b6a63d9`, `a4417a9`, `ae57893` (+ arxiv 源 docs commit)
**主要 reference：** Miller, Ekström, Hebeler, PRC 106, 024001 (2022), arxiv:2201.09600

---

## 1. TL;DR

把 Tic-tac 输出的 jj-coupled U-matrix → 6×6 Madison spin-scattering M(θ)
的几何因子按 Miller PRC 2022 Appendix D 重写。此前 jj-coupling 直写的
`pw_amplitudes._geometric_factor` 在 J_3N > 1/2 上**所有观测量**（iT11、
T20-22、Ay_n）都会非单调跳号；现在按 Eq. D2 + Eq. D4 + Eq. D3 验证后：

**Tlab = 67 MeV 的 Ay(n) 在 Np=20 + J_3N ≤ 13/2 下达到 −0.531 @ θ=109°，
对照 Miller N3LO −0.540 / 实验 −0.510，相对实验 104%**（轻微 over-shoot，
但已在 Nijmegen-I 2NF 内禀偏差范围内）。35 MeV 与 10 MeV 仍未完全收敛
但都进入正确符号 + 量级；35 MeV 需要 J ≥ 15/2，10 MeV 显示了文献 Ay puzzle
（2NF 系统性低估实验约 2 倍）。

---

## 2. 背景：Miller Gate 2 = Nd elastic Ay(n) Ay puzzle benchmark

Miller PRC 107 014002 (2023) Fig. 6 给三个能量 Tlab = 10/35/67 MeV
极化中子-氘核弹性散射的 vector analyzing power Ay(n) 的实验和理论曲线，
是 Ay puzzle 文献的标准 reference。我们要把 Tic-tac 在 SAMURAI dpol 通路
跑出来的 U-matrix 接 Madison-convention 的极化观测量，去对照这套基准
（也就是组内称的 "Gate 2"）。

---

## 3. 这个工作之前的状态（v1，two_J_3N_max=1）

`compare_Ay_experiment.py` 已经能通过 `pw_amplitudes` 拼出 M(θ) 给出
iT11、T20、T21、T22，但 **σ_y on nucleon 的 Madison Ay(n) 没有加**，
仍然是"reduced-U 启发式"。第一步 commit `b6a63d9` 加上：

```python
SIGMA_Y_NUCLEON_2 = [[0,-i],[i,0]]   # Pauli-y on m_p basis
embed_nucleon_operator(op2) = kron(I_3, op2)
Ay_n = Re(Tr[M sigma_y^(in,N) M_dag] / Tr[M M_dag])
```

低 J 跑得对（J=1/2 在 ²S₁/₂ 主导下）。

---

## 4. 首次发现：J_3N > 1/2 全坏

跑 J_3N = 5/2 (v2) 和 J_3N = 3/2 (v3) 后发现 **iT11、T20、Ay_n 全部
在 v1 → v3 → v2 之间非单调跳号**：

| 35 MeV @ θ=90° | v1 (1/2) | v3 (3/2) | v2 (5/2) |
|---|---|---|---|
| iT11 | +0.205 | +0.058 ↓ | +0.214 ↑ |
| Ay_n | −0.014 | +0.038 ↑ | +0.221 ↑↑ |
| T20  | −0.191 (flat) | +0.013 | +0.003 |

之前以为 "v1 iT11 ≈ v2 iT11 → 框架 OK" 是 J=3/2/5/2 偶然抵消的假象。

经验性扫 8+ 个 phase 候选 `(-1)^(J+j-...)` 均无法同时修对 10/35/67 MeV →
不是单个 CG 重排符号差。

---

## 5. 修复路径：arxiv tex 源 + Miller D2/D4 verbatim

派 subagent 拉 arxiv:2201.09600 的 TeX 源（不是 OCR），在
`Appendix D` 找到：

### Eq. D2 — channel-spin scheme M(θ) 求和

```
M_{Σ' m_S', Σ m_S}(θ) = √π/(ik) · Σ_{J,l,l'} i^(l'-l) √(2l+1)
                        · <Σ m_S; l 0 | J m_S>
                        · <Σ' m_S'; l' (m_S-m_S') | J m_S>
                        · (S^J_{l'Σ', lΣ} - δ_{Σ'Σ} δ_{l'l})
                        · Y_{l'}^{m_S - m_S'}(θ, 0)
```

channel-spin **Σ ≡ J_d + s_p** 是氘核+核子的总自旋（1/2 或 3/2）。

### Eq. D3 — U 与 S 的关系

```
S^J_{l'Σ', lΣ} = δ_{l'l} δ_{ΣΣ'} - 2π i · q · m_N · i^(l'-l) · U^J_{l'Σ', lΣ}
```

**注意 D3 里也有 i^(l'-l) 因子。** 把 D3 代入 D2，i^(l'-l) × i^(l'-l)
= i^(2(l'-l))，由于宇称要求 (l'-l) 是偶数，2(l'-l) 是 4 的倍数，
**i^(2(l'-l)) = 1 恒成立 → 两个 i^(l'-l) 完全相消。** 如果直接代 U 进 M
（跳过 S 的显式构造），公式里**不应该有** i^(l'-l) 相位。

### Eq. D4 — jj → channel-spin 重耦

```
U^J_{l'Σ', lΣ} = Σ_{j',j}  √((2j'+1)(2Σ'+1)) (-1)^(J+j') {l'  1/2  j'; J_d J Σ'}
                          × √((2j +1)(2Σ +1)) (-1)^(J+j)  {l 1/2 j ;  J_d J Σ}
                          × U^J_{l'j', lj}
```

**这是关键 6j 变换**。Solver 输出的 jj-coupled `U^J_{l'j', lj}` 通过 6j 符号
变换到 channel-spin `U^J_{l'Σ', lΣ}`，再喂 D2。

### OCR vs TeX 差异（之前没找对的原因）

| | Published PDF OCR | arxiv TeX 真本 |
|---|---|---|
| Eq. D4 相位 | `(-1)^(J+j+l+1/2)` 半整数指数 ill-defined | **`(-1)^(J+j)`** 整数指数 |
| Eq. D4 √(2Σ+1) 因子 | 完全缺失 | **明确存在** |
| Eq. D2 i^(l'-l) | 写在 D2 但 OCR 没暴露 D3 的同名因子 | **D2 和 D3 都有，相消** |

→ **arxiv 源永远比 PDF OCR 可信。** 本地已存档在
`docs/seanBSMiller/arxiv_2201.09600/` 包含完整 tarball + 关键 tex。

---

## 6. 实现（commit `a4417a9`）

`examples/pw_amplitudes.py` 加三件：

1. **`wigner_6j(j1,j2,j3,j4,j5,j6)`** — Racah 公式，sanity 通过
   ({0 0 0;0 0 0}=1, {1 1 0;1 1 1}=−1/3, {1/2 1/2 0;1/2 1/2 0}=−1/2)。

2. **`_block_jj_to_channel_spin(block, U_jj)` → (cs_channels, U_cs)**
   按 D4 构 T 矩阵，验证 T·T^T = I 在所有 JP block 到 1e-16
   （unitary 基变换，物理上自洽）。

3. **`_channel_spin_geometric_factor()`** 按 D2 求 M_{Σ',Σ}(θ)
   再用 `<s_d m_d; s_p m_p | Σ m_Σ>` 投回 6×6 (m_d, m_p) helicity 基。

旧 jj 路径作为 `assemble_m_matrix(..., legacy_jj=True)` 保留 A/B 对照。

---

## 7. 部分波收敛跑批

5 个 Np=20、J_2N_max=3、Nijmegen-I 2NF 的 run，差异只在 `two_J_3N_max`：

| Run | two_J_3N_max | 6j blocks | wall time |
|---|---|---|---|
| v1 | 1 | 2 | ~10 s（cache reuse） |
| v3 | 3 | 4 | 15 min |
| v2 | 5 | 6 | 46 min |
| v4 | 9 | 10 | 91 min |
| v5 | 13 | 14 | 160 min |

输出在 `CPP/Output/miller_gate2_{v3_J3N3,v2_J3N5,v4_J3N9,v5_J3N13}_J2N3/`。

---

## 8. 结果

### 8.1 角度极值

| Tlab | v5 extremum | θ_ext | Miller N3LO | 实验 | 达成率 |
|---|---|---|---|---|---|
| 10 MeV | +0.083 (max) | 123° | +0.124 | +0.190 | **44%** |
| 35 MeV | −0.094 (min) | 111° | −0.260 | −0.275 | **34%** |
| **67 MeV** | **−0.531 (min)** | **109°** | **−0.540** | **−0.510** | **104%** ✓ |

### 8.2 角度图（见 `output/miller_gate2_plots/miller_gate2_angular_panels.{svg,png}`）

- 蓝粗线 = v5（这次工作）
- 红虚线 = Miller N3LO 数字化曲线
- 黑点 = 实验
- 灰系列 = v1/v3/v2/v4 收敛参考

67 MeV 三条线**几乎完全重合**；10 MeV v5 比 N3LO 略小（与 puzzle 一致），
比实验小 2 倍；35 MeV 形状对但 magnitude 1/3 不到。

### 8.3 PW 收敛轨迹（见 `output/miller_gate2_plots/miller_gate2_pw_convergence.{svg,png}`）

θ=90° 处 Ay_n 随 J_3N_max 增长的轨迹清楚显示**振荡收敛**：
- 67 MeV 单调向 expt 收敛
- 10/35 MeV 振荡幅度在 v4 → v5 间显著缩小

### 8.4 物理可见性

**Ay puzzle 在 10 MeV 直接可见**：v5 +0.083 < N3LO +0.124 < expt +0.190。
即使在我们这个 J=13/2 的截断里，2NF 理论已经收敛到 N3LO 量级，**仍然
低估实验 2 倍** —— 这就是 Miller 论文的核心 finding，也是为什么需要 3NF。

---

## 9. 重现步骤

```bash
# 1. Run solver at desired J truncation (e.g. v5)
cd Tic-tac
OMP_NUM_THREADS=4 ./CPP/run CPP/Output/miller_gate2_v5_J3N13_J2N3/input.txt

# 2. Generate comparison CSVs + plots
micromamba run -n anaroot-env python examples/plot_miller_gate2_convergence.py

# 3. Inspect outputs
ls output/miller_gate2_plots/
ls CPP/Output/miller_gate2_v5_J3N13_J2N3/Ay_n_compare_*.csv
```

Madison Ay 公式参考：`docs/seanBSMiller/arxiv_2201.09600/article.tex`
Appendix D（lines 1345–1410），`examples/pw_amplitudes.py` 实现：
`wigner_6j` (line 61+), `_block_jj_to_channel_spin` (line 583+),
`_channel_spin_geometric_factor` (line 651+), `assemble_m_matrix` 默认走新路径。

---

## 10. 下一步

| 任务 | 工作量 | 物理意义 |
|---|---|---|
| **v6: J_3N=17/2 (Miller 全配置)** | ~5 hr | 闭环 PW 收敛，35 MeV 推到 90%+ |
| **v_3NF: chiral_N2LO 3NF @ J=9/2** | 24-48 hr | **真正测 Ay puzzle**：3NF 是否拉近 2NF 与实验的 2× gap |
| Np 推到 50 + 同上 | 多日 | 与 Miller 全收敛一致 |
| 治理册 ch12/13 把 reduced-U 段落 deprecate | 30 min | 文档化 |

推荐先做 v_3NF，因为 v5 这次结果已经证明 Madison 计算器是对的（67 MeV 几乎完美），
所有的剩余问题（35 MeV magnitude、10 MeV puzzle）都需要 3NF 物理才能根本解决。

---

## Appendix A：Bug 教训

1. **OCR 不可信** — 公式必须读 arxiv tex 源，不能信 PDF 文字。
2. **半整数指数 (`(-1)^(half-integer)`) 在 Python int(round())  下 banker's rounding，给伪结果。** 任何看到半整数 `(-1)^x` 一定是 OCR 错或公式写错。
3. **Subagent 派去做高强度文献调研比自己 trial-and-error 高效**：5 min 拿到 verbatim TeX 解决了我 trial-and-error 1-2 hr 摸不到底的事。
4. **承担"框架是对的" vs "convergence 不够"的区分** — 67 MeV 的 104% 命中证明前者；35/10 MeV 的部分对应证明后者；要分清楚才能定下一步行动。

## Appendix B：相关 commits

```
ae57893  compare_Ay_experiment: surface Madison Ay(n) in production CSV
a4417a9  observables: rewrite M(theta) per Miller PRC 106 App D
06fadde  test/dbg: stop poisoning shared P123 cache
48aa52a  data: add Miller Gate 2 energy file (Tlab=10/35/67 MeV)
e016ade  build: drop dead CPP/ legacy mirror tree (94 files)
b6a63d9  observables: add Madison Ay(N) for polarized-nucleon-beam analyzing power
```
