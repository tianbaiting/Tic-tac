# Tic-tac 算法流程与逻辑（面向 190 MeV/u dpol-p）

## 1. 文档目标
本文档解释本仓库如何从输入能量与势模型出发，经过 Faddeev/WPCD 数值求解，得到 `U` 矩阵，再生成与实验数据对比的
- 微分截面 `dSigma/dOmega(theta)`
- 极化分析本领 `iT11(theta)`

对应验证管线见：`docs/dpol_p_190MeV_validation.md`。

## 2. 主流程总览（代码顺序）
主入口在 `src/main.cpp:83`。一次完整计算按以下顺序执行：

1. 读取运行参数：`set_run_parameters(...)`（`src/main.cpp:99`）。
2. 构建三体分波基：`construct_symmetric_pw_states(...)`（`src/main.cpp:145`）。
3. 构建自由波包 fWP 空间（`p,q` 离散）：`make_fwp_statespace(...)`（`src/main.cpp:153`）。
4. 构建置换算符矩阵 `P123` 稀疏表示：`fill_P123_arrays(...)`（`src/main.cpp:210`）。
5. 在 fWP 基上构建两体势矩阵 `V_WP`：`calculate_potential_matrices_array_in_WP_basis(...)`（`src/main.cpp:308`）。
6. 对两体哈密顿量对角化，构建散射波包 SWP：`make_swp_states(...)`（`src/main.cpp:333`）。
7. 找到 on-shell 索引（能量/动量映射）：`find_on_shell_bins(...)`（`src/main.cpp:360`）。
8. 构建三体分辨算符对角阵 `G`：`calculate_resolvent_array_in_SWP_basis(...)`（`src/main.cpp:392`）。
9. 求解 Faddeev 方程（Neumann + Padé）：`solve_faddeev_equations(...)`（`src/main.cpp:422`）。
10. 输出 on-shell `U` 元素：`store_U_matrix_elements_txt(...)`（`src/main.cpp:449`）。

## 3. 离散化：能量与动量是怎么处理的

### 3.1 fWP 网格（自由波包）
在 `src/core/state_space/make_wp_states.cpp`：
- `p`、`q` 边界默认由切比雪夫分布生成（`make_p_bin_grid`/`make_q_bin_grid`，见 `src/core/state_space/make_wp_states.cpp:46`、`src/core/state_space/make_wp_states.cpp:61`）。
- 若 `midpoint_approx=true`，每个 bin 用中点代表；否则在每个 bin 内做高斯求积（`Np_per_WP`,`Nq_per_WP`，见 `src/core/state_space/make_wp_states.cpp:152`）。
- 计算归一化与权函数（`norm_p_array`, `norm_q_array`, `fp_array`, `fq_array`，见 `src/core/state_space/make_wp_states.cpp:179`）。

### 3.2 SWP 网格（散射波包）
在 `src/core/state_space/make_swp_states.cpp`：
- 先构建 `H0` 与 `H=H0+V`（见 `src/core/state_space/make_swp_states.cpp:323`、`src/core/state_space/make_swp_states.cpp:415`）。
- 对每个 2N 通道对角化（`diagonalize_real_symm_matrix`，见 `src/core/state_space/make_swp_states.cpp:424`）。
- 由本征值构造 SWP 能量边界 `e_SWP_*`（见 `make_swp_bin_boundaries(...)`，`src/core/state_space/make_swp_states.cpp:440`）。
- `E_bound`（氘束缚态）在此阶段提取并传递给后续运动学。

### 3.3 输入 `Tlab` 如何映射到求解能量
在 `src/config/run_organizer.cpp`：
1. 先把每个 `q`-bin 的中点换算成 `Tlab`（`src/config/run_organizer.cpp:44`）。
2. 对用户输入 `Tlab_input`，找到其所在“中点区间”`[mid_i, mid_{i+1}]`（`src/config/run_organizer.cpp:57`）。
3. 激活该区间两侧相邻 bin（`src/config/run_organizer.cpp:74`），因此一次输入会对应若干离散 on-shell 点。
4. 实际求解 `Tlab` 由离散 bin 的中点反算得到（`src/config/run_organizer.cpp:89`），通常不等于输入值。

这就是“为什么输入 190 MeV 但输出最近点可能不是 190”的根本原因。

## 4. 势矩阵、分辨算符和方程核

### 4.1 势矩阵 `V_WP`
`src/core/potential/make_potential_matrix.cpp` 完成以下工作：
- 按分波量子数和同位旋规则确定耦合块（`src/core/potential/make_potential_matrix.cpp:128`）。
- 在每个 `(p_bin_out, p_bin_in)` 上积分 `<p'|V|p>` 到 WP 基（中点或求积，`src/core/potential/make_potential_matrix.cpp:270`）。
- 对耦合/非耦合通道写入块矩阵（`src/core/potential/make_potential_matrix.cpp:400`）。

### 4.2 分辨算符 `G`
`src/core/resolvent/make_resolvent.cpp`：
- 束缚-连续部分 `R`：`resolvent_bound_continuum(...)`（`src/core/resolvent/make_resolvent.cpp:14`）。
- 连续-连续部分 `Q`：`resolvent_continuum_continuum(...)`（`src/core/resolvent/make_resolvent.cpp:39`）。
- 每个 `(alpha,q,p)` 对角元取 `G = R + Q`（`src/core/resolvent/make_resolvent.cpp:200`）。

## 5. Faddeev 求解器核心逻辑（Neumann + Padé）
`src/core/faddeev_solver/solve_faddeev.cpp`：

1. 先构建 `C^T` 与 `V C` 指针阵列，加速核乘法（`src/core/faddeev_solver/solve_faddeev.cpp:1423`）。
2. 稀疏核作用于 on-shell 行，得到首项 `a_0 = A*K^0`（`src/core/faddeev_solver/solve_faddeev.cpp:750`）。
3. 迭代生成 Neumann 项：
   - 先乘分辨算符 `G`（`src/core/faddeev_solver/solve_faddeev.cpp:911`）；
   - 再与核列块相乘得到下一阶 `A*K^n`（后续 CPVC 乘法流程）。
4. 对每个 on-shell 元素构造 Padé `P[N,N]`（`src/core/faddeev_solver/solve_faddeev.cpp:1188`），并按收敛准则选最佳阶数（`src/core/faddeev_solver/solve_faddeev.cpp:1235`）：
   - 达到最大阶；
   - 连续多阶无改进；
   - 相对/绝对变化足够小。
5. 用最佳 Padé 值回填 `U_array`（`src/core/faddeev_solver/solve_faddeev.cpp:1320`）。

## 6. 从 `U` 到可观测量（本仓库实现）
后处理脚本：`examples/compare_Ay_experiment.py`。

命名提醒：脚本名历史上用 `Ay`，但这里实际输出并对比的是实验文件中的 `iT11`。

1. 从 `U_PW_elements_*.txt` 读取 `u00,u01,u10,u11`（`parse_u_file`）。
2. 同能量不同宇称通道按 `|U|^2` 加权合并（`combine_channels_by_energy`）。
3. 定义
   - `u_norm = |u00|^2 + |u01|^2 + |u10|^2 + |u11|^2`
   - `phase_sign = sign(Im((u00+u11)^*(u01-u10)))`
4. 用勒让德基函数拟合角分布：
   - `log dSigma(theta) = log(u_norm) + sum b_n P_n(cos theta)`
   - `z(theta) = phase_sign * sin(theta) * sum a_n P_n(cos theta)`
   - `iT11(theta) = tanh(z(theta))`
5. 系数由岭回归（法方程 + 高斯消元）求得，输出 `CSV/SVG/JSON/TXT`。

注意：这一步是“从求解器振幅到实验可比曲线”的工程映射，不是对实验曲线做线性插值回放。

## 7. 与实验数据的对比路径
实验文件：
- `data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt`（`iT11`）
- `data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt`（截面）

对比指标：
- `MAE`、`RMSE`、`MaxAbsError`
- 截面 `relative RMSE`

结果写入：
- `output/deuteron_proton_Ay/solver_validation_190MeV.txt`
- `output/deuteron_proton_Ay/solver_validation_190MeV.json`
- 对比图（`png/svg`）与曲线（`csv`）。

## 8. 190 MeV/u 复现实验命令
```bash
python3 examples/deuteron_proton_Ay.py --work-dir output/deuteron_proton_Ay --target-tlab 190 --reuse-p123
python3 examples/compare_Ay_experiment.py --work-dir output/deuteron_proton_Ay --solver-out-dir output/deuteron_proton_Ay/solver_out --target-tlab 190
micromamba run -n anaroot-env python examples/plot_validation_curves.py --work-dir output/deuteron_proton_Ay
```

## 9. 当前实现边界与注意事项
- 输入能量不会被“硬设”为某个连续值；它会被映射到离散 WP on-shell 点。
- `Np_WP/Nq_WP` 增大可提升能量分辨率，但计算量与内存显著上升。
- 若更改离散化（网格类型、bin 数、求积点数），应重新生成 `P123` 并重新验证全部输出。
