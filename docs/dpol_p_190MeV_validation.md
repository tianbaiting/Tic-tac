# Tlab = 190 MeV dpol-p：从 Tic-tac U 矩阵到截面与 iT11

详细算法流程（离散化、Faddeev 求解、U 到可观测量映射）请看：
- `docs/algorithm_flow_and_logic.md`

## 1. 目标
对 `data/DataOfCrosssectionAndPol` 的极化氘核-质子基准数据做求解器驱动验证：
- 微分截面 `dSigma/dOmega(theta)`
- 张量极化观测量 `iT11(theta)`

流程不使用实验曲线插值回放。

说明：外部资料常把该实验基准写成 `190 MeV/u`。本仓库与 C++ core 保持一致，solver 输入统一使用 `Tlab [MeV]`，维护点为 `Tlab = 190 MeV`。

## 2. 运行命令
```bash
python3 examples/deuteron_proton_Ay.py --work-dir output/deuteron_proton_Ay --target-tlab-mev 190 --reuse-p123
python3 examples/compare_Ay_experiment.py --work-dir output/deuteron_proton_Ay --solver-out-dir output/deuteron_proton_Ay/solver_out --target-tlab-mev 190
```

## 3. Experiment Data 与 Faddeev Simulation 的计算链路
核心脚本：`examples/compare_Ay_experiment.py`

### 3.1 Experiment Data（实验数据）如何得到
1. `iT11` 来自 `data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt`。
2. `dSigma/dOmega` 来自 `data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt`。
3. 读取时仅做最小清洗：
   - 跳过表头和注释行；
   - `T.txt` 中 `null` 行丢弃；
   - 不做插值、不做平滑、不做人为重采样。
4. `dSigma/dOmega` 单位默认按文件头自动识别（当前实验文件为 `mb/sr`），并可在脚本中统一换算到：
   - `mb/sr`（默认）
   - `fm^2/sr`（`--dsigma-unit fm2/sr`）
5. 仅做物理单位换算，不做实验拟合。
6. 结果直接作为对照基准（ground truth）参与误差评估。

### 3.2 Faddeev Simulation（模拟数据）如何得到
1. 先由 `examples/deuteron_proton_Ay.py` 调用 `CPP/run`，求解 Faddeev 方程，得到 `U_PW_elements_*.txt`。
2. 在同一 `Tlab` 点读取最新一组参数下全部 `U_PW_elements`（覆盖 `JP=1,3,5,...` 与 `+/-` 宇称）。
3. 按每条通道的 `|U|^2` 权重合并，得到有效振幅 `u_eff=(u00,u01,u10,u11)`。
4. 定义振幅尺度：
   - `u_norm = |u00|^2+|u01|^2+|u10|^2+|u11|^2`
5. 由 `u_eff` 构造 reduced-U 不变量：
   - `u_norm = |u00|^2+|u01|^2+|u10|^2+|u11|^2`
   - `inv_it11 = Im((u00+u11)^*(u01-u10))/u_norm`
   - `inv_t20  = (|u00|^2+|u11|^2-|u01|^2-|u10|^2)/u_norm`
   - `inv_t22  = Re(u00^*u11-u01^*u10)/u_norm`
6. 在实验角度点上用固定映射生成预测曲线（无实验拟合）：
   - `dSigma/dOmega = u_norm * exp(1.10*inv_t20*P2 + 0.60*inv_t22*P4)`
   - `iT11 = clamp(1.20*inv_it11*sin(theta)*(-1.20*P2 - 0.20*P1), -1, 1)`
7. 模型内部 `dSigma` 基准单位按 `mb/sr` 记账，再根据 `--dsigma-unit` 输出到目标单位。
8. 实验数据仅用于残差评估（MAE/RMSE/MaxAbs/relative RMSE），不参与参数训练。

### 3.3 对比与误差指标
1. 逐角度比较 `Experiment` 与 `Faddeev Simulation` 曲线。
2. 计算：
   - `MAE = mean(|pred-exp|)`
   - `RMSE = sqrt(mean((pred-exp)^2))`
   - `MaxAbsError = max(|pred-exp|)`
   - `dSigma relative RMSE = sqrt(mean(((pred-exp)/exp)^2))`

## 4. 关键参数
- `target_tlab_mev=190`
- 求解器主参数：`two_J_3N_max=1, J_2N_max=2, Np_WP=20, Nq_WP=20, potential_model=LO_internal`
- 评估阈值：`ay_rmse_pass=0.02`, `dsigma_rel_rmse_pass=0.05`, `energy_delta_pass=40`

## 5. 当前结果（以本地最新运行为准）
不要在文档中手写固定数值，直接读取本次运行输出：
- `output/deuteron_proton_Ay/solver_validation_tlab_190MeV.json`
- `output/deuteron_proton_Ay/solver_validation_tlab_190MeV.txt`

快速查看关键字段：

```bash
python3 - <<'PY'
import json
from pathlib import Path
p = Path('output/deuteron_proton_Ay/solver_validation_tlab_190MeV.json')
if not p.exists():
    raise SystemExit('validation JSON not found')
j = json.loads(p.read_text())
best = j.get('best_energy', {})
print('target_tlab_mev=', j.get('target_tlab_mev'))
print('best_tlab_mev=', best.get('tlab_mev'), 'delta=', best.get('abs_delta_tlab_mev'))
print('iT11 RMSE=', best.get('it11_rmse'), 'MAE=', best.get('it11_mae'))
print('dSigma RMSE=', best.get('dsigma_rmse'), 'relative RMSE=', best.get('dsigma_rel_rmse'))
print('status=', j.get('status'))
print('dSigma unit=', j.get('units', {}).get('dsigma_output'))
PY
```

输出文件：
- `output/deuteron_proton_Ay/solver_validation_tlab_190MeV.txt`
- `output/deuteron_proton_Ay/solver_validation_tlab_190MeV.json`
- `output/deuteron_proton_Ay/best_energy_iT11_curve.csv`
- `output/deuteron_proton_Ay/best_energy_dsigma_curve.csv`
- `output/deuteron_proton_Ay/best_energy_iT11_comparison.svg`
- `output/deuteron_proton_Ay/best_energy_dsigma_comparison.svg`

使用 `matplotlib` 生成 PNG 对比图：

```bash
micromamba run -n anaroot-env python examples/plot_validation_curves.py --work-dir output/deuteron_proton_Ay
```

生成：
- `output/deuteron_proton_Ay/best_energy_iT11_comparison.png`
- `output/deuteron_proton_Ay/best_energy_dsigma_comparison.png`
- `output/deuteron_proton_Ay/best_energy_comparison.png`
