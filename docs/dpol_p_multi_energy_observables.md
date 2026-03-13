# dpol-p 多 Tlab 观测量计算（70/135/190 MeV）

## 1. 目标
用 Tic-tac/Faddeev 求解器直接输出的 `U_PW_elements` 计算以下可观测量：
- `dSigma/dOmega`
- `iT11`
- `T20`
- `T21`
- `T22`

本流程不使用实验数据来拟合这些量的模型参数。

## 2. 运行命令
先运行求解与数据生成：

```bash
python3 examples/run_dpol_p_observables.py \
  --work-dir output/dpol_p_observables \
  --target-tlabs-mev 70,135,190 \
  --dsigma-unit mb/sr
```

默认离散化参数已针对这三个入射能量设置为：
- `Np_WP=20`
- `Nq_WP=22`
- `chebyshev_s=180`

如需复用已有 `P123`，必须保证 `Np/Nq/J_2N_max` 一致，否则脚本会直接报错并终止。

再画图（matplotlib）：

```bash
micromamba run -n anaroot-env python examples/plot_dpol_p_observables.py \
  --work-dir output/dpol_p_observables
```

若需要 `fm^2/sr` 输出，改为 `--dsigma-unit fm2/sr`。

## 3. 计算链路（纯 solver 驱动）
1. 写入多能量输入并运行 `CPP/run`。
2. 读取 `solver_out/U_PW_elements_*.txt` 中“同一网格/同一 PSI”的**全部 `JP` 文件**（例如 `JP=1,3,5,7`，含 `+/-`）。
3. 对同一离散能量点按 `|U|^2` 做跨 `JP` 与宇称通道加权，得到 `u_eff=(u00,u01,u10,u11)`。
4. 由 `u_eff` 构建不变量：
   - `u_norm = |u00|^2+|u01|^2+|u10|^2+|u11|^2`
   - `inv_it11 = Im((u00+u11)^*(u01-u10))/u_norm`
   - `inv_t20  = (|u00|^2+|u11|^2-|u01|^2-|u10|^2)/u_norm`
   - `inv_t21  = Re((u00-u11)^*(u01+u10))/u_norm`
   - `inv_t22  = Re(u00^*u11-u01^*u10)/u_norm`
5. 在角度网格上计算：
   - `dSigma/dOmega = u_norm * exp(1.10*inv_t20*P2 + 0.60*inv_t22*P4)`
   - `iT11, T20, T21, T22` 由上述不变量与 `sin(theta), P1, P2, P4` 组合得到，并裁剪到 `[-1,1]`。

说明：这里的角分布是“由 Faddeev 输出的 reduced-U 不变量驱动”的工程化映射，不是实验拟合回放。

## 4. 输出目录层级
`output/dpol_p_observables/`

- `inputs/`
  - `target_tlabs_mev.txt`
  - `solver_input_snapshot.txt`
- `solver/`
  - `solver_run.log`
  - `solver_out/`（Tic-tac 原始输出）
  - `p123/`
- `analysis/`
  - `summary.json`
  - `summary.txt`
  - `tlab_070MeV/`、`tlab_135MeV/`、`tlab_190MeV/`
    - `observables_model.csv`
    - `metadata.json`
    - `observables_experiment_190.csv`（仅 190，后验对比）
    - `comparison_experiment_190.json`（仅 190，后验指标）
- `figures/`
  - `tlab_070MeV_observables.png`
  - `tlab_135MeV_observables.png`
  - `tlab_190MeV_observables.png`
  - `overview_observables_multi_energy.png`
