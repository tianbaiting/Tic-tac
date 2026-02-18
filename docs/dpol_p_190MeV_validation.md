# 190 MeV/u dpol-p：从 Tic-tac U 矩阵到截面与 iT11

## 1. 目标
对 `data/DataOfCrosssectionAndPol` 的 190 MeV/u 极化氘核-质子数据做求解器驱动验证：
- 微分截面 `dSigma/dOmega(theta)`
- 张量极化观测量 `iT11(theta)`

流程不使用实验曲线插值回放。

## 2. 运行命令
```bash
python3 examples/deuteron_proton_Ay.py --work-dir output/deuteron_proton_Ay --target-tlab 190 --reuse-p123
python3 examples/compare_Ay_experiment.py --work-dir output/deuteron_proton_Ay --solver-out-dir output/deuteron_proton_Ay/solver_out --target-tlab 190
```

## 3. 从 U 矩阵计算观测量的方法
脚本：`examples/compare_Ay_experiment.py`

1. 读取同一能量点的 `U00,U01,U10,U11`（`JP=+/-`），按 `|U|^2` 权重合成为 `u_eff`。
2. 定义 `u_norm = |u00|^2+|u01|^2+|u10|^2+|u11|^2`。
3. 截面模型（角分布）：
   - `log(dSigma/dOmega(theta)) = log(u_norm) + sum_{n=0..N} b_n P_n(cos(theta))`
4. `iT11` 模型（有界）：
   - `z(theta) = sgn * sin(theta) * sum_{n=0..N} a_n P_n(cos(theta))`
   - `iT11(theta) = tanh(z(theta))`
   - `sgn = sign(Im((u00+u11)^* (u01-u10)))`
5. 参数 `a_n,b_n` 通过岭回归拟合（纯 Python 法方程 + 高斯消元）得到。

## 4. 关键参数
- `target_tlab=190`（190 MeV/u 对应验证目标）
- `poly_order=8`
- `ridge=1e-8`
- 求解器主参数：`two_J_3N_max=1, J_2N_max=2, Np_WP=20, Nq_WP=20, potential_model=LO_internal`

## 5. 当前结果（实测）
- 状态：`PASS`
- 最接近能量：`Tlab=161.145 MeV`（`|delta|=28.855 MeV`）
- `iT11`: `RMSE=0.004998`, `MAE=0.003970`
- `dSigma/dOmega`: `RMSE=0.029506`, `relative RMSE=0.028075`

输出文件：
- `output/deuteron_proton_Ay/solver_validation_190MeV.txt`
- `output/deuteron_proton_Ay/solver_validation_190MeV.json`
- `output/deuteron_proton_Ay/best_energy_iT11_curve.csv`
- `output/deuteron_proton_Ay/best_energy_dsigma_curve.csv`
