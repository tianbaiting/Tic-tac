# 190 MeV dpol-p：Tic-tac 求解器与实验数据对比

## 1. 目标
用 Tic-tac 求解器真实输出（`U_PW_elements`）验证 190 MeV dpol-p 数据，不再用实验数据插值去“构造通过”。

## 2. 实验数据
- 张量观测量：`data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt`
  - 用于对比 `iT11`
- 微分截面：`data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt`

## 3. 执行命令
```bash
python3 examples/deuteron_proton_Ay.py --work-dir output/deuteron_proton_Ay --reuse-p123
python3 examples/compare_Ay_experiment.py --work-dir output/deuteron_proton_Ay --solver-out-dir output/deuteron_proton_Ay/solver_out
```

说明：
- `deuteron_proton_Ay.py` 会调用 `CPP/run`（Tic-tac 求解器），产出 `U_PW_elements_*.txt`。
- `compare_Ay_experiment.py` 只从求解器输出计算代理量并对比实验，不走实验曲线插值。

## 4. 当前结果（最近一次实测）
- 状态：`FAIL`
- 最接近目标能量：`Tlab = 161.145 MeV`（目标 `135.6 MeV`，差值 `25.545 MeV`）
- `Ay` 代理量与 `iT11` 对比：
  - `MAE = 0.162860`
  - `RMSE = 0.181008`
  - `MaxAbsError = 0.342954`
- 截面代理量（标量）对 `dSigma/dOmega(θ)` 角分布：
  - `shape RMSE = 1.497132`

完整报告：
- `output/deuteron_proton_Ay/solver_validation_190MeV.txt`
- `output/deuteron_proton_Ay/solver_validation_190MeV.json`

## 5. 目前如何从求解器输出计算对比量
在 `examples/compare_Ay_experiment.py` 中，对每个能量点读取：`U00, U01, U10, U11`，定义
- `f_no_flip = U00 + U11`
- `f_flip = U01 + U10`
- `Ay_proxy = Im(f_no_flip* * f_flip) / (|f_no_flip|^2 + |f_flip|^2)`
- `dSigma_proxy = |f_no_flip|^2 + |f_flip|^2`

然后按宇称通道用 `dSigma_proxy` 加权合并为单个能量点，再与实验数据做误差统计。

## 6. 结论
- 现在的流程已经是“求解器驱动验证”，不是实验插值回放。
- 以当前参数（`two_J_3N_max=1, J_2N_max=2, Np=Nq=20, LO_internal`）与当前代理定义，尚不能匹配 190 MeV dpol-p 实验数据。
- 下一步应提高分波截断、网格分辨率和势模型精度，并实现完整的角分布观测量重建，而不是标量代理。
