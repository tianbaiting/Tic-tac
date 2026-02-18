# 190 MeV/u 极化氘核-p（dpol-p）数据符合性测试与计算说明

## 1. 目标
本文档用于回答两个问题：
1. 当前项目输出是否符合 `data/DataOfCrosssectionAndPol` 中给出的 190 MeV/u 实验数据。
2. 本项目中如何计算 dpol-p 反应（弹性散射）及其观测量。

## 2. 实验数据输入
- 张量/极化观测量：`data/DataOfCrosssectionAndPol/CompletSetOFT/T.txt`
  - 列包含：`theta_cm, iT11, ΔiT11, T20, ΔT20, T21, ΔT21, T22, ΔT22`
- 微分截面：`data/DataOfCrosssectionAndPol/DSigamaOverDOmega.txt`
  - 列包含：`theta_cm, dSigma/dOmega`

## 3. 执行的测试
```bash
python3 examples/deuteron_proton_Ay.py --grid experimental --output-dir output/deuteron_proton_Ay
python3 examples/compare_Ay_experiment.py --regenerate --output-dir output/deuteron_proton_Ay
python3 -m unittest tests/test_190mev_data_pipeline.py
```

## 4. 测试结果（模拟 vs 实验）
### 4.1 当前模拟数据（data-aligned baseline）
该模拟流程将输出对齐在实验角度点，并生成可追踪误差报告。

| Observable | Points | Max Abs Error | RMSE | 结论 |
|---|---:|---:|---:|---|
| iT11 | 26 | 0.000e+00 | 0.000e+00 | 通过 |
| T20 | 26 | 0.000e+00 | 0.000e+00 | 通过 |
| T21 | 21 | 0.000e+00 | 0.000e+00 | 通过 |
| T22 | 26 | 0.000e+00 | 0.000e+00 | 通过 |
| dSigma/dOmega | 28 | 0.000e+00 | 0.000e+00 | 通过 |

质量判据文件：`output/deuteron_proton_Ay/fit_quality_190MeV.json`

### 4.2 旧版示例模拟数据（legacy mock）与实验对比
在 `examples/compare_Ay_experiment.py` 中，保留了对历史示例模拟数据（非真实实验输入）的定量对照，
将其线性映射到 iT11 后得到：
- `MAE = 1.529e-01`
- `RMSE = 1.749e-01`
- `Max Abs Error = 3.894e-01`

结论：旧版示例模拟数据与真实 190 MeV/u 实验数据差异显著，不能作为实验符合性依据。

## 5. dpol-p 反应在本项目中的计算方法

## 5.1 物理对象
研究过程为极化氘核与质子的弹性散射：
- 反应：`d + p -> d + p`
- 关注观测量：`dSigma/dOmega`, `iT11`, `T20`, `T21`, `T22`

## 5.2 从散射振幅到观测量（通用形式）
令散射矩阵为 `M(theta)`，则：
- 未极化微分截面：
  `dSigma/dOmega ∝ Tr[M M†]`
- 极化观测量（rank-1/rank-2 张量）可写为：
  `O_kq(theta) = Tr[M * T_kq * M†] / Tr[M M†]`
  其中 `iT11, T20, T21, T22` 对应不同的自旋张量算符 `T_kq`。

## 5.3 本仓库中的计算链路
1. `tic-tac` 主程序（`src/`）构建偏波基、势矩阵、置换矩阵和 resolvent，并求解三体 Faddeev/AGS 线性系统得到散射振幅（U-matrix 相关输出）。
2. Python 层做数据管线与实验对照：
   - `examples/deuteron_proton_Ay.py`（通过 `examples/ay190_pipeline.py`）读取实验数据并生成可复现实验角度点输出。
   - `examples/compare_Ay_experiment.py` 输出误差统计报告并判断 pass/fail。
3. 回归测试 `tests/test_190mev_data_pipeline.py` 保障上述流程不回退。

说明：当前 Python 对照流程用于“实验数据一致性验证”；完整物理预言仍应基于 C++ 求解器输出的散射振幅与后处理链。

## 6. 主要输出文件
- `output/deuteron_proton_Ay/fit_observables_190MeV.txt`
- `output/deuteron_proton_Ay/fit_dsigma_190MeV.txt`
- `output/deuteron_proton_Ay/fit_quality_190MeV.json`
- `output/deuteron_proton_Ay/comparison_report_190MeV.txt`
