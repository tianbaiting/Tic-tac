# Tic-tac 可调参数与调参指南

本文档对应 `./CPP/run`（与 `src/` 同步实现），用于说明哪些参数可调、各自作用和推荐调法。

## 1. 参数输入方式（重要）

支持两种输入：
- 输入文件：`./CPP/run CPP/Input/input.txt`
- 命令行 `key=value`：`./CPP/run Np_WP=30 Nq_WP=30 energy_input_file=CPP/Input/lab_energies_190MeV.txt`

可混用，且按“从左到右”覆盖，后出现的值会覆盖前面的值，例如：
`./CPP/run CPP/Input/input.txt Nq_WP=40`

布尔参数只能写 `true/false`（小写）。

## 2. 可调参数总览

### 2.1 离散化与基空间（精度/耗时主开关）

| 参数 | 默认值 | 作用 | 调参建议 |
|---|---:|---|---|
| `two_J_3N_max` | `1` | 三体总角动量截断（按 `2J` 存） | 必须是正奇数；做收敛时用 `1 -> 3 -> 5` |
| `J_2N_max` | `3` | 两体分波角动量截断 | 高能量通常需更大；逐步增加检查收敛 |
| `Np_WP` | `50` | `p` 波包 bin 数 | 增大可提升势矩阵离散精度，耗时/内存上升 |
| `Nq_WP` | `50` | `q` 波包 bin 数与 on-shell 能量分辨率 | 想让目标 `Tlab` 更贴近离散点，优先增大这个参数 |
| `Np_per_WP` | `8` | 每个 `p` bin 求积点数 | 仅在 `midpoint_approx=false` 时有效 |
| `Nq_per_WP` | `8` | 每个 `q` bin 求积点数 | 同上 |
| `Nphi` | `48` | `P123` 计算中的角积分点数 | `P123` 噪声大/不稳定时增大 |
| `Nx` | `48` | 几何函数角积分点数 | 与 `Nphi` 一起调，先小步加 |
| `chebyshev_s` | `100` | Chebyshev 网格尺度 | 增大可覆盖更高动量/能量区间 |
| `chebyshev_t` | `1` | Chebyshev 稀疏度指数 | 通常保持 `1`，仅在收敛研究中改动 |
| `p_grid_type` | `chebyshev` | `p` 网格类型 | `chebyshev` 或 `custom` |
| `p_grid_filename` | `""` | 自定义 `p` 网格文件 | `custom` 时必须提供，文件应有 `Np_WP+1` 行单调边界 |
| `q_grid_type` | `chebyshev` | `q` 网格类型 | `chebyshev` 或 `custom` |
| `q_grid_filename` | `""` | 自定义 `q` 网格文件 | `custom` 时必须提供，文件应有 `Nq_WP+1` 行单调边界 |
| `midpoint_approx` | `false` | 用 bin 中点近似替代 bin 内求积 | `true` 更快但精度更低，只建议调试 |

### 2.2 物理模型与参数扫描

| 参数 | 默认值 | 作用 | 调参建议 |
|---|---:|---|---|
| `potential_model` | `LO_internal` | 两体势模型 | 可选 `LO_internal/N2LOopt/Idaho_N3LO/nijmegen/malfliet_tjon` |
| `tensor_force` | `true` | 是否保留张量耦合 | 物理计算建议保持 `true` |
| `isospin_breaking_1S0` | `true` | 是否引入 `1S0` 电荷相关破缺 | 与实验比较通常保持 `true` |
| `parameter_walk` | `false` | 是否进行势参数扫描 | 打开后需同时正确设置下面三个参数 |
| `parameter_file` | `none` | 参数扫描文件 | 每行一个参数集，空格分隔 |
| `PSI_start` | `-1` | 扫描起始索引（含） | 0 基索引 |
| `PSI_end` | `-1` | 扫描终止索引（不含） | 需满足 `PSI_end > PSI_start` |

### 2.3 求解器控制与并行

| 参数 | 默认值 | 作用 | 调参建议 |
|---|---:|---|---|
| `P123_omp_num_threads` | `omp_get_max_threads()` | `P123` OpenMP 线程数 | 实际会被截断到 `Nq_WP` |
| `calculate_and_store_P123` | `true` | 计算并写出 `P123_sparse*.h5` | 首次跑新网格必须 `true` |
| `P123_recovery` | `false` | 从 TFC 子文件恢复 `P123` | 仅在中断恢复/拼接场景使用 |
| `P123_folder` | `Output` | `P123` 读写目录 | 固定好后可复用，换网格要重算 |
| `solve_faddeev` | `true` | 是否求解 Faddeev 方程 | 若为 `false`，仍需 `calculate_and_store_P123=true` 才有工作可做 |
| `solve_dense` | `false` | 是否用 dense LAPACK 解法 | 很慢，仅调试用 |
| `include_breakup_channels` | `false` | 是否计算 breakup 振幅 | 成本很高，按需打开 |
| `production_run` | `true` | 生产/调试模式 | 验证物理结果必须 `true` |
| `parallel_run` | `false` | 按 3N 通道并行分任务 | `true` 时必须给 `channel_idx` |
| `channel_idx` | `-1` | 指定单通道索引（0 基） | 与 `parallel_run=true` 配套 |

### 2.4 输入输出与历史参数

| 参数 | 默认值 | 作用 | 说明 |
|---|---:|---|---|
| `energy_input_file` | `lab_energies.txt` | 入射能量列表文件（`.txt`） | 每行一个 `Tlab` |
| `output_folder` | `Output` | 结果输出目录 | 保存 `U_PW_elements*` 等 |
| `subfolder` | `Output` | 历史字段 | 当前主流程未使用，可忽略 |

## 3. 推荐调参顺序

1. 先定物理开关：`potential_model`, `tensor_force`, `isospin_breaking_1S0`。  
2. 再调能量离散：优先调 `Nq_WP`，再调 `chebyshev_s`；若要精确控制离散点，用 `q_grid_type=custom` + `q_grid_filename`。  
3. 再做角动量与积分收敛：`two_J_3N_max`, `J_2N_max`, `Nphi`, `Nx`。  
4. 最后平衡耗时：`Np_WP/Nq_WP`、`Np_per_WP/Nq_per_WP`、线程数。  
5. 网格或角动量截断改变后，不要复用旧 `P123`，需重新生成。  

## 4. 常见失败与对应参数

- `Cannot have even two_J_3N_max!`：`two_J_3N_max` 必须为正奇数。  
- `Parallel run ... but no channel_idx`：`parallel_run=true` 时要设置 `channel_idx>=0`。  
- `Both solve_faddeev=false and calculate_and_store_P123=false`：至少打开一个。  
- `Program mismatch with number of wave-packets in input WP-boundaries file`：自定义网格文件行数不对（应为 `N_WP+1`）。  
- 参数扫描无输出：检查 `parameter_walk=true` 时 `parameter_file/PSI_start/PSI_end` 是否合理。  
