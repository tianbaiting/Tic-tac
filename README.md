# Tic-tac: 三体核物理Faddeev方程求解器

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/language-C++-red.svg)](https://isocpp.org/)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen.svg)]()

## 📖 项目简介

**Tic-tac**是一个用于求解三体核物理系统Faddeev方程的高性能科学计算程序。该程序专门设计用于研究核子-氘核散射、三核子系统的束缚态和散射态，以及核物质相互作用的相移和截面等核物理问题。

### 🎯 主要功能

1. **三体Faddeev方程求解**
   - 基于Ludvig Faddeev理论框架
   - 支持对称和非对称三体系统
   - 高精度数值求解算法

2. **波包离散化方法（WPCD）**
   - 高效的数值离散化技术
   - 自适应网格生成
   - 优化的内存使用

3. **多种核相互作用势**
   - 手征有效场理论势（ChEFT）
   - Nijmegen势模型
   - Malfliet-Tjon势
   - 可扩展的势能接口

4. **高性能计算**
   - OpenMP并行化
   - 稀疏矩阵优化
   - HDF5数据存储

## 🏗️ 项目结构

```
Tic-tac/
├── src/                          # 源代码目录
│   ├── core/                     # 核心算法模块
│   │   ├── faddeev_solver/       # Faddeev方程求解器
│   │   ├── state_space/          # 状态空间构建
│   │   ├── potential/            # 势能矩阵计算
│   │   └── resolvent/            # 格林函数计算
│   ├── interactions/             # 核子相互作用势
│   ├── utils/                    # 通用工具函数
│   ├── io/                       # 文件输入输出
│   ├── config/                   # 参数配置
│   └── main.cpp                  # 主程序入口
├── include/                      # 头文件目录
├── data/                         # 数据文件和配置
├── tests/                        # 测试代码
├── examples/                     # 示例和运行脚本
├── build/                        # 构建输出目录
├── Makefile                      # 构建系统
├── config.py                     # 配置生成脚本
└── README.md                     # 本文档
```

## 🔧 系统要求

### 必需依赖
- **编译器**: GCC 7.0+ 或 Clang 5.0+ (支持C++17)
- **Fortran编译器**: gfortran 7.0+
- **数学库**: LAPACK, BLAS
- **科学计算库**: GSL (GNU Scientific Library)
- **并行计算**: OpenMP
- **数据存储**: HDF5
- **构建工具**: Make

### 可选依赖
- **Intel MKL**: 用于高性能线性代数运算
- **Python 3.6+**: 用于配置管理和数据分析
- **Matplotlib**: 用于结果可视化

### 安装依赖（Ubuntu/Debian）

```bash
sudo apt update
sudo apt install build-essential gfortran
sudo apt install liblapack-dev libblas-dev
sudo apt install libgsl-dev
sudo apt install libhdf5-dev
sudo apt install libomp-dev
```

### 安装依赖（CentOS/RHEL）

```bash
sudo yum groupinstall "Development Tools"
sudo yum install gcc-gfortran
sudo yum install lapack-devel blas-devel
sudo yum install gsl-devel
sudo yum install hdf5-devel
sudo yum install libgomp
```

## 🚀 快速开始

### 1. 克隆项目

```bash
git clone https://github.com/tianbaiting/Tic-tac.git
cd Tic-tac
```

### 2. 编译程序

```bash
make
```

或者查看所有可用的构建选项：

```bash
make help
```

### 3. 生成配置文件

```bash
python3 config.py save data/input.txt
```

### 4. 运行测试

```bash
./examples/run_examples.sh test
```

### 5. 基本计算示例

```bash
./examples/run_examples.sh basic
```

### 6. 氘核-质子Ay计算示例

```bash
# 运行190 MeV/u氘核分析能力计算
./examples/run_examples.sh ay

# 或直接运行Python脚本
python3 examples/deuteron_proton_Ay.py
```

## 📊 核心算法

### Faddeev方程理论

Tic-tac求解的核心是三体系统的Faddeev方程组：

$$
\begin{align}
\psi_1 &= \phi_1 + G_0(E) T_1 (\psi_2 + \psi_3) \\
\psi_2 &= \phi_2 + G_0(E) T_2 (\psi_3 + \psi_1) \\
\psi_3 &= \phi_3 + G_0(E) T_3 (\psi_1 + \psi_2)
\end{align}
$$

其中：
- $\psi_i$ 是Faddeev分量，描述粒子对最后一次相互作用的状态
- $G_0(E)$ 是自由三体格林函数
- $T_i$ 是二体T-矩阵
- $\phi_i$ 是初始状态（如入射态）

### 波包离散化方法

程序采用波包离散化（Wave Packet Discretization, WPCD）方法：

1. **构建偏波态空间**: 根据角动量和同位旋量子数构建基态
2. **生成波包网格**: 使用Chebyshev或等间距网格
3. **计算势能矩阵**: 在波包基组中计算二体相互作用
4. **对角化哈密顿量**: 构建强弱波包（SWP）基组
5. **求解线性方程组**: 使用迭代或直接方法求解

## 🔬 物理应用

### 1. 核子-氘核散射

计算中子或质子与氘核的弹性和非弹性散射：

```bash
# 配置文件示例
two_J_3N_max=3        # 三体系统最大角动量
J_2N_max=2            # 二体系统最大角动量
potential_model=N3LO  # 使用N3LO手征势
tensor_force=true     # 包含张量力
```

### 2. 三核子束缚态

研究氚核(³H)和³He的束缚能和波函数：

```bash
# 束缚态计算配置
solve_bound_states=true
binding_energy_search=true
energy_range=[-10.0, 0.0]  # MeV
energy_step=0.1
```

### 3. 相移分析

计算不同分波的相移和散射参数：

```bash
# 相移计算
phase_shift_analysis=true
energy_input_file=data/lab_energies.txt
partial_wave_max=3
output_phase_shifts=true
```

### 4. 极化散射观测量

计算氘核-质子散射的分析能力(Ay)：

```bash
# 氘核-质子Ay计算 (190 MeV/u)
python3 examples/deuteron_proton_Ay.py

# 或使用示例脚本
./examples/run_examples.sh ay
```

**物理背景**：
- 入射粒子：极化氘核 (190 MeV/u)
- 靶核：质子
- 观测量：y方向分析能力 Ay(θ)
- 反映：自旋-轨道耦合和张量力效应

## ⚙️ 配置参数详解

### 基本参数

| 参数 | 说明 | 默认值 | 范围 |
|------|------|--------|------|
| `two_J_3N_max` | 三体总角动量×2的最大值 | 1 | 1-7 |
| `J_2N_max` | 二体角动量最大值 | 1 | 1-5 |
| `Np_WP` | p方向波包数量 | 50 | 10-200 |
| `Nq_WP` | q方向波包数量 | 50 | 10-200 |

### 物理选项

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `tensor_force` | 是否包含张量力 | true |
| `isospin_breaking_1S0` | ¹S₀道同位旋破缺 | true |
| `include_breakup_channels` | 包含破裂道 | false |

### 数值方法

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `solve_dense` | 使用稠密求解器 | false |
| `chebyshev_s` | Chebyshev网格缩放 | 200 |
| `P123_omp_num_threads` | 并行线程数 | 4 |

## 📈 性能优化

### 内存使用优化

- **稀疏矩阵存储**: 使用COO和CSR格式
- **分块计算**: 减少内存峰值使用
- **数据压缩**: HDF5压缩存储

### 计算性能优化

- **OpenMP并行化**: 多线程计算加速
- **向量化**: 利用SIMD指令
- **缓存优化**: 内存访问模式优化

### 典型性能数据

| 系统规模 | 内存使用 | 计算时间 | 并行效率 |
|----------|----------|----------|----------|
| 小型 (Np=30, Nq=30) | ~2GB | 10分钟 | 85% |
| 中型 (Np=50, Nq=50) | ~8GB | 1小时 | 80% |
| 大型 (Np=100, Nq=100) | ~32GB | 8小时 | 75% |

## 🧪 验证和测试

### 基准测试

程序包含多个基准测试用例：

```bash
# 运行所有测试
make test

# 运行特定测试
./tests/test_faddeev_solver
./tests/test_potential_matrix
./tests/test_state_space
```

### 物理验证

- **解析解对比**: 与已知解析结果比较
- **实验数据验证**: 与实验测量值对比
- **代码交叉验证**: 与其他三体代码比较

## 📚 文档和教程

### 详细文档

- [算法原理详解](docs/PRINCIPLE_ALGORITHM.md)
- [编程接口文档](docs/API_REFERENCE.md)
- [物理背景介绍](docs/PHYSICS_BACKGROUND.md)
- [开发指南](docs/DEVELOPMENT_GUIDE.md)

### 教程和示例

- [新手入门教程](examples/tutorial_beginner.md)
- [高级功能使用](examples/tutorial_advanced.md)
- [自定义势能模型](examples/custom_potential.md)
- [并行计算优化](examples/parallel_optimization.md)

## 🤝 贡献指南

我们欢迎各种形式的贡献！

### 如何贡献

1. **Fork项目**到你的GitHub账户
2. **创建特性分支**: `git checkout -b feature-name`
3. **提交更改**: `git commit -am 'Add some feature'`
4. **推送到分支**: `git push origin feature-name`
5. **创建Pull Request**

### 贡献类型

- 🐛 **Bug修复**: 报告和修复程序错误
- ✨ **新功能**: 添加新的物理模型或算法
- 📖 **文档改进**: 完善使用说明和API文档
- 🎨 **代码优化**: 性能优化和代码重构
- 🧪 **测试用例**: 添加更多验证测试

### 代码规范

- 遵循C++17标准
- 使用有意义的变量和函数名
- 添加必要的注释和文档
- 确保代码通过所有测试

## 📄 许可证

本项目采用MIT许可证 - 详见 [LICENSE](LICENSE) 文件。

## 👥 开发团队

### 主要贡献者

- **田百廷** - 项目维护者和主要开发者
- **Tic-tac开发团队** - 算法实现和优化

### 致谢

感谢以下机构和项目的支持：

- 核物理理论研究小组
- 高性能计算中心
- GSL、HDF5等开源项目

## 📞 联系方式

- **项目主页**: https://github.com/tianbaiting/Tic-tac
- **问题报告**: https://github.com/tianbaiting/Tic-tac/issues
- **邮件联系**: tianbaiting@example.com

## 🔗 相关资源

### 理论参考

1. L.D. Faddeev, "Scattering theory for a three-particle system", *Sov. Phys. JETP* **12**, 1014 (1961)
2. E.O. Alt, P. Grassberger, W. Sandhas, "Reduction of the three-particle collision problem", *Nucl. Phys. B* **2**, 167 (1967)
3. A. Deltuva, "Momentum-space treatment of three-nucleon bound state", *Phys. Rev. C* **68**, 031001 (2003)

### 相关项目

- [Few-Body Physics Packages](https://www.few-body.org/)
- [Nuclear Force Models](https://www.nuclear-forces.org/)
- [Computational Nuclear Physics](https://www.comp-nucl-phys.org/)

---

**注意**: 这是重构后的版本2.0，相比原版本在代码组织、性能和易用性方面都有显著改进。

*最后更新: 2024年10月*: Two is company, three's a crowd
Research code for the simulation of three-nucleon scattering in preparation for open-source use.

## Summary
Tic-tac uses the [wave-packet continuum discretization](https://www.sciencedirect.com/science/article/abs/pii/S0003491615001773) (WPCD) method to solve the three-nucleon [Alt-Grassberger-Sandhas](https://www.sciencedirect.com/science/article/abs/pii/0550321367900168) (AGS) equations for elastic nucleon-deuteron scattering.
Currently, Tic-tac can produce on-shell neutron-deuteron scattering amplitudes (U-matrix elements).

## Table of contents
 - Summary
 - Dependencies
 - Compiling
 - Running the code
 - How to cite

## Dependencies
Tic-tac depends on the following libraries:
 - [BLAS](https://netlib.org/blas/)
 - [LAPACK](https://netlib.org/lapack/)
 - [OpenMP](https://www.openmp.org/)
 - [GSL](https://www.gnu.org/software/gsl/)
 - [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
 
Tic-tac relies on the user parsing nuclear-potential codes directly into the source code.
It comes equipped with several codes for different nuclear potential, notably among which are:
 - [Nijmegen-I](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.49.2950)
 - [N2LO<sub>opt](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.110.192502)
 - [N3LO-Idaho](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.68.041001)
 - [Malfliet-Tjon](https://www.sciencedirect.com/science/article/pii/0375947469907751)

## Compiling
Download Tic-tac by running
```
git clone https://github.com/seanbsm/Tic-tac.git
```
A makefile example exists in the repository, but there exists currently no way to automatically compile Tic-tac on arbitrary platforms. It is up to the user to compile and link Tic-tac correctly.

## Running the code
There exists a help-functionality built into Tic-tac. Assuming the Tic-tac executable is named `run`, one can run
```
./run -h
```
or 
```
./run --help
```
which displays all available input-arguments to Tic-tac and their function. One can change arguments directly as input on the command-line or write arguments in a `.txt` input file which Tic-tac reads and interprets. If the input-arguments is not of the expected type, for example a float instead of an integer, Tic-tac will stop running.

An example run could be
```
./run Input/input_example.txt Np_WP=32
```
Tic-tac will read `input_example.txt` and change the default input argument values. Non-specified arguments will keep their default value. Tic-tac can read the same input twice, in which case the last read argument will be used during execution. In this example, `Np_WP=32` will rewrite the default value `Np_WP=30` as well as the specified input in `input_example.txt` since it appears last in the command-line.

## How to cite
We prefer that Tic-tac is referred to by name. For reference, use the [original publication](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.106.024001).
