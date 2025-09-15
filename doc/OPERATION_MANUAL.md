# Tic-tac 代码操作手册

## 1. 简介

`Tic-tac` 是一个用于核物理研究的计算程序，专门用于求解三核子（3N）系统的Faddeev方程。它可以计算在给定能量下，一个核子与一个束缚的核子对（如氘核）发生散射后的各种可观测量，例如散射截面和极化观测量。

本文档旨在指导一名典型的物理系毕业生完成代码的编译、运行和基本配置。

## 2. 环境要求

在编译和运行此代码之前，请确保您的系统已安装以下软件和库：

- **C++ 编译器**: 需要支持 C++17 标准，例如 `g++` (版本 >= 7)。
- **Fortran 编译器**: 例如 `gfortran`。
- **Make 工具**: 用于自动化编译过程。
- **BLAS (Basic Linear Algebra Subprograms)**: 提供基本的向量和矩阵运算。
- **LAPACK (Linear Algebra PACKage)**: 提供更高级的线性代数例程，如矩阵分解和求解线性方程组。
- **LAPACKE**: LAPACK 的 C 语言接口。
- **GSL (GNU Scientific Library)**: 提供广泛的数学函数库。
- **HDF5 (Hierarchical Data Format 5)**: 用于存储和管理大量科学数据。

在基于 Debian/Ubuntu 的系统上，您可以使用以下命令安装大部分依赖项：
```bash
sudo apt-get update
sudo apt-get install build-essential gfortran libblas-dev liblapack-dev liblapacke-dev libgsl-dev libhdf5-serial-dev
```

## 3. 目录结构

代码库的主要结构如下：

```
.
├── CPP/                  # C++ 源代码主目录
│   ├── General_functions/  # 通用数学和矩阵例程
│   ├── Input/              # 输入文件目录
│   │   └── lab_energies.txt  # 定义散射能量
│   ├── Interactions/       # 定义核子-核子相互作用势
│   ├── Output/             # 默认输出目录（需手动创建）
│   ├── makefile            # 编译脚本
│   ├── main.cpp            # 程序主入口
│   ├── run_organizer.cpp   # 组织和驱动计算流程
│   ├── set_run_parameters.cpp # 解析命令行参数
│   └── ...                 # 其他核心计算模块
└── README.md
```

- `CPP/`: 所有核心逻辑和计算都在这个目录中。
- `CPP/makefile`: 控制整个项目的编译流程。
- `CPP/Input/lab_energies.txt`: 一个重要的输入文件，用于指定入射核子的实验室能量（单位：MeV）。
- `CPP/Output/`: 程序运行时会在此目录下生成结果文件。如果该目录不存在，程序可能会报错。

## 4. 编译与运行

### 4.1 编译

1.  **进入源代码目录**:
    ```bash
    cd /home/tbt/workspace/deutron_compute/F_caculate_3b/Tic-tac/CPP/
    ```

2.  **清理旧文件 (可选但推荐)**:
    如果之前有编译过，建议先清理所有旧的目标文件和可执行文件。
    ```bash
    make cleanall
    ```

3.  **编译项目**:
    运行 `make` 命令开始编译。
    ```bash
    make
    ```
    编译器会逐一编译所有 C++ 和 Fortran 源文件，并最终链接成一个名为 `run` 的可执行文件。在编译过程中，您可能会看到一些警告信息，但只要没有 `Error` 标志，通常可以忽略。

### 4.2 运行

1.  **创建输出目录**:
    程序需要一个名为 `Output` 的目录来存放结果。如果它不存在，请手动创建。
    ```bash
    mkdir -p Output
    ```

2.  **准备输入文件**:
    确保 `Input/lab_energies.txt` 文件存在且包含有效的能量值。例如，要计算 140 MeV 的散射，文件内容应为：
    ```
    140.0
    ```

3.  **执行程序**:
    在 `CPP` 目录下运行可执行文件：
    ```bash
    ./run
    ```
    程序将开始运行，并在终端上打印详细的计算进度，包括基空间构建、置换矩阵计算、求解Faddeev方程等步骤。

## 5. 输入与输出

### 5.1 输入

- **`Input/lab_energies.txt`**: 定义一个或多个入射核子的实验室能量（MeV）。每行一个能量值。程序会为每个能量值进行一次完整的Faddeev方程求解。
- **命令行参数**: 程序支持通过命令行参数进行详细配置。运行 `./run -h` 或 `./run --help` 可以查看所有可用选项。

### 5.2 输出

- **标准输出 (Terminal)**: 程序运行时会将详细的日志信息打印到终端，包括：
    - 当前运行的参数配置。
    - 3N 分波基的构建情况。
    - 波包（WP）基的构建情况。
    - 置换矩阵 `P123` 的计算进度和内存占用。
    - Faddeev 方程的迭代求解过程。
- **`Output/` 目录**:
    - `run_parameters.txt`: 记录本次运行所使用的所有参数。
    - `P123_sparse_*.h5`: 以 HDF5 格式存储的稀疏置换矩阵。这是一个关键的中间结果，计算非常耗时，因此会被存储起来以便重用。
    - 其他 `.csv` 或 `.h5` 文件：存储了波包边界、SWP（散射波包）能量等信息。

## 6. 高级配置

`Tic-tac` 提供了丰富的命令行参数来微调计算的各个方面。详细的参数列表和解释请参考 `ADVANCED_CONFIGURATION.md` 文档。

## 7. 常见问题解答

- **编译错误: `fatal error: some_header.h: No such file or directory`**
  > 这通常意味着缺少某个依赖库的头文件。请根据错误信息中提到的文件名，确认您已正确安装了对应的开发库（通常以 `-dev` 或 `-devel` 结尾）。

- **链接错误: `undefined reference to 'cblas_dgemm'` 或 `undefined reference to 'LAPACKE_dspevd'`**
  > 这意味着链接器找不到 BLAS 或 LAPACK/LAPACKE 库。请检查 `makefile` 中的 `LDLIBS` 变量，确保 `-lblas`、`-llapack` 和 `-llapacke` 标志存在且顺序正确。

- **运行错误: `Unable to open file Output/run_parameters.txt`**
  > 程序无法写入输出文件。这通常是因为 `Output` 目录不存在。请在 `CPP` 目录下手动创建该目录：`mkdir Output`。

- **运行错误: `Unable to open file lab_energies.txt`**
  > 程序找不到能量输入文件。请确保 `Input/lab_energies.txt` 文件存在于正确的位置。

- **计算速度慢/内存占用高**
  > 三体计算本质上是计算密集型和内存密集型的。特别是 `P123` 置换矩阵的计算，其大小与基空间维度（`Np_WP`, `Nq_WP`, `Nphi`, `Nx`）高度相关。如果遇到性能问题，可以尝试：
  > 1.  在 `makefile` 中调整 `CPPFLAGS` 以启用更高等级的优化（如 `-O3`）。
  > 2.  通过命令行参数减小基空间维度（例如，减小 `Np_WP`, `Nq_WP` 或 `J_2N_max`）。
  > 3.  在计算 `P123` 时，使用 `P123_omp_num_threads` 参数开启 OpenMP 并行计算。

