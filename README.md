# 数值求解偏微分方程 (Numerical Solution of Partial Differential Equation)

本项目旨在实现二维区域上偏微分方程的数值求解，特别是泊松(Poisson)方程。通过构建高效的稀疏矩阵数据结构和求解器，支持大规模线性系统的快速计算。

## 项目概述

本项目提供了一个完整的数值求解偏微分方程的框架，包含以下核心组件：

1. **稀疏矩阵库** - 基于CSR(Compressed Sparse Row)格式实现的高效稀疏矩阵类
2. **线性求解器** - 实现GMRES等迭代法求解大型稀疏线性系统
3. **应用示例** - 多边形区域上泊松方程的五点差分求解器，支持Dirichlet和Neumann边界条件
4. **可视化工具** - Python脚本用于结果可视化和分析

## 核心功能

### 稀疏矩阵库 (Sparse Matrix Library)
位于 `sparse_matrix_lib` 目录下，包含以下主要功能：

- CSR格式稀疏矩阵的创建、读写、运算与属性判断
- 稠密/稀疏矩阵转换、向量运算支持
- 支持COO文件格式读写
- 矩阵运算：加减、数乘、转置、范数计算、稀疏矩阵-向量乘法(SpMV)
- 线性求解器：GMRES迭代法求解大型稀疏线性系统

### 应用程序 (Applications)
位于 `applications` 目录下，包含以下应用示例：

#### 1. 多边形区域Dirichlet边界条件泊松方程求解器
目录：`applications/polygon_dirichlet_poisson_solver`

该程序实现了在二维多边形区域内求解泊松方程，具有以下特点：
- 支持任意形状的多边形区域（顶点按逆时针顺序给出）
- 使用五点差分格式离散化
- 处理Dirichlet边界条件
- 自动处理不规则网格点（靠近边界点）
- 计算精度高，支持相对误差分析

#### 2. 多边形区域Neumann边界条件泊松方程求解器
目录：`applications/polygon_neumann_poisson_solver`

该程序实现了在二维多边形区域内求解泊松方程，具有以下特点：
- 支持任意形状的多边形区域（顶点按逆时针顺序给出）
- 使用五点差分格式离散化
- 处理Neumann边界条件
- 自动处理不规则网格点（靠近边界点）
- 采用特殊处理确保Neumann问题解的唯一性

#### 3. 工具函数库
目录：`applications/utils`

包含通用工具函数和可视化脚本：
- 网格生成工具
- 文件操作工具
- Python可视化脚本（绘图、矩阵稀疏结构显示等）

### 可视化工具
提供Python脚本用于结果可视化：
- `plot_grid.py` - 绘制网格点和区域边界
- `plot_solution.py` - 绘制数值解和真解
- `spy_CSR_matrix.py` - 显示稀疏矩阵的非零元素结构

## 编译和运行

### 环境要求
- C++17兼容编译器
- CMake 3.10或更高版本
- Python 3（用于可视化）
- matplotlib（Python绘图库）

### 编译步骤
```bash
# 创建构建目录
mkdir build && cd build

# 生成Makefile
cmake ..

# 编译所有应用程序
cmake --build .
# 或使用 make
```

### 运行示例
编译完成后，可执行文件将生成在build目录下：

```bash
# 运行Dirichlet边界条件泊松方程求解器
./applications/polygon_dirichlet_poisson_solver/2D_polygonal_Dirichlet_Poisson_problem

# 运行Neumann边界条件泊松方程求解器
./applications/polygon_neumann_poisson_solver/2D_polygonal_Neumann_Poisson_problem
```

### 可视化结果
运行求解器后，可使用Python脚本可视化结果：

```bash
# 绘制网格
python ./applications/utils/plot_grid.py polygon_dirichlet_poisson_solver

# 绘制解
python ./applications/utils/plot_solution.py polygon_dirichlet_poisson_solver grid_X grid_Y U_exact exact_solution

# 显示矩阵稀疏结构
python ./applications/utils/spy_CSR_matrix.py polygon_dirichlet_poisson_solver matrix_A spy_matrix
```

## 项目结构

```
.
├── sparse_matrix_lib/           # 稀疏矩阵库
│   ├── include/                 # 头文件
│   ├── src/                     # 源文件
│   └── test/                    # 测试文件
├── applications/                # 应用程序
│   ├── polygon_dirichlet_poisson_solver/  # Dirichlet边界条件求解器
│   ├── polygon_neumann_poisson_solver/    # Neumann边界条件求解器
│   └── utils/                   # 工具函数和可视化脚本
├── docs/                       # 文档
├── build/                      # 构建目录
├── CMakeLists.txt              # CMake构建配置
└── README.md                   # 项目说明文件
```

## 文档

项目使用Doxygen生成代码文档。可在`docs/html/index.html`查看完整文档。

## 测试

项目包含多个测试程序，用于验证稀疏矩阵库的正确性和性能：
- TestCore - 测试稀疏矩阵基本功能
- TestSpMV - 测试稀疏矩阵向量乘法
- TestPoissonMatrix - 测试泊松矩阵生成

## 许可证

本项目仅供学习和研究使用。
