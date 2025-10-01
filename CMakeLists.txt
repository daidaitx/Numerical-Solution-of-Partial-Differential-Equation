# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)  # 指定CMake最低版本要求（确保兼容性）
project(NPDE_SparseMatrixLib CXX)     # 定义项目名称和语言（CXX表示C++）
set(CMAKE_CXX_STANDARD 17)            # 强制使用C++17标准编译

# 核心库：创建静态库目标
add_library(SparseMatrixLib STATIC    # 生成名为SparseMatrixLib的静态库
  sparse_matrix_lib/src/SparseMatrixCSR.cpp           # 库源文件1：稀疏矩阵CSR格式实现
  sparse_matrix_lib/src/SpecificMatrixGenerators.cpp  # 库源文件2：特殊矩阵生成器实现
  sparse_matrix_lib/src/FileIO.cpp                    # 库源文件3：文件输入输出实现
)

# 设置库的头文件搜索路径
target_include_directories(SparseMatrixLib PUBLIC sparse_matrix_lib/include) 
# PUBLIC表示：
#   - 库自身编译时需要include目录
#   - 链接此库的**其他目标**（如测试程序）也需要include目录

# 测试套件：每个测试文件生成独立可执行文件
add_executable(TestCore                          # 生成名为TestCore的可执行文件
  sparse_matrix_lib/test/TestCore.cpp)           # 测试主程序源文件
add_executable(TestPoissonMatrix                 # 生成名为TestPoissonMatrix的可执行文件
  sparse_matrix_lib/test/TestPoissonMatrix.cpp)  # 测试Poisson矩阵生成器源文件
add_executable(TestSpMV                          # 生成名为TestSpMV的可执行文件
  sparse_matrix_lib/test/TestSpMV.cpp)           # 测试稀疏矩阵向量乘法源文件

# 统一链接库
foreach(test_exe IN ITEMS TestCore TestPoissonMatrix TestSpMV)
  target_link_libraries(${test_exe} PRIVATE SparseMatrixLib)
endforeach()
# 使可执行文件能够访问SparseMatrixLib中的函数和类

# 注释说明：
# 1. 静态库(.a/.lib)会将所有源文件编译后打包成一个文件，适合封装通用功能
# 2. PUBLIC关键字确保头文件路径能被传递到依赖目标（测试程序），避免手动添加路径
# 3. 测试程序通过链接库实现模块化测试，符合"编译一次，多处运行"的测试理念
# 4. 此配置支持跨平台（Windows/Linux/macOS）构建，CMake会根据系统选择编译器