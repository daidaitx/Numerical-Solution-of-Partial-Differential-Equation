# CMakeLists.txt
cmake_minimum_required(VERSION 3.10)  # 指定CMake最低版本要求（确保兼容性）
project(NPDE_SparseMatrixLib CXX)     # 定义项目名称和语言（CXX表示C++）
set(CMAKE_CXX_STANDARD 17)            # 强制使用C++17标准编译

# 核心库：创建静态库目标
add_library(SparseMatrixLib STATIC    # 生成名为SparseMatrixLib的静态库
  sparse_matrix_lib/src/SparseMatrixCSR.cpp           # 库源文件1：稀疏矩阵CSR格式实现
  sparse_matrix_lib/src/VectorOperations.cpp          # 库源文件2：向量运算实现
  sparse_matrix_lib/src/IterativeSolvers.cpp          # 库源文件3：迭代求解器实现
  sparse_matrix_lib/src/SpecificMatrixGenerators.cpp  # 库源文件4：特殊矩阵生成器实现
)

# 设置库的头文件搜索路径
target_include_directories(SparseMatrixLib PUBLIC sparse_matrix_lib/include) 
# PUBLIC表示：
#   - 库自身编译时需要include目录
#   - 链接此库的**其他目标**（如测试程序）也需要include目录

# 自动查找所有测试文件
file(GLOB TEST_SOURCES "sparse_matrix_lib/test/Test*.cpp")
# 自动查找所有应用程序
file(GLOB APP_SOURCES "applications/*/*.cpp")

# 自动创建测试可执行文件
foreach(test_source ${TEST_SOURCES})
    get_filename_component(test_name ${test_source} NAME_WE)
    add_executable(${test_name} ${test_source})
    target_link_libraries(${test_name} PRIVATE SparseMatrixLib)
    list(APPEND ALL_EXECUTABLES ${test_name})
endforeach()

# 自动创建应用程序可执行文件
foreach(app_source ${APP_SOURCES})
    get_filename_component(app_name ${app_source} NAME_WE)
    add_executable(${app_name} ${app_source})
    target_link_libraries(${app_name} PRIVATE SparseMatrixLib)
    list(APPEND ALL_EXECUTABLES ${app_name})
endforeach()