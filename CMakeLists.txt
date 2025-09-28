# CMakeLists.txt 示例核心配置
project(PDE_SparseMatrixLib CXX)
set(CMAKE_CXX_STANDARD 17)

# 核心库
add_library(SparseMatrixLib STATIC
  src/SparseMatrixCSR.cpp
  src/MatrixGenerators.cpp
)
target_include_directories(SparseMatrixLib PUBLIC include)

# 测试套件
add_executable(MatrixTests test/TestCore.cpp)
target_link_libraries(MatrixTests SparseMatrixLib)