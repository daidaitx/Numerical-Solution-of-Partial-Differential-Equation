#pragma once
#include "SparseMatrixCSR.hpp"
#include <vector>

// 前向声明，避免包含整个头文件
class SparseMatrixCSR;

/**
 * @brief 稀疏矩阵求解器命名空间
 */
namespace SparseMatrixSolvers {

/**
 * @brief 使用GMRES法求解线性方程组 Ax = b
 * @param A 系数矩阵 (CSR格式)
 * @param b 右端项向量
 * @param x0 初始解（可选）
 * @param restart 重启步数
 * @param max_restarts 最大重启次数
 * @param tol 收敛容差
 * @return 解向量
 */
std::vector<double> solveGMRES(const SparseMatrixCSR& A,
							const std::vector<double>& b,
							const std::vector<double>& x0 = {},
							size_t restart = 30,
							size_t max_restarts = 100,
							double tol = 1e-8);

} // namespace SparseMatrixSolvers