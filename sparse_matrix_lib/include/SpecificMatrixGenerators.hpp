#pragma once
#include "SparseMatrixCSR.hpp"
#include <vector>
#include <functional>
#include <string>

/**
 * @brief 特定矩阵生成器命名空间
 * @details 用于生成各种特定问题的稀疏矩阵
 */
namespace MatrixGenerators {

/**
 * @brief 二维Poisson问题生成器
 */
namespace Poisson2D {

// 前向声明主要结构体
struct PoissonProblem2D;

/**
 * @brief 二维Poisson问题配置参数
 */
struct PoissonConfig2D {
	size_t M;                    ///< x方向网格数
	size_t N;                    ///< y方向网格数
	double x_min, x_max;         ///< x方向边界
	double y_min, y_max;         ///< y方向边界
	std::string discretization = "5-point"; ///< 离散格式
	
	// 网格步长
	double hx() const { return (x_max - x_min) / M; }
	double hy() const { return (y_max - y_min) / N; }
};

/**
 * @brief 二维Poisson问题完整结果
 */
struct PoissonProblem2D {
	SparseMatrixCSR A;                   ///< 系数矩阵
	std::vector<double> b;               ///< 右端项
	std::vector<std::vector<double>> solution_2d; ///< 二维格式的解
	PoissonConfig2D config;              ///< 配置参数
	
	// 实用函数
	void saveToFiles(const std::string& basename) const;
	void printInfo() const;
};

// 函数重载：接受函数句柄
PoissonProblem2D generatePoisson2D(
	const PoissonConfig2D& config,
	std::function<double(double, double)> f,  // 源项函数
	std::function<double(double, double)> g   // 边界条件函数
);

// 函数重载：接受离散点阵
PoissonProblem2D generatePoisson2D(
	const PoissonConfig2D& config,
	const std::vector<std::vector<double>>& f_discrete,  // 离散源项
	const std::vector<std::vector<double>>& g_discrete   // 离散边界条件
);

// 仅生成矩阵（不求解）
SparseMatrixCSR generatePoissonMatrix2D(const PoissonConfig2D& config);

// 仅生成右端项
std::vector<double> generatePoissonRHS2D(
	const PoissonConfig2D& config,
	std::function<double(double, double)> f,
	std::function<double(double, double)> g
);

std::vector<double> generatePoissonRHS2D(
	const PoissonConfig2D& config,
	const std::vector<std::vector<double>>& f_discrete,
	const std::vector<std::vector<double>>& g_discrete
);

// 工具函数：向量拉直和反拉直
std::vector<double> flattenMatrix(const std::vector<std::vector<double>>& matrix_2d);
std::vector<std::vector<double>> unflattenVector(
	const std::vector<double>& vec_1d, 
	size_t M, 
	size_t N
);

} // namespace Poisson2D

} // namespace MatrixGenerators