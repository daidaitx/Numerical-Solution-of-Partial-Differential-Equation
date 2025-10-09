#include "SparseMatrixCSR.hpp"
#include "VectorOperations.hpp"
#include "SpecificMatrixGenerators.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;
using namespace MatrixGenerators::Poisson2D;

int main() {
	cout << "当前运行的是 TestPoissonMatrix 测试程序。" << endl;
	cout << "=========================================" << endl;
	cout << "测试1 - 以函数输入生成简单二维 Poisson 方程的稀疏矩阵和右端项。" << endl;
	// 定义网格尺寸
	const size_t M_1 = 5;
	const size_t N_1 = 5;
	const double x_min_1 = 0.0;
	const double x_max_1 = 1.0;
	const double y_min_1 = 0.0;
	const double y_max_1 = 1.0;
	// 定义 Poisson 方程的配置
	const PoissonConfig2D config_1 = { M_1, N_1, x_min_1, x_max_1, y_min_1, y_max_1 };
	// 预设解 u(x,y) = sin(2*pi*x)*sin(2*pi*y)
	auto u_exact_1 = [](double x, double y) { return sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y); };
	// 定义源项函数 f(x,y) 和边界条件函数 g(x,y)
	// f(x,y) =  -Δu   = 8*pi^2*sin(2*pi*x)*sin(2*pi*y)
	// g(x,y) = u(x,y) = sin(2*pi*x)*sin(2*pi*y)
	auto f_1 = [](double x, double y) { return 8.0 * M_PI * M_PI * sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y); };
	auto g_1 = [](double x, double y) { return sin(2.0 * M_PI * x) * sin(2.0 * M_PI * y); };
	// 生成 Poisson 问题
	PoissonProblem2D problem_1 = generatePoisson2D(config_1, f_1, g_1);
	// 输出矩阵和右端项信息
	problem_1.printInfo();
	// // 保存矩阵和右端项到文件
	// problem_1.saveToFiles("poisson_2d_1");
	// // 读取矩阵和右端项的文件
	// SparseMatrixCSR A_1;

	cout << "=========================================" << endl;
	cout << "TestPoissonMatrix 测试结束。" << endl;
	return 0;
}