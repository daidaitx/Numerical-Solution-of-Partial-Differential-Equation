#include "SparseMatrixCSR.hpp"              // CSR稀疏矩阵类
#include "VectorOperations.hpp"             // 稠密向量/矩阵操作
#include "application_utils.hpp"            // 公共工具函数
#include <iostream>
#include <vector>
#include <utility>                          // 用于 pair
#include <cstddef>                          // 用于 size_t
#include <fstream>                          // 文件读写
#include <cmath>                            // 数学函数
#include <cstdlib>                          // 系统函数（用于调用 python 绘图脚本）
#include <functional>                       // 函数对象
#include <algorithm>                        // 排序算法（用于find函数）
#include <stdexcept>                        // 异常处理
#include <optional>                         // 可选值类型
#include <filesystem>                       // 文件系统（用于创建目录）

using namespace std;
using namespace VectorOps;

/**
 * @brief 构造线性方程组的CSR稀疏矩阵 A 和右端项 b
 * @param polygon_vertices 多边形的顶点坐标
 * @param hx 网格在 x 方向的步长
 * @param hy 网格在 y 方向的步长
 * @param grid_points 网格点的坐标列表
 * @param grid_indices 网格点的整数坐标索引列表
 * @param f 源项函数
 * @param g Dirichlet边界条件函数
 * @return 线性方程组的矩阵 A 和右端项 b
 */
pair<SparseMatrixCSR, vector<double>> ConstructLinearSystem(
	const vector<pair<double, double>>& polygon_vertices, 
	double hx, 
	double hy, 
	const vector<pair<double, double>>& grid_points, 
	const vector<pair<size_t, size_t>>& grid_indices, 
	const function<double(double, double)>& f, 
	const function<double(double, double)>& g)
	{
		// 获取网格点的数量
		size_t total_points = grid_points.size();

		// 初始化右端项 b
		vector<double> b(total_points, 0.0);
		
		// 初始化CSR稀疏矩阵 A 的三要素
		vector<double> values(total_points * 5);      // 预分配空间，最多每行5个非零元素
		vector<size_t> col_indices(total_points * 5); // 预分配空间，最多每行5个非零元素
		vector<size_t> row_ptrs(total_points + 1, 0); // row_ptrs 的大小固定为 total_points + 1
	
		// 按x方向构造CSR稀疏矩阵 A 和右端项 b
		for (size_t row_idx = 0; row_idx < total_points; ++row_idx) { // row_idx 既是矩阵行索引，也是当前中心网格点的拉直索引
			/**********************************************************************************************
			 * 检查每一行，记录：                                                                          *
			 * nnz_in_one_row           每行非零元素数   {1, 2, 3, 4, 5}                                  *
			 * idx_*_* < total_points   是否存在四个内部邻近点 {True, False}                              *
			 * idx_*_*                  四个内部邻近点拉直索引 {0, 1, ... , total_points-1, total_points} *
			 * l_*_*                    到四个邻居的距离 [0, h]                                          *
			 * g_*_*                    边界条件        (-∞, +∞)                                        *
			 * f_center                 中心点的源项值   (-∞, +∞)                                       *
			 ******************************************************************************************/

			// 重置 nnz_in_one_row ： 既是矩阵A的当前行的非零元素数，也是当前中心网格点的邻近网格点数，包括中心网格点自己，该值不超过5
			size_t nnz_in_one_row = 0;
			
			// 检查-x邻居
			size_t idx_minus_x = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first - 1, grid_indices[row_idx].second)) - grid_indices.begin(); // idx_minus_x 是 -x 邻居的拉直索引；如果不存在，则等于 total_points
			bool has_minus_x = (idx_minus_x < total_points); // 是否存在-x内部邻近点
			double l_minus_x = hx; // 到-x邻居的距离，默认为 hx，即规则网格距离
			double g_minus_x = 0.0; // -x邻居的边界值，默认为 0.0，即不贡献边界值
			if (has_minus_x) { // 存在-x内部邻近点
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在-x内部邻近点，左侧靠近边界
				// 计算到左侧边界邻交点的距离l_minus_x、左侧边界邻交点的值g_minus_x
				tuple<double, double, pair<double, double>> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "left", hx, g);
				l_minus_x = get<0>(boundary_info);
				g_minus_x = get<1>(boundary_info);
			}

			// 检查-y邻居
			size_t idx_minus_y = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first, grid_indices[row_idx].second - 1)) - grid_indices.begin(); // idx_minus_y 是 -y 邻居的拉直索引；如果不存在，则等于 total_points
			bool has_minus_y = (idx_minus_y < total_points); // 是否存在-y内部邻近点
			double l_minus_y = hy; // 到-y邻居的距离，默认为 hy，即规则网格距离
			double g_minus_y = 0.0; // -y邻居的边界值，默认为 0.0，即不贡献边界值
			if (has_minus_y) { // 存在-y内部邻近点
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在-y内部邻近点，底部靠近边界
				// 计算到底部边界邻交点的距离l_minus_y、底部边界邻交点的值g_minus_y
				tuple<double, double, pair<double, double>> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "bottom", hy, g);
				l_minus_y = get<0>(boundary_info);
				g_minus_y = get<1>(boundary_info);
			}

			// 中心点
			nnz_in_one_row++; // 一定有值，本行非零元素数++
			double f_center = f(grid_points[row_idx].first, grid_points[row_idx].second); // 中心点的源项值

			// 检查+y邻居
			size_t idx_plus_y = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first, grid_indices[row_idx].second + 1)) - grid_indices.begin(); // idx_plus_y 是 +y 邻居的拉直索引；如果不存在，则等于 total_points
			bool has_plus_y = (idx_plus_y < total_points); // 是否存在+y内部邻近点
			double l_plus_y = hy; // 到+y邻居的距离，默认为 hy，即规则网格距离
			double g_plus_y = 0.0; // +y邻居的边界值，默认为 0.0，即不贡献边界值
			if (has_plus_y) { // 存在+y内部邻近点
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在+y内部邻近点，顶部靠近边界
				// 计算到顶部边界邻交点的距离l_plus_y、顶部边界邻交点的值g_plus_y
				tuple<double, double, pair<double, double>> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "top", hy, g);
				l_plus_y = get<0>(boundary_info);
				g_plus_y = get<1>(boundary_info);
			}

			// 检查+x邻居
			size_t idx_plus_x = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first + 1, grid_indices[row_idx].second)) - grid_indices.begin(); // idx_plus_x 是 +x 邻居的拉直索引；如果不存在，则等于 total_points
			bool has_plus_x = (idx_plus_x < total_points); // 是否存在+x内部邻近点
			double l_plus_x = hx; // 到+x邻居的距离，默认为 hx，即规则网格距离
			double g_plus_x = 0.0; // +x邻居的边界值，默认为 0.0，即不贡献边界值
			if (has_plus_x) { // 存在+x内部邻近点
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在+x内部邻近点，右侧靠近边界
				// 计算到右侧边界邻交点的距离l_plus_x、右侧边界邻交点的值g_plus_x
				tuple<double, double, pair<double, double>> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "right", hx, g);
				l_plus_x = get<0>(boundary_info);
				g_plus_x = get<1>(boundary_info);
			}
			
			/*****************************************
			 * 更新 row_ptrs, col_indices, values, b *
			 ****************************************/

			// 递归更新 row_ptrs
			row_ptrs[row_idx + 1] = row_ptrs[row_idx] + nnz_in_one_row;
			
			// 更新 col_indices, values, b
			size_t elem_ptr = row_ptrs[row_idx];

			// -x内部邻近点
			if (has_minus_x) { // 存在-x内部邻近点
				col_indices[elem_ptr] = idx_minus_x; // 列索引
				values[elem_ptr] = - 2.0 / (l_minus_x * (l_minus_x + l_plus_x)); // 元素值
				elem_ptr++;
			}
			else { // 不存在-x内部邻近点，左侧靠近边界
				b[row_idx] += (2.0 * g_minus_x) / (l_minus_x * (l_minus_x + l_plus_x)); // 右端项贡献
			}

			// -y内部邻近点
			if (has_minus_y) { // 存在-y内部邻近点
				col_indices[elem_ptr] = idx_minus_y; // 列索引
				values[elem_ptr] = - 2.0 / (l_minus_y * (l_minus_y + l_plus_y)); // 元素值
				elem_ptr++;
			}
			else { // 不存在-y内部邻近点，底部靠近边界
				b[row_idx] += (2.0 * g_minus_y) / (l_minus_y * (l_minus_y + l_plus_y)); // 右端项贡献
			}

			// 中心元素
			col_indices[elem_ptr] = row_idx; // 列索引 = 行索引，即位于矩阵对角线位置
			values[elem_ptr] = 2.0 / (l_minus_x * l_plus_x) + 2.0 / (l_minus_y * l_plus_y); // 元素值
			elem_ptr++;
			b[row_idx] += f_center; // 右端项贡献

			// +y内部邻近点
			if (has_plus_y) { // 存在+y内部邻近点
				col_indices[elem_ptr] = idx_plus_y; // 列索引
				values[elem_ptr] = - 2.0 / (l_plus_y * (l_minus_y + l_plus_y)); // 元素值
				elem_ptr++;
			}
			else { // 不存在+y内部邻近点，顶部靠近边界
				b[row_idx] += (2.0 * g_plus_y) / (l_plus_y * (l_minus_y + l_plus_y)); // 右端项贡献
			}

			// +x内部邻近点
			if (has_plus_x) { // 存在+x内部邻近点
				col_indices[elem_ptr] = idx_plus_x; // 列索引
				values[elem_ptr] = - 2.0 / (l_plus_x * (l_minus_x + l_plus_x)); // 元素值
				elem_ptr++;
			}
			else { // 不存在+x内部邻近点，右侧靠近边界
				b[row_idx] += (2.0 * g_plus_x) / (l_plus_x * (l_minus_x + l_plus_x)); // 右端项贡献
			}
		}
		
		// 调整values和col_indices的大小
		size_t total_nnz = row_ptrs[total_points];
		values.resize(total_nnz);
		col_indices.resize(total_nnz);

		// 构造和输出
		SparseMatrixCSR A(total_points, total_points, values, col_indices, row_ptrs);
		pair<SparseMatrixCSR, vector<double>> result = {A, b};
		return result;
	}

/**
 * @brief 2维多边形区域Dirichlet边界Poisson问题求解程序
 */
int main() {
	cout << "当前运行的是 2维多边形区域Dirichlet边界Poisson问题 求解程序。" << endl;
	cout << "==============================================================" << endl;
	const string filepath = ".\\applications\\polygon_dirichlet_poisson_solver\\";
	const string datapath = filepath + "data\\";
	if (CheckFolder(filepath, datapath)) return -1;
	cout << "文件路径：" << filepath << endl;
	cout << "数据路径：" << datapath << endl;
	
	cout << "1. 定义多边形区域" << endl;
	cout << "1.1. 逆时针顺序的顶点坐标：" << endl;
	vector<pair<double, double>> polygon_vertices = {
		{0.0, 1.0}, {0.0, -2.0}, {2.0, -2.0}, {1.0, -1.0}, {2.0, 1.0}, {1.0, 2.0}
	};
	cout << "       ";
	for (auto& v : polygon_vertices) {
		cout << "(" << v.first << ", " << v.second << ") ";
	}
	cout << endl;
	cout << "1.2. 将顶点坐标存入文件：" << endl;
	ofstream file(datapath + "polygon_vertices.txt");
	for (auto& v : polygon_vertices) {
		file << v.first << " " << v.second << endl;
	}
	file.close();
	cout << "       已将顶点坐标存入文件：" << datapath + "polygon_vertices.txt" << endl;

	cout << "2. 定义网格" << endl;
	cout << "2.1. 定义网格步长：" << endl;
	const double hx = 0.1;
	const double hy = 0.1;
	cout << "       hx = " << hx << ", hy = " << hy << endl;
	cout << "2.2. 生成网格：" << endl;
	pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> grid_points_and_indices = GenerateInteriorGrid(hx, hy, polygon_vertices);
	vector<pair<double, double>> grid_points = grid_points_and_indices.first;
	vector<pair<size_t, size_t>> grid_indices = grid_points_and_indices.second;
	size_t grid_points_size = grid_points.size();
	cout << "       已生成网格点数量：" << grid_points_size << endl;
	cout << "2.3. 将网格点坐标存入文件：" << endl;
	ofstream fileX(datapath + "grid_X.txt");
	for (auto& v : grid_points) {
		fileX << v.first << " ";
	}
	fileX.close();
	ofstream fileY(datapath + "grid_Y.txt");
	for (auto& v : grid_points) {
		fileY << v.second << " ";
	}
	fileY.close();
	cout << "       已将网格点坐标存入文件：" << datapath + "grid_X.txt" << " 和 " << datapath + "grid_Y.txt" << endl;
	
	cout << "3. 定义解函数、源项函数和边界条件" << endl;
	cout << "3.1. 定义真解：" << endl;
	cout << "3.1.1. 定义真解函数 u(x,y)" << endl;
	auto u_exact = [](double x, double y) { return sin(2.0 * M_PI * x) * sin(M_PI * y); };
	cout << "       已定义真解函数 u(x,y) = sin(2*pi*x)*sin(pi*y)" << endl;
	cout << "3.1.2. 离散真解函数 U(x_i,y_i)：" << endl;
	vector<double> U_exact = DiscretizeExactSolution(u_exact, grid_points);
	cout << "       已离散化真解函数 U(x_i,y_i)" << endl;
	cout << "3.1.3. 将离散真解 U(x_i,y_i) 存入文件：" << endl;
	ofstream fileU_exact(datapath + "U_exact.txt");
	for (double u : U_exact) {
		fileU_exact << u << " ";
	}
	fileU_exact.close();
	cout << "       已将离散真解 U(x_i,y_i) 存入文件：" << datapath + "U_exact.txt" << endl;
	cout << "3.2. 定义源项函数 f(x,y)：" << endl;
	auto f = [](double x, double y) { return 5.0 * M_PI * M_PI * sin(2.0 * M_PI * x) * sin(M_PI * y); };
	cout << "       已定义源项函数 f(x,y) = 5*pi^2*sin(2*pi*x)*sin(pi*y)" << endl;
	cout << "3.3. 定义Dirichlet边界条件 g(x,y)：" << endl;
	auto g = [](double x, double y) { return sin(2.0 * M_PI * x) * sin(M_PI * y); };
	cout << "	    已定义Dirichlet边界条件 g(x,y) = u(x,y)" << endl;

	cout << "4. 构造矩阵 A 和右端项 b" << endl;
	pair<SparseMatrixCSR, vector<double>> linear_system = 
		ConstructLinearSystem(polygon_vertices, hx, hy, grid_points, grid_indices, f, g);
	cout << "4.1. 构造CSR稀疏矩阵：" << endl;
	cout << "4.1.1. 构造矩阵 A：" << endl;
	SparseMatrixCSR A = linear_system.first;
	cout << "         已构造矩阵 A，大小为 " << A.getRows() << " x " << A.getCols() << "，非零元数：" << A.getNNZ() << endl;
	cout << "4.1.2. 将矩阵 A 存入文件：" << endl;
	A.saveToFile(datapath + "matrix_A.coo");
	cout << "         已将矩阵 A 存入文件：" << datapath + "matrix_A.coo" << endl;
	cout << "4.2. 构造右端项 b：" << endl;
	vector<double> b = linear_system.second;
	cout << "       已构造右端项 b，大小为 " << b.size() << endl;

	cout << "5. 求解线性方程组 A U = b" << endl;
	cout << "5.1. 计算线性方程组的解 U_numeric：" << endl;
	vector<double> U_numeric = A.solve(b);
	double epsilon = norm(A * U_numeric - b) / norm(b);
	cout << "       已计算线性方程组的解 U_numeric，大小为 " << U_numeric.size() << "，求解线性方程组的精度为" << epsilon << endl;
	cout << "5.2. 计算相对误差：" << endl;
	double error = norm(U_numeric - U_exact) / norm(U_exact);
	cout << "       在离散步长 hx = " << hx << " , hy = " << hy << " 下，相对误差为：" << error << endl;
	cout << "5.3. 将线性方程组的解 U_numeric 存入文件：" << endl;
	ofstream fileU_numeric(datapath + "U_numeric.txt");
	for (double u : U_numeric) {
		fileU_numeric << u << " ";
	}
	fileU_numeric.close();
	cout << "       已将线性方程组的解 U_numeric 存入文件：" << datapath + "U_numeric.txt" << endl;
	
	cout << "==============================================================" << endl;
	cout << "2维多边形区域Dirichlet边界Poisson问题求解程序 运行结束。" << endl;

	return 0;

	/**
	
	绘制多边形区域及网格： python .\applications\utils\plot_grid.py      polygon_dirichlet_poisson_solver
	绘制真解：            python .\applications\utils\plot_solution.py  polygon_dirichlet_poisson_solver grid_X grid_Y U_exact exact_solution_506
	绘制矩阵 A 的稀疏结构：python .\applications\utils\spy_CSR_matrix.py polygon_dirichlet_poisson_solver matrix_A spy_matrix_A_506
	绘制计算解：          python .\applications\utils\plot_solution.py  polygon_dirichlet_poisson_solver grid_X grid_Y U_numeric numerical_solution_506

	**/
}