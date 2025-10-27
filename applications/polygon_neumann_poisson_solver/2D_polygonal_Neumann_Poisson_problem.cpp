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
#include <numeric>                          // accumulate 函数

using namespace std;
using namespace VectorOps;

/**
 * @brief 计算多边形区域的边界单位外法向量
 * @param polygon_vertices 多边形的顶点坐标
 * @return 多边形区域的边界单位外法向量
 * @attention 多边形顶点坐标必须按逆时针顺序排列
 */
vector<pair<double, double>> ComputePolygonBoundaryNormals(const vector<pair<double, double>>& polygon_vertices) {
	vector<pair<double, double>> normals;
	for (size_t i = 0; i < polygon_vertices.size(); ++i) {
		size_t j = (i + 1) % polygon_vertices.size();
		double dx = polygon_vertices[j].first - polygon_vertices[i].first;
		double dy = polygon_vertices[j].second - polygon_vertices[i].second;
		double normal_x = dy;
		double normal_y = -dx;
		double normal_norm = sqrt(normal_x * normal_x + normal_y * normal_y);
		normals.push_back({normal_x / normal_norm, normal_y / normal_norm});
	}
	return normals;
}

/**
 * @brief 判断点是否在线段上（含端点）
 * @param x 点的横坐标
 * @param y 点的纵坐标
 * @param A 线段端点A的坐标
 * @param B 线段端点B的坐标
 * @return 点是否在线段上
 * @note 本函数使用 1e-10 的容差值处理浮点数精度问题
 */
bool isAtSegment(double x, double y, const pair<double, double>& A, const pair<double, double>& B) {
	double dx = B.first - A.first;
	double dy = B.second - A.second;
	double len2 = dx * dx + dy * dy;
	double t = ((x - A.first) * dx + (y - A.second) * dy) / len2;
	if (t < 0.0 - 1e-10 || t > 1.0 + 1e-10) { // 点不在线段上
		return false;
	}
	double x_proj = A.first + t * dx;
	double y_proj = A.second + t * dy;
	double dist2 = (x - x_proj) * (x - x_proj) + (y - y_proj) * (y - y_proj);
	if (dist2 > 1e-10) { // 点不在线段上
		return false;
	}
	return true;
}

/**
 * @brief 构造线性方程组的CSR稀疏矩阵 A 和右端项 b
 * @param polygon_vertices 多边形的顶点坐标
 * @param hx 网格在 x 方向的步长
 * @param hy 网格在 y 方向的步长
 * @param grid_points 网格点的坐标列表
 * @param grid_indices 网格点的整数坐标索引列表
 * @param f 源项函数
 * @param g Neumann边界条件函数
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
			double g_minus_x = 0.0; // -x邻居的Neumann边界值，默认为 0.0，即不贡献边界值
			double nu_x_minus_x = 0.0; // -x邻居的单位外法向量（横坐标），默认为 0.0，表示-x内部邻近点存在，-x邻居不在边界上
			double nu_y_minus_x = 0.0; // -x邻居的单位外法向量（纵坐标），默认为 0.0，表示-x内部邻近点存在，-x邻居不在边界上
			if (has_minus_x) { // 存在-x内部邻近点
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在-x内部邻近点，左侧靠近边界
				// 计算到左侧边界邻交点的距离l_minus_x、左侧边界邻交点的Neumann边界值g_minus_x、单位外法向量
				tuple<double, double, pair<double, double>> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "left", hx, g);
				l_minus_x = get<0>(boundary_info);
				g_minus_x = get<1>(boundary_info);
				nu_x_minus_x = get<2>(boundary_info).first;
				nu_y_minus_x = get<2>(boundary_info).second;
			}

			// 检查-y邻居
			size_t idx_minus_y = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first, grid_indices[row_idx].second - 1)) - grid_indices.begin(); // idx_minus_y 是 -y 邻居的拉直索引；如果不存在，则等于 total_points
			bool has_minus_y = (idx_minus_y < total_points); // 是否存在-y内部邻近点
			double l_minus_y = hy; // 到-y邻居的距离，默认为 hy，即规则网格距离
			double g_minus_y = 0.0; // -y邻居的边界值，默认为 0.0，即不贡献边界值
			double nu_x_minus_y = 0.0; // -y邻居的单位外法向量（横坐标），默认为 0.0，表示-y内部邻近点存在，-y邻居不在边界上
			double nu_y_minus_y = 0.0; // -y邻居的单位外法向量（纵坐标），默认为 0.0，表示-y内部邻近点存在，-y邻居不在边界上
			if (has_minus_y) { // 存在-y内部邻近点
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在-y邻居，底部靠近边界
				// 计算到底部边界邻交点的距离l_minus_y、底部边界邻交点的Neumann边界值g_minus_y、单位外法向量
				tuple<double, double, pair<double, double>> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "bottom", hy, g);
				l_minus_y = get<0>(boundary_info);
				g_minus_y = get<1>(boundary_info);
				nu_x_minus_y = get<2>(boundary_info).first;
				nu_y_minus_y = get<2>(boundary_info).second;
			}

			// 中心点
			nnz_in_one_row++; // 一定有值，本行非零元素数++
			double f_center = f(grid_points[row_idx].first, grid_points[row_idx].second); // 中心点的源项值

			// 检查+y邻居
			size_t idx_plus_y = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first, grid_indices[row_idx].second + 1)) - grid_indices.begin(); // idx_plus_y 是 +y 邻居的拉直索引；如果不存在，则等于 total_points
			bool has_plus_y = (idx_plus_y < total_points); // 是否存在+y内部邻近点
			double l_plus_y = hy; // 到+y邻居的距离，默认为 hy，即规则网格距离
			double g_plus_y = 0.0; // +y邻居的边界值，默认为 0.0，即不贡献边界值
			double nu_x_plus_y = 0.0; // +y邻居的单位外法向量（横坐标），默认为 0.0，表示+y内部邻近点存在，+y邻居不在边界上
			double nu_y_plus_y = 0.0; // +y邻居的单位外法向量（纵坐标），默认为 0.0，表示+y内部邻近点存在，+y邻居不在边界上
			if (has_plus_y) { // 存在+y内部邻近点
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在+y邻居，顶部靠近边界
				// 计算到顶部边界邻交点的距离l_plus_y、顶部边界邻交点的Neumann边界值g_plus_y、单位外法向量
				tuple<double, double, pair<double, double>> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "top", hy, g);
				l_plus_y = get<0>(boundary_info);
				g_plus_y = get<1>(boundary_info);
				nu_x_plus_y = get<2>(boundary_info).first;
				nu_y_plus_y = get<2>(boundary_info).second;
			}

			// 检查+x邻居
			size_t idx_plus_x = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first + 1, grid_indices[row_idx].second)) - grid_indices.begin(); // idx_plus_x 是 +x 邻居的拉直索引；如果不存在，则等于 total_points
			bool has_plus_x = (idx_plus_x < total_points); // 是否存在+x内部邻近点
			double l_plus_x = hx; // 到+x邻居的距离，默认为 hx，即规则网格距离
			double g_plus_x = 0.0; // +x邻居的边界值，默认为 0.0，即不贡献边界值
			double nu_x_plus_x = 0.0; // +x邻居的单位外法向量（横坐标），默认为 0.0，表示+x内部邻近点存在，+x邻居不在边界上
			double nu_y_plus_x = 0.0; // +x邻居的单位外法向量（纵坐标），默认为 0.0，表示+x内部邻近点存在，+x邻居不在边界上
			if (has_plus_x) { // 存在+x内部邻近点
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在+x邻居，右侧靠近边界
				// 计算到右侧边界邻交点的距离l_plus_x、右侧边界邻交点的Neumann边界值g_plus_x、单位外法向量
				tuple<double, double, pair<double, double>> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "right", hx, g);
				l_plus_x = get<0>(boundary_info);
				g_plus_x = get<1>(boundary_info);
				nu_x_plus_x = get<2>(boundary_info).first;
				nu_y_plus_x = get<2>(boundary_info).second;
			}
			
			/*****************************************
			 * 更新 row_ptrs, col_indices, values, b *
			 ****************************************/

			// 递归更新 row_ptrs
			row_ptrs[row_idx + 1] = row_ptrs[row_idx] + nnz_in_one_row;
			
			// 更新 col_indices, values, b
			size_t elem_ptr = row_ptrs[row_idx];

			size_t elem_ptr_minus_x; // -x邻居的元素指针
			size_t elem_ptr_minus_y; // -y邻居的元素指针
			size_t elem_ptr_center; // 中心元素的元素指针
			size_t elem_ptr_plus_y; // +y邻居的元素指针
			size_t elem_ptr_plus_x; // +x邻居的元素指针

			// 存在-x内部邻近点
			if (has_minus_x) {
				col_indices[elem_ptr] = idx_minus_x; // 列索引
				values[elem_ptr] = - 2.0 / (l_minus_x * (l_minus_x + l_plus_x)); // 元素值
				elem_ptr_minus_x = elem_ptr;
				elem_ptr++;
			}

			// 存在-y内部邻近点
			if (has_minus_y) {
				col_indices[elem_ptr] = idx_minus_y; // 列索引
				values[elem_ptr] = - 2.0 / (l_minus_y * (l_minus_y + l_plus_y)); // 元素值
				elem_ptr_minus_y = elem_ptr;
				elem_ptr++;
			}

			// 中心元素
			col_indices[elem_ptr] = row_idx; // 列索引 = 行索引，即位于矩阵对角线位置
			values[elem_ptr] = 2.0 / (l_minus_x * l_plus_x) + 2.0 / (l_minus_y * l_plus_y); // 元素值
			elem_ptr_center = elem_ptr;
			elem_ptr++;
			b[row_idx] += f_center; // 右端项贡献

			// 存在+y内部邻近点
			if (has_plus_y) {
				col_indices[elem_ptr] = idx_plus_y; // 列索引
				values[elem_ptr] = - 2.0 / (l_plus_y * (l_minus_y + l_plus_y)); // 元素值
				elem_ptr_plus_y = elem_ptr;
				elem_ptr++;
			}

			// 存在+x内部邻近点
			if (has_plus_x) {
				col_indices[elem_ptr] = idx_plus_x; // 列索引
				values[elem_ptr] = - 2.0 / (l_plus_x * (l_minus_x + l_plus_x)); // 元素值
				elem_ptr_plus_x = elem_ptr;
				elem_ptr++;
			}

			/*****************
			 * 处理非正则内点 *
			 ****************/

			double extra_minus_x = 0.0;
			double extra_minus_y = 0.0;
			
			// 存在-x边界邻交点
			if (!has_minus_x) { // -x内部邻近点不存在，即-x邻居在边界上，存在-x边界邻交点
				double nu_l_l = nu_x_minus_x * (l_minus_x + l_plus_x);
				if (has_minus_y) { // -y内部邻近点存在
					b[row_idx] -= 2.0 * g_minus_x / nu_l_l;
					values[elem_ptr_minus_y] += 2.0 * nu_y_minus_x / (l_minus_y * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_x_minus_x * l_minus_y + nu_y_minus_x * l_minus_x) / (l_minus_y * l_minus_x * nu_l_l);
				}
				else if (has_plus_y) { // -y内部邻近点不存在，+y内部邻近点存在
					b[row_idx] -= 2.0 * g_minus_x / nu_l_l;
					values[elem_ptr_plus_y] += -2.0 * nu_y_minus_x / (l_plus_y * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_x_minus_x * l_plus_y - nu_y_minus_x * l_minus_x) / (l_plus_y * l_minus_x * nu_l_l);
				}
				else { // -y和+y内部邻近点都不存在
					// values[elem_ptr_center] -= 2.0 / (l_minus_x * (l_minus_x + l_plus_x));
					b[row_idx] -= 2.0 * g_minus_x / nu_l_l;
					extra_minus_y += 2.0 * nu_y_minus_x / (l_minus_y * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_x_minus_x * l_minus_y + nu_y_minus_x * l_minus_x) / (l_minus_y * l_minus_x * nu_l_l);
				}
			}

			// 存在-y边界邻交点
			if (!has_minus_y) { // -y内部邻近点不存在，即-y邻居在边界上，存在-y边界邻交点
				double nu_l_l = nu_y_minus_y * (l_minus_y + l_plus_y);
				if (has_minus_x) { // -x内部邻近点存在
					b[row_idx] -= 2.0 * g_minus_y / nu_l_l;
					values[elem_ptr_minus_x] += 2.0 * nu_x_minus_y / (l_minus_x * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_y_minus_y * l_minus_x + nu_x_minus_y * l_minus_y) / (l_minus_x * l_minus_y * nu_l_l);
				}
				else if (has_plus_x) { // -x内部邻近点不存在，+x内部邻近点存在
					b[row_idx] -= 2.0 * g_minus_y / nu_l_l;
					values[elem_ptr_plus_x] += -2.0 * nu_x_minus_y / (l_plus_x * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_y_minus_y * l_plus_x - nu_x_minus_y * l_minus_y) / (l_plus_x * l_minus_y * nu_l_l);
				}
				else { // -x和+x内部邻近点都不存在，作零梯度假设
					// values[elem_ptr_center] -= 2.0 / (l_minus_y * (l_minus_y + l_plus_y));
					b[row_idx] -= 2.0 * g_minus_y / nu_l_l;
					extra_minus_x += 2.0 * nu_x_minus_y / (l_minus_x * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_y_minus_y * l_minus_x + nu_x_minus_y * l_minus_y) / (l_minus_x * l_minus_y * nu_l_l);
				}
			}

			// 存在+y边界邻交点
			if (!has_plus_y) { // +y内部邻近点不存在，即+y邻居在边界上，存在+y边界邻交点
				double nu_l_l = nu_y_plus_y * (l_minus_y + l_plus_y);
				if (has_minus_x) { // -x内部邻近点存在
					b[row_idx] -= -2.0 * g_plus_y / nu_l_l;
					values[elem_ptr_minus_x] += -2.0 * nu_x_plus_y / (l_minus_x * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_y_plus_y * l_minus_x - nu_x_plus_y * l_plus_y) / (l_minus_x * l_plus_y * nu_l_l);
				}
				else if (has_plus_x) { // -x内部邻近点不存在，+x内部邻近点存在
					b[row_idx] -= -2.0 * g_plus_y / nu_l_l;
					values[elem_ptr_plus_x] += 2.0 * nu_x_plus_y / (l_plus_x * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_y_plus_y * l_plus_x + nu_x_plus_y * l_plus_y) / (l_plus_x * l_plus_y * nu_l_l);
				}
				else { // -x和+x内部邻近点都不存在，作零梯度假设
					// values[elem_ptr_center] -= 2.0 / (l_plus_y * (l_minus_y + l_plus_y));
					b[row_idx] -= -2.0 * g_plus_y / nu_l_l;
					extra_minus_x += -2.0 * nu_x_plus_y / (l_minus_x * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_y_plus_y * l_minus_x - nu_x_plus_y * l_plus_y) / (l_minus_x * l_plus_y * nu_l_l);
				}
			}

			// 存在+x边界邻交点
			if (!has_plus_x) { // +x内部邻近点不存在，即+x邻居在边界上，存在+x边界邻交点
				double nu_l_l = nu_x_plus_x * (l_minus_x + l_plus_x);
				if (has_minus_y) { // -y内部邻近点存在
					b[row_idx] -= -2.0 * g_plus_x / nu_l_l;
					values[elem_ptr_minus_y] += -2.0 * nu_y_plus_x / (l_minus_y * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_x_plus_x * l_minus_y - nu_y_plus_x * l_plus_x) / (l_minus_y * l_plus_x * nu_l_l);
				}
				else if (has_plus_y) { // -y内部邻近点不存在，+y内部邻近点存在
					b[row_idx] -= -2.0 * g_plus_x / nu_l_l;
					values[elem_ptr_plus_y] += 2.0 * nu_y_plus_x / (l_plus_y * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_x_plus_x * l_plus_y + nu_y_plus_x * l_plus_x) / (l_plus_y * l_plus_x * nu_l_l);
				}
				else { // -y和+y内部邻近点都不存在，作零梯度假设
					// values[elem_ptr_center] -= 2.0 / (l_plus_x * (l_minus_x + l_plus_x));
					b[row_idx] -= -2.0 * g_plus_x / nu_l_l;
					extra_minus_y += -2.0 * nu_y_plus_x / (l_minus_y * nu_l_l);
					values[elem_ptr_center] -= 2.0 * (nu_x_plus_x * l_minus_y - nu_y_plus_x * l_plus_x) / (l_minus_y * l_plus_x * nu_l_l);
				}
			}

			/*****************
			 * 处理角点额外项 *
			 ****************/

			if (extra_minus_x != 0.0) { 
				double nu_l_l = nu_x_minus_x * (l_minus_x + l_plus_x);
				double coeff = -2.0 / (l_minus_x * (l_minus_x + l_plus_x));
				double scale = extra_minus_x / coeff;
				if (has_minus_y) { // -y内部邻近点存在
					b[row_idx] -= scale * 2.0 * g_minus_x / nu_l_l;
					values[elem_ptr_minus_y] += scale * 2.0 * nu_y_minus_x / (l_minus_y * nu_l_l);
					values[elem_ptr_center] -= scale * 2.0 * (nu_x_minus_x * l_minus_y + nu_y_minus_x * l_minus_x) / (l_minus_y * l_minus_x * nu_l_l);
				}
				else if (has_plus_y) { // -y内部邻近点不存在，+y内部邻近点存在
					b[row_idx] -= scale * 2.0 * g_minus_x / nu_l_l;
					values[elem_ptr_plus_y] += scale * -2.0 * nu_y_minus_x / (l_plus_y * nu_l_l);
					values[elem_ptr_center] -= scale * 2.0 * (nu_x_minus_x * l_plus_y - nu_y_minus_x * l_minus_x) / (l_plus_y * l_minus_x * nu_l_l);
				}
				else { // 四个内部邻近点都不存在，抛出异常
					throw invalid_argument("网格剖分过于稀疏！");
				}
			}

			if (extra_minus_y != 0.0) {
				double nu_l_l = nu_y_minus_y * (l_minus_y + l_plus_y);
				double coeff = -2.0 / (l_minus_y * (l_minus_y + l_plus_y));
				double scale = extra_minus_y / coeff;
				if (has_minus_x) { // -x内部邻近点存在
					b[row_idx] -= scale * 2.0 * g_minus_y / nu_l_l;
					values[elem_ptr_minus_x] += scale * 2.0 * nu_x_minus_y / (l_minus_x * nu_l_l);
					values[elem_ptr_center] -= scale * 2.0 * (nu_y_minus_y * l_minus_x + nu_x_minus_y * l_minus_y) / (l_minus_x * l_minus_y * nu_l_l);
				}
				else if (has_plus_x) { // -x内部邻近点不存在，+x内部邻近点存在
					b[row_idx] -= scale * 2.0 * g_minus_y / nu_l_l;
					values[elem_ptr_plus_x] += scale * -2.0 * nu_x_minus_y / (l_plus_x * nu_l_l);
					values[elem_ptr_center] -= scale * 2.0 * (nu_y_minus_y * l_plus_x - nu_x_minus_y * l_minus_y) / (l_plus_x * l_minus_y * nu_l_l);
				}
				else { // 四个内部邻近点都不存在，抛出异常
					throw invalid_argument("网格剖分过于稀疏！");
				}
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
 * @brief 2维多边形区域Neumann边界Poisson问题求解程序
 */
int main() {
	cout << "当前运行的是 2维多边形区域Neumann边界Poisson问题 求解程序。" << endl;
	cout << "==============================================================" << endl;
	const string filepath = ".\\applications\\polygon_Neumann_poisson_solver\\";
	const string datapath = filepath + "data\\";
	if (CheckFolder(filepath, datapath)) return -1;
	cout << "文件路径：" << filepath << endl;
	cout << "数据路径：" << datapath << endl;
	
	cout << "1. 定义多边形区域" << endl;
	cout << "1.1. 逆时针顺序的顶点坐标：" << endl;
	vector<pair<double, double>> polygon_vertices = {
		{0.0, 1.0}, {0.0, -2.0}, {2.0, -2.0}, {1.0, -1.0}, {2.0, 1.0}, {1.0, 2.0}
	};
	// vector<pair<double, double>> polygon_vertices = {
	// 	{0.0, 0.0}, {2.0, 0.0}, {2.0, 1.0}, {0.0, 1.0}
	// };
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
	pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> grid_points_and_indices = GenerateInteriorGridWithAdaptiveBoundary(hx, hy, polygon_vertices);
	// pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> grid_points_and_indices = GenerateInteriorGrid(hx, hy, polygon_vertices);
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
	cout << "3.3. 定义Neumann边界条件 g(x,y)：" << endl;
	cout << "3.3.1. 计算区域边界的单位外法向量：" << endl;
	vector<pair<double, double>> boundary_normals = ComputePolygonBoundaryNormals(polygon_vertices);
	cout << "         已计算区域边界的单位外法向量：" << endl;
	cout << "         ";
	for (auto& n : boundary_normals) {
		cout << "(" << n.first << ", " << n.second << ") ";
	}
	cout << endl;
	cout << "3.3.2. 计算函数 u(x,y) 梯度：" << endl;
	auto dudx = [](double x, double y) { return 2.0 * M_PI * cos(2.0 * M_PI * x) * sin(M_PI * y); };
	auto dudy = [](double x, double y) { return M_PI * sin(2.0 * M_PI * x) * cos(M_PI * y); };
	cout << "         已定义函数 u(x,y) 梯度 dudx(x,y) = 2*pi*cos(2*pi*x)*sin(pi*y), dudy(x,y) = pi*sin(2*pi*x)*cos(pi*y)" << endl;
	cout << "3.3.3. 计算逐边界法向导数：" << endl;
	auto g = [ polygon_vertices, boundary_normals, dudx, dudy ](double x, double y)
	{ // 对每一个边界线段，利用函数 isAtSegment() 判断输入点(x,y)是否在边界线段上，如果在，则计算单位外法向量和梯度的内积；如果不在边界上，抛出异常。
		for (size_t i = 0; i < polygon_vertices.size(); i++) {
			if (isAtSegment(x, y, polygon_vertices[i], polygon_vertices[(i + 1) % polygon_vertices.size()])) {
				return dudx(x, y) * boundary_normals[i].first + dudy(x, y) * boundary_normals[i].second; // 外法向量
			}
		}
		throw invalid_argument("输入点(x,y)不在多边形区域内！");
	};
	cout << "	      已定义逐边界Neumann边界条件 g(x,y)." << endl;

	cout << "4. 构造矩阵 A 和右端项 b" << endl;
	pair<SparseMatrixCSR, vector<double>> linear_system = 
		ConstructLinearSystem(polygon_vertices, hx, hy, grid_points, grid_indices, f, g);
	cout << "4.1. 构造CSR稀疏矩阵：" << endl;
	cout << "4.1.1. 构造矩阵 A：" << endl;
	SparseMatrixCSR A = linear_system.first;
	cout << "         已构造矩阵 A，大小为 " << A.getRows() << " x " << A.getCols() << "，非零元数：" << A.getNNZ() << endl;
	cout << "4.1.2. 检验奇异性 A * 1 = 0：" << endl;
	vector<double> ones(grid_points_size, 1.0);
	vector<double> sigular = A * ones;
	double sigular_norm = norm(sigular);
	cout << "         向量 A * 1 的 L2 范数：" << sigular_norm << endl;
	cout << "4.1.3. 将矩阵 A 存入文件：" << endl;
	A.saveToFile(datapath + "matrix_A.coo");
	cout << "         已将矩阵 A 存入文件：" << datapath + "matrix_A.coo" << endl;
	cout << "4.2. 构造右端项 b：" << endl;
	vector<double> b = linear_system.second;
	cout << "       已构造右端项 b，大小为 " << b.size() << endl;

	cout << "5. 求解线性方程组 A U = b" << endl;
	cout << "5.1. 计算线性方程组的解 U_numeric：" << endl;
	vector<double> U_numeric = A.solve(b);
	cout << "       已计算线性方程组的解 U_numeric，大小为 " << U_numeric.size() << endl;
	cout << "5.2. 计算矫正系数 c：" << endl;
	vector<double> deltaU = U_numeric - U_exact;
	double c = -accumulate(deltaU.begin(), deltaU.end(), 0.0) / size(U_numeric);
	cout << "       已计算矫正系数 c = " << c << endl;
	cout << "5.3. 矫正 U_numeric 并计算求解精度：" << endl;
	U_numeric = U_numeric + c * ones;
	double epsilon = norm(A * U_numeric - b) / norm(b);
	cout << "       线性方程组的求解精度为" << epsilon << endl;
	cout << "5.4. 计算相对误差：" << endl;
	double error = norm(U_numeric - U_exact) / norm(U_exact);
	cout << "       在离散步长 hx = " << hx << " , hy = " << hy << " 下，相对误差为：" << error << endl;
	cout << "5.5. 将线性方程组的解 U_numeric 存入文件：" << endl;
	ofstream fileU_numeric(datapath + "U_numeric.txt");
	for (double u : U_numeric) {
		fileU_numeric << u << " ";
	}
	fileU_numeric.close();
	cout << "       已将线性方程组的解 U_numeric 存入文件：" << datapath + "U_numeric.txt" << endl;
	
	cout << "==============================================================" << endl;
	cout << "2维多边形区域Neumann边界Poisson问题求解程序 运行结束。" << endl;

	return 0;

	/**
	
	绘制多边形区域及网格： python .\applications\utils\plot_grid.py      polygon_neumann_poisson_solver
	绘制真解：            python .\applications\utils\plot_solution.py  polygon_neumann_poisson_solver grid_X grid_Y U_exact Neumann_exact_solution_504
	绘制矩阵 A 的稀疏结构：python .\applications\utils\spy_CSR_matrix.py polygon_neumann_poisson_solver matrix_A Neumann_polygon_cut_A_504
	绘制计算解：          python .\applications\utils\plot_solution.py  polygon_neumann_poisson_solver grid_X grid_Y U_numeric Neumann_cut_solution_504

	**/

	/**
	
	绘制多边形区域及网格： python .\applications\utils\plot_grid.py      polygon_neumann_poisson_solver
	绘制真解：            python .\applications\utils\plot_solution.py  polygon_neumann_poisson_solver grid_X grid_Y U_exact Neumann_exact_solution_171
	绘制矩阵 A 的稀疏结构：python .\applications\utils\spy_CSR_matrix.py polygon_neumann_poisson_solver matrix_A Neumann_rectangle_A_171
	绘制计算解：          python .\applications\utils\plot_solution.py  polygon_neumann_poisson_solver grid_X grid_Y U_numeric Neumann_rectangle_solution_171

	**/
}