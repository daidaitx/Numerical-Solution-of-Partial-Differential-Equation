#include "SparseMatrixCSR.hpp"              // CSR稀疏矩阵类
#include "VectorOperations.hpp"             // 稠密向量/矩阵操作
#include <iostream>
#include <vector>
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
 * @brief 检查文件夹是否存在，若不存在则创建
 * @param filepath 主文件路径
 * @param datapath 数据文件路径
 * @return 0 成功，-1 失败
 */
int CheckFolder(const string& filepath, const string& datapath) {
	// 检查路径是否存在
	namespace fs = std::filesystem;
	
	// 检查主文件路径是否存在
	if (!fs::exists(filepath)) {
		cerr << "错误: 文件路径不存在: " << filepath << endl;
		cerr << "请确保程序在正确的目录中运行。" << endl;
		return -1;
	}
	
	// 检查数据路径是否存在，如果不存在则创建
	if (!fs::exists(datapath)) {
		cout << "数据路径不存在，正在创建: " << datapath << endl;
		try {
			if (fs::create_directories(datapath)) {
				cout << "成功创建数据目录。" << endl;
			} else {
				cerr << "警告: 无法创建数据目录，但将继续运行。" << endl;
			}
		} catch (const fs::filesystem_error& e) {
			cerr << "错误: 创建数据目录失败: " << e.what() << endl;
			return -1;
		}
	}

	return 0;
}

/**
 * @brief 判断点 (x, y) 是否在多边形内部（不包括边界）
 * @param x 点的 x 坐标
 * @param y 点的 y 坐标
 * @param polygon_vertices 多边形的顶点坐标
 * @return true 点在多边形内部，false 点在多边形外部或边界上
 */
bool isInsidePolygon(const double x, const double y, const vector<pair<double, double>>& polygon_vertices) {
	int n = polygon_vertices.size();

	// 首先检查点是否在多边形的边界上
	for (int i = 0, j = n - 1; i < n; j = i++) {
		double x1 = polygon_vertices[j].first, y1 = polygon_vertices[j].second;
		double x2 = polygon_vertices[i].first, y2 = polygon_vertices[i].second;
		
		// 检查点是否在线段上
		if (fabs((y2 - y1) * (x - x1) - (y - y1) * (x2 - x1)) < 1e-10) {
			if (min(x1, x2) - 1e-10 <= x && x <= max(x1, x2) + 1e-10 &&
				min(y1, y2) - 1e-10 <= y && y <= max(y1, y2) + 1e-10) {
				return false; // 点在边界上，视为不在多边形内部
			}
		}
	}

	// 使用射线法判断内部/外部
	bool inside = false;
	for (int i = 0, j = n - 1; i < n; j = i++) {
		double xi = polygon_vertices[i].first, yi = polygon_vertices[i].second;
		double xj = polygon_vertices[j].first, yj = polygon_vertices[j].second;
		
		if (((yi > y) != (yj > y)) && 
			(x < (xj - xi) * (y - yi) / (yj - yi) + xi)) {
			inside = !inside;
		}
	}
	return inside;
}

/**
 * @brief 生成多边形区域内的网格点
 * @param hx 网格在 x 方向的步长
 * @param hy 网格在 y 方向的步长
 * @param polygon_vertices 多边形的顶点坐标
 * @return 网格点的坐标列表
 */
pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> GenerateGrid(double hx, double hy, const vector<pair<double, double>>& polygon_vertices) {
	vector<pair<double, double>> grid_points;
	vector<pair<size_t, size_t>> grid_indices;
	
	// 计算多边形的边界框
	double x_min = polygon_vertices[0].first;
	double x_max = polygon_vertices[0].first;
	double y_min = polygon_vertices[0].second;
	double y_max = polygon_vertices[0].second;
	for (const auto& vertex : polygon_vertices) {
		if (vertex.first < x_min) x_min = vertex.first;
		if (vertex.first > x_max) x_max = vertex.first;
		if (vertex.second < y_min) y_min = vertex.second;
		if (vertex.second > y_max) y_max = vertex.second;
	}

	// 生成网格点（真实坐标grid_points、整数坐标索引grid_indices）
	int num_x = static_cast<int>((x_max - x_min) / hx) + 1;
	int num_y = static_cast<int>((y_max - y_min) / hy) + 1;

	for (size_t i = 0; i < num_x; i++) {
		double x = x_min + i * hx;
		for (size_t j = 0; j < num_y; j++) {
			double y = y_min + j * hy;
			if (isInsidePolygon(x, y, polygon_vertices)) {
				grid_points.emplace_back(x, y);
				grid_indices.emplace_back(i, j);
			}
		}
	}

	// 如果网格点数量为0，则抛出异常
	if (grid_points.empty()) {
		throw runtime_error("网格点数量为0，可能需要缩小网格步长。");
	}

	pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> result = {grid_points, grid_indices};
	return result;
}

/**
 * @brief 离散真解函数 u(x,y) 到网格点上的离散值 U(x_i,y_i)
 * @param u_exact 真解函数
 * @param grid_points 网格点的坐标列表
 * @return 网格点上离散值 U(x_i,y_i)
 */
vector<double> DiscretizeExactSolution(const function<double(double, double)>& u_exact, const vector<pair<double, double>>& grid_points) {
	vector<double> U_exact(grid_points.size());
	for (size_t i = 0; i < grid_points.size(); i++) {
		U_exact[i] = u_exact(grid_points[i].first, grid_points[i].second);
	}
	return U_exact;
}

/**
 * @brief 计算两条线段的交点坐标，以及交点到其中一个端点 A1 的距离
 * @param A_1 线段1的端点A坐标
 * @param B_1 线段1的端点B坐标
 * @param A_2 线段2的端点A坐标
 * @param B_2 线段2的端点B坐标
 * @return 两条线段的交点坐标，以及交点到A1的距离。如果不存在交点，则
 */
optional<pair<pair<double, double>, double>> segmentIntersection(pair<double, double> A_1, pair<double, double> B_1, pair<double, double> A_2, pair<double, double> B_2) {
	// 检查端点是否重合
	if (A_1 == A_2 || A_1 == B_2) return make_pair(A_1, 0.0);
	if (B_1 == A_2 || B_1 == B_2) return make_pair(B_1, sqrt(pow(B_1.first - A_1.first, 2) + pow(B_1.second - A_1.second, 2)));
	
	// 计算线段1的向量
	double dx1 = B_1.first - A_1.first;
	double dy1 = B_1.second - A_1.second;

	// 计算线段2的向量
	double dx2 = B_2.first - A_2.first;
	double dy2 = B_2.second - A_2.second;

	// 计算分母
	double denom = dx1 * dy2 - dy1 * dx2;
	if (fabs(denom) < 1e-10) {
		// 线段平行或重合，无交点
		return nullopt;
	}

	// 计算参数方程的解t和u
	double t = ((A_2.first - A_1.first) * dy2 - (A_2.second - A_1.second) * dx2) / denom;
	double u = ((A_2.first - A_1.first) * dy1 - (A_2.second - A_1.second) * dx1) / denom;

	const double EPS = 1e-12;  // 更严格的容差
	// 在比较时使用容差
	if (t >= -EPS && t <= 1.0 + EPS && u >= -EPS && u <= 1.0 + EPS) {
		// 如果接近端点，进行修正
		if (t < 0.0) t = 0.0;
		if (t > 1.0) t = 1.0;
		
		double x_intersect = A_1.first + t * dx1;
		double y_intersect = A_1.second + t * dy1;
		double distance_to_A1 = sqrt(pow(x_intersect - A_1.first, 2) + pow(y_intersect - A_1.second, 2));

		return make_pair(make_pair(x_intersect, y_intersect), distance_to_A1);
	}
	else {
		// 交点不在两条线段上
		return nullopt;
	}
}

/**
 * @brief 计算非正则网格点的边界信息
 * @param polygon_vertices 多边形的顶点坐标
 * @param point 非正则网格点的坐标
 * @param boundary_type 边界类型，'left'、'bottom'、'top'、'right'
 * @param h 沿着 boundary_type 的网格步长
 * @param g Dirichlet边界条件函数
 * @return 边界点到网格点的距离和边界值
 */
pair<double, double> irregularCalculator(
	const vector<pair<double, double>>& polygon_vertices, 
	const pair<double, double>& point, 
	const string& boundary_type, 
	const double h, 
	const function<double(double, double)>& g)
	{
		// 初始化
		double l = 0.0; // 边界点到网格点的距离
		double g_value = 0.0; // 边界值
		pair<double, double> outside_point; // 区域外相邻点
		pair<double, double> boundary_point; // 边界点

		// 找到区域外相邻点
		if (boundary_type == "left") {
			outside_point = {point.first - h, point.second};
		} else if (boundary_type == "bottom") {
			outside_point = {point.first, point.second - h};
		} else if (boundary_type == "top") {
			outside_point = {point.first, point.second + h};
		} else if (boundary_type == "right") {
			outside_point = {point.first + h, point.second};
		} else {
			throw invalid_argument("无效的边界类型！");
		}

		// 寻找边界交点
		vector<pair<pair<double, double>, double>> boundary_points;
		for (size_t i = 0; i < polygon_vertices.size(); ++i) {
			pair<double, double> B1 = polygon_vertices[i];
			pair<double, double> B2 = polygon_vertices[(i + 1) % polygon_vertices.size()];
			optional<pair<pair<double, double>, double>> intersection_info = segmentIntersection(point, outside_point, B1, B2);
			if (intersection_info.has_value()) { // 存在交点
				boundary_points.push_back(intersection_info.value());
			}
		}

		// 为边界交点排序，选择距离最近的交点
		if (!boundary_points.empty()) {
			sort(boundary_points.begin(), boundary_points.end(), 
			[](const pair<pair<double, double>, double>& a, const pair<pair<double, double>, double>& b) {
				return a.second < b.second;
			});
			boundary_point = boundary_points[0].first;
			l = boundary_points[0].second;
		}
		else { // 无边界交点，靠近边界
			throw runtime_error("无法找到网格步长内的边界交点！");
		}

		// 计算边界值
		g_value = g(boundary_point.first, boundary_point.second);

		return make_pair(l, g_value);
	}

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

		// 初始化
		vector<double> b(total_points, 0.0);
		double hx2 = hx * hx;
		double hy2 = hy * hy;
	
		vector<double> values(total_points * 5);      // 预分配空间，最多每行5个非零元素
		vector<size_t> col_indices(total_points * 5); // 预分配空间，最多每行5个非零元素
		vector<size_t> row_ptrs(total_points + 1, 0);
	
		// 按x方向构造CSR稀疏矩阵 A 和右端项 b
		for (size_t row_idx = 0; row_idx < total_points; ++row_idx) { // row_idx 既是矩阵行索引，也是当前中心网格点的拉直索引
			
			// 重置 nnz_in_one_row ： 既是矩阵A的当前行的非零元素数，也是当前中心网格点的邻近网格点数，包括中心网格点自己，该值不超过5
			size_t nnz_in_one_row = 0;
			
			// 检查-x邻居
			size_t idx_minus_x = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first - 1, grid_indices[row_idx].second)) - grid_indices.begin(); // idx_minus_x 是 -x 邻居的拉直索引；如果不存在，则等于 total_points
			double l_minus_x = hx; // 到-x邻居的距离，默认为 hx，即规则网格距离
			double g_minus_x = 0.0; // -x邻居的边界值，默认为 0.0，即不贡献边界值
			if (idx_minus_x < total_points) { // 存在-x邻居
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在-x邻居，左侧靠近边界
				// 计算到左侧边界点的距离l_minus_x、左侧边界点的值g_minus_x
				pair<double, double> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "left", hx, g);
				l_minus_x = boundary_info.first;
				g_minus_x = boundary_info.second;
			}

			// 检查-y邻居
			size_t idx_minus_y = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first, grid_indices[row_idx].second - 1)) - grid_indices.begin(); // idx_minus_y 是 -y 邻居的拉直索引；如果不存在，则等于 total_points
			double l_minus_y = hy; // 到-y邻居的距离，默认为 hy，即规则网格距离
			double g_minus_y = 0.0; // -y邻居的边界值，默认为 0.0，即不贡献边界值
			if (idx_minus_y < total_points) { // 存在-y邻居
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在-y邻居，底部靠近边界
				// 计算到底部边界点的距离l_minus_y、底部边界点的值g_minus_y
				pair<double, double> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "bottom", hy, g);
				l_minus_y = boundary_info.first;
				g_minus_y = boundary_info.second;
			}

			// 中心点
			nnz_in_one_row++; // 一定有值，本行非零元素数++
			double f_center = f(grid_points[row_idx].first, grid_points[row_idx].second); // 中心点的源项值

			// 检查+y邻居
			size_t idx_plus_y = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first, grid_indices[row_idx].second + 1)) - grid_indices.begin(); // idx_plus_y 是 +y 邻居的拉直索引；如果不存在，则等于 total_points
			double l_plus_y = hy; // 到+y邻居的距离，默认为 hy，即规则网格距离
			double g_plus_y = 0.0; // +y邻居的边界值，默认为 0.0，即不贡献边界值
			if (idx_plus_y < total_points) { // 存在+y邻居
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在+y邻居，顶部靠近边界
				// 计算到顶部边界点的距离l_plus_y、顶部边界点的值g_plus_y
				pair<double, double> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "top", hy, g);
				l_plus_y = boundary_info.first;
				g_plus_y = boundary_info.second;
			}

			// 检查+x邻居
			size_t idx_plus_x = find(grid_indices.begin(), grid_indices.end(), make_pair(grid_indices[row_idx].first + 1, grid_indices[row_idx].second)) - grid_indices.begin(); // idx_plus_x 是 +x 邻居的拉直索引；如果不存在，则等于 total_points
			double l_plus_x = hx; // 到+x邻居的距离，默认为 hx，即规则网格距离
			double g_plus_x = 0.0; // +x邻居的边界值，默认为 0.0，即不贡献边界值
			if (idx_plus_x < total_points) { // 存在+x邻居
				nnz_in_one_row++; // 本行非零元素数++
			}
			else { // 不存在+x邻居，右侧靠近边界
				// 计算到右侧边界点的距离l_plus_x、右侧边界点的值g_plus_x
				pair<double, double> boundary_info = irregularCalculator(polygon_vertices, grid_points[row_idx], "right", hx, g);
				l_plus_x = boundary_info.first;
				g_plus_x = boundary_info.second;
			}
			
			// 递归更新 row_ptrs
			row_ptrs[row_idx + 1] = row_ptrs[row_idx] + nnz_in_one_row;
			
			// 更新 col_indices, values, b
			size_t elem_ptr = row_ptrs[row_idx];

			// -x邻居
			if (idx_minus_x < total_points) { // 存在-x邻居
				col_indices[elem_ptr] = idx_minus_x; // 列索引
				values[elem_ptr] = - 2.0 / (l_minus_x * (l_minus_x + l_plus_x)); // 元素值
				elem_ptr++;
			}
			else { // 不存在-x邻居，左侧靠近边界
				b[row_idx] += (2.0 * g_minus_x) / (l_minus_x * (l_minus_x + l_plus_x)); // 右端项贡献
			}

			// -y邻居
			if (idx_minus_y < total_points) { // 存在-y邻居
				col_indices[elem_ptr] = idx_minus_y; // 列索引
				values[elem_ptr] = - 2.0 / (l_minus_y * (l_minus_y + l_plus_y)); // 元素值
				elem_ptr++;
			}
			else { // 不存在-y邻居，底部靠近边界
				b[row_idx] += (2.0 * g_minus_y) / (l_minus_y * (l_minus_y + l_plus_y)); // 右端项贡献
			}

			// 中心元素
			col_indices[elem_ptr] = row_idx; // 列索引 = 行索引，即位于矩阵对角线位置
			values[elem_ptr] = 2.0 / (l_minus_x * l_plus_x) + 2.0 / (l_minus_y * l_plus_y); // 元素值
			elem_ptr++;
			b[row_idx] += f_center; // 右端项贡献

			// +y邻居
			if (idx_plus_y < total_points) { // 存在+y邻居
				col_indices[elem_ptr] = idx_plus_y; // 列索引
				values[elem_ptr] = - 2.0 / (l_plus_y * (l_minus_y + l_plus_y)); // 元素值
				elem_ptr++;
			}
			else { // 不存在+y邻居，顶部靠近边界
				b[row_idx] += (2.0 * g_plus_y) / (l_plus_y * (l_minus_y + l_plus_y)); // 右端项贡献
			}

			// +x邻居
			if (idx_plus_x < total_points) { // 存在+x邻居
				col_indices[elem_ptr] = idx_plus_x; // 列索引
				values[elem_ptr] = - 2.0 / (l_plus_x * (l_minus_x + l_plus_x)); // 元素值
				elem_ptr++;
			}
			else { // 不存在+x邻居，右侧靠近边界
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
	pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> grid_points_and_indices = GenerateGrid(hx, hy, polygon_vertices);
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

	cout << "==============================================================" << endl;
	cout << "绘图部分" << endl;
	cout << "2.4. 绘制多边形区域及网格：" << endl;
	// python .\applications\polygon_dirichlet_poisson_solver\plot_grid.py
	cout << "       已调用 python 绘制多边形区域。" << endl;
	cout << "3.1.4. 绘制真解：" << endl;
	// python .\applications\polygon_dirichlet_poisson_solver\plot_solution.py grid_X grid_Y U_exact exact_solution_506
	cout << "       已调用 python 绘制真解。" << endl;
	cout << "4.1.3. 绘制矩阵 A 的稀疏结构：" << endl;
	// python .\applications\polygon_dirichlet_poisson_solver\spy_CSR_sparse_matrix.py matrix_A.coo spy_matrix_A_506.png
	cout << "         已调用 python 绘制矩阵 A 的稀疏结构，保存为 spy_matrix_A_" << grid_points_size << ".png" << endl;
	cout << "5.4. 绘制计算解：" << endl;
	// python .\applications\polygon_dirichlet_poisson_solver\plot_solution.py grid_X grid_Y U_numeric numerical_solution_506
	cout << "       已调用 python 绘制计算解。" << endl;

	/**
	
	python .\applications\utils\plot_grid.py      polygon_dirichlet_poisson_solver
	python .\applications\utils\plot_solution.py  polygon_dirichlet_poisson_solver grid_X grid_Y U_exact exact_solution_506
	python .\applications\utils\spy_CSR_matrix.py polygon_dirichlet_poisson_solver matrix_A spy_matrix_A_506
	python .\applications\utils\plot_solution.py  polygon_dirichlet_poisson_solver grid_X grid_Y U_numeric numerical_solution_506

	**/
}