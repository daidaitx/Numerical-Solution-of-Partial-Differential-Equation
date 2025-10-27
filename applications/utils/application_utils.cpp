#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
#include <utility>            // 用于 pair
#include <stdexcept>          // 异常处理
#include <cstddef>            // 用于 size_t
#include <cmath>              // 数学函数
#include <functional>         // 函数对象
#include <algorithm>          // 算法
#include <optional>           // 可选类型

using namespace std;

int CheckFolder(const string& filepath, const string& datapath) {
	// 检查主文件路径是否存在
	if (!filesystem::exists(filepath)) {
		cerr << "错误: 文件路径不存在: " << filepath << endl;
		cerr << "请确保程序在正确的目录中运行。" << endl;
		return -1;
	}
	
	// 检查数据路径是否存在，如果不存在则创建
	if (!filesystem::exists(datapath)) {
		cout << "数据路径不存在，正在创建: " << datapath << endl;
		try {
			if (filesystem::create_directories(datapath)) {
				cout << "成功创建数据目录。" << endl;
			} else {
				cerr << "警告: 无法创建数据目录，但将继续运行。" << endl;
			}
		} catch (const filesystem::filesystem_error& e) {
			cerr << "错误: 创建数据目录失败: " << e.what() << endl;
			return -1;
		}
	}

	return 0;
}

namespace { // 内部命名空间，该函数仅在本文件中使用
	/**
	 * @brief 判断点 (x, y) 是否在多边形内部（不包括边界）
	 * @param x 点的 x 坐标
	 * @param y 点的 y 坐标
	 * @param polygon_vertices 多边形的顶点坐标
	 * @return true 点在多边形内部，false 点在多边形外部或边界上
	 * @note 本函数使用 1e-10 的容差值处理浮点数精度问题
	 * @internal
	 * @attention 这是内部辅助函数，只在当前源文件中使用
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
}

pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> GenerateInteriorGrid(double hx, double hy, const vector<pair<double, double>>& polygon_vertices) {
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

namespace {
/**
 * @brief 计算非正则内点的阶数
 * @param grid_indices 网格索引向量
 * @param index 中心点拉直索引
 * @return 非正则内点的阶数
 */
int NonRegularOrder(const vector<pair<size_t, size_t>>& grid_indices, size_t index) {
	int count = 0;
	size_t i = grid_indices[index].first;
	size_t j = grid_indices[index].second;
	if(find(grid_indices.begin(), grid_indices.end(), make_pair(i-1, j)) == grid_indices.end()) count++;
	if(find(grid_indices.begin(), grid_indices.end(), make_pair(i+1, j)) == grid_indices.end()) count++;
	if(find(grid_indices.begin(), grid_indices.end(), make_pair(i, j-1)) == grid_indices.end()) count++;
	if(find(grid_indices.begin(), grid_indices.end(), make_pair(i, j+1)) == grid_indices.end()) count++;
	return count;
}
}

pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> GenerateInteriorGridWithAdaptiveBoundary(double hx, double hy, const vector<pair<double, double>>& polygon_vertices) {
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

	// 循环删去所有3、4阶非正则奇点
	bool change = true;
	while (change) { 
		change = false;
		for (int i = 0; i < grid_points.size(); ++i) {
			if (NonRegularOrder(grid_indices, i) >= 3) {
				change = true;
				grid_points.erase(grid_points.begin() + i);
				grid_indices.erase(grid_indices.begin() + i);
				--i;
			}
		}
	}

	pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> result = {grid_points, grid_indices};
	return result;
}

vector<double> DiscretizeExactSolution(const function<double(double, double)>& u_exact, const vector<pair<double, double>>& grid_points) {
	vector<double> U_exact(grid_points.size());
	for (size_t i = 0; i < grid_points.size(); i++) {
		U_exact[i] = u_exact(grid_points[i].first, grid_points[i].second);
	}
	return U_exact;
}

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

tuple<double, double, pair<double, double>> irregularCalculator(
	const vector<pair<double, double>>& polygon_vertices, 
	const pair<double, double>& point, 
	const string& boundary_type, 
	const double h, 
	const function<double(double, double)>& g)
	{
		// 初始化
		double l = 0.0; // 边界邻交点到中心点的距离
		double g_value = 0.0; // 边界g值
		pair<double, double> normal; // 边界单位外法向量
		pair<double, double> outside_point; // 外部邻近点
		pair<double, double> boundary_point; // 边界邻交点

		// 找到外部邻近点（方向）
		double kh = h * 10;
		if (boundary_type == "left") {
			outside_point = {point.first - kh, point.second};
		} else if (boundary_type == "bottom") {
			outside_point = {point.first, point.second - kh};
		} else if (boundary_type == "top") {
			outside_point = {point.first, point.second + kh};
		} else if (boundary_type == "right") {
			outside_point = {point.first + kh, point.second};
		} else {
			throw invalid_argument("无效的边界类型！");
		}

		// 寻找边界邻交点
		vector<tuple<pair<double, double>, double, pair<double, double>>> boundary_points; // 边界邻交点坐标，到边界邻交点的距离，所属边界的单位外法向量
		for (size_t i = 0; i < polygon_vertices.size(); ++i) {
			pair<double, double> B1 = polygon_vertices[i];
			pair<double, double> B2 = polygon_vertices[(i + 1) % polygon_vertices.size()];
			optional<pair<pair<double, double>, double>> intersection_info = segmentIntersection(point, outside_point, B1, B2);
			if (intersection_info.has_value()) { // 存在交点
				// 计算当前边界的单位外法向量
				double normal_x = (B2.second - B1.second);
				double normal_y = -(B2.first - B1.first);
				double normal_norm = sqrt(normal_x * normal_x + normal_y * normal_y);
				pair<double, double> normal = {normal_x / normal_norm, normal_y / normal_norm};
				boundary_points.push_back(make_tuple(intersection_info.value().first, intersection_info.value().second, normal));
			}
		}
		// 为边界邻交点排序，选择距离最近的交点
		if (!boundary_points.empty()) {
			sort(boundary_points.begin(), boundary_points.end(), 
				[](const tuple<pair<double, double>, double, pair<double, double>>& a, const tuple<pair<double, double>, double, pair<double, double>>& b) {
					return get<1>(a) < get<1>(b);
				});
			l = get<1>(boundary_points[0]);
			boundary_point = get<0>(boundary_points[0]);
			g_value = g(boundary_point.first, boundary_point.second);
			normal = get<2>(boundary_points[0]);
		}
		else { // 无边界邻交点
			throw runtime_error("无法找到网格步长内的边界邻交点！");
		}

		return make_tuple(l, g_value, normal);
	}


