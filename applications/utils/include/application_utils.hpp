#ifndef APPLICATION_UTILS_HPP
#define APPLICATION_UTILS_HPP

#include <iostream>
#include <string>
#include <filesystem>
#include <vector>
#include <utility>            // 用于 pair
#include <stdexcept>          // 异常处理
#include <cstddef>            // 用于 size_t
#include <cmath>              // 数学函数
#include <functional>         // 函数对象
#include <optional>           // 可选类型

using namespace std;

/**
 * @brief 检查文件夹是否存在，若不存在则创建
 * @param filepath 主文件路径
 * @param datapath 数据文件路径
 * @return 0 成功，-1 失败
 */
int CheckFolder(const string& filepath, const string& datapath);

/**
 * @brief 生成多边形区域内的网格点
 * @param hx 网格在 x 方向的步长
 * @param hy 网格在 y 方向的步长
 * @param polygon_vertices 多边形的顶点坐标
 * @return 网格点的坐标列表，前者为网格点的真实double类型坐标，后者为网格点的size_t类型坐标索引
 * @note 本函数使用 1e-10 的容差值处理浮点数精度问题
 */
pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> GenerateInteriorGrid(double hx, double hy, const vector<pair<double, double>>& polygon_vertices);

/**
 * @brief 生成多边形区域内的网格点，且自适应调整边界网格步长
 * @param hx 网格在 x 方向的步长
 * @param hy 网格在 y 方向的步长
 * @param polygon_vertices 多边形的顶点坐标
 * @return 网格点的坐标列表，前者为网格点的真实double类型坐标，后者为网格点的size_t类型坐标索引
 * @note 本函数使用 1e-10 的容差值处理浮点数精度问题
 */
pair<vector<pair<double, double>>, vector<pair<size_t, size_t>>> GenerateInteriorGridWithAdaptiveBoundary(double hx, double hy, const vector<pair<double, double>>& polygon_vertices);

/**
 * @brief 离散真解函数 u(x,y) 到网格点上的离散值 U(x_i,y_i)
 * @param u_exact 真解函数
 * @param grid_points 网格点的坐标列表
 * @return 网格点上离散值 U(x_i,y_i)
 */
vector<double> DiscretizeExactSolution(const function<double(double, double)>& u_exact, const vector<pair<double, double>>& grid_points);

/**
 * @brief 计算两条线段的交点坐标，以及交点到其中一个端点 A1 的距离
 * @param A_1 线段1的端点A坐标
 * @param B_1 线段1的端点B坐标
 * @param A_2 线段2的端点A坐标
 * @param B_2 线段2的端点B坐标
 * @return 两条线段的交点坐标，以及交点到A1的距离。如果不存在交点，则
 */
optional<pair<pair<double, double>, double>> segmentIntersection(pair<double, double> A_1, pair<double, double> B_1, pair<double, double> A_2, pair<double, double> B_2);

/**
 * @brief 计算非正则网格点的边界信息
 * @param polygon_vertices 多边形的顶点坐标
 * @param point 非正则网格点的坐标
 * @param boundary_type 边界类型，'left'、'bottom'、'top'、'right'
 * @param h 沿着 boundary_type 的网格步长
 * @param g Neumann边界条件函数
 * @return 边界点到网格点的距离、边界g值、边界单位外法向量
 */
tuple<double, double, pair<double, double>> irregularCalculator(const vector<pair<double, double>>& polygon_vertices, const pair<double, double>& point, const string& boundary_type, double h, const function<double(double, double)>& g);


#endif // APPLICATION_UTILS_HPP