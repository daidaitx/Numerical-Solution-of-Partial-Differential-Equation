#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <stdexcept>

/**
 * @brief 向量运算工具库命名空间
 */
namespace VectorOps {

// ======================= 基本运算 ======================= //

/**
 * @brief 计算向量范数
 * @param v 输入向量
 * @param normType 范数类型，"1"表示1-范数，"2"表示2-范数，"inf"表示无穷范数，默认为"2"
 * @return 向量范数
 */
double norm(const std::vector<double>& v, const std::string& normType = "2");

/**
 * @brief 计算向量点积
 * @param a 第一个向量
 * @param b 第二个向量
 * @return 点积结果
 */
double dot(const std::vector<double>& a, const std::vector<double>& b);

// ======================= 算术运算符 ======================= //

/**
 * @brief 向量加法
 * @param a 第一个向量
 * @param b 第二个向量
 * @return 向量加法结果
 */
std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b);

/**
 * @brief 向量减法
 * @param a 第一个向量
 * @param b 第二个向量
 * @return 向量减法结果
 */
std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b);

/**
 * @brief 向量数乘（右侧）
 * @param v 向量
 * @param scalar 标量
 * @return 向量数乘结果
 */
std::vector<double> operator*(const std::vector<double>& v, double scalar);

/**
 * @brief 向量数乘（左侧）
 * @param scalar 标量
 * @param v 向量
 * @return 向量数乘结果
 */
std::vector<double> operator*(double scalar, const std::vector<double>& v);

/**
 * @brief 向量数除
 * @param v 向量
 * @param scalar 标量
 * @return 向量数除结果
 */
std::vector<double> operator/(const std::vector<double>& v, double scalar);

/**
 * @brief 向量单目负
 * @param v 向量
 * @return 向量单目负结果
 */
std::vector<double> operator-(const std::vector<double>& v);

// ======================= 复合赋值运算符 ======================= //

/**
 * @brief 向量累加
 * @param a 第一个向量
 * @param b 第二个向量
 * @return 向量累加结果
 */
std::vector<double>& operator+=(std::vector<double>& a, const std::vector<double>& b);

/**
 * @brief 向量累减
 * @param a 第一个向量
 * @param b 第二个向量
 * @return 向量累减结果
 */
std::vector<double>& operator-=(std::vector<double>& a, const std::vector<double>& b);

/**
 * @brief 向量数乘赋值
 * @param v 向量
 * @param scalar 标量
 * @return 向量数乘赋值结果
 */
std::vector<double>& operator*=(std::vector<double>& v, double scalar);

/**
 * @brief 向量数除赋值
 * @param v 向量
 * @param scalar 标量
 * @return 向量数除赋值结果
 */
std::vector<double>& operator/=(std::vector<double>& v, double scalar);

// ======================= 输出和工具函数 ======================= //

/**
 * @brief 向量输出到流
 * @param os 输出流
 * @param v 向量
 * @return 输出流
 */
std::ostream& operator<<(std::ostream& os, const std::vector<double>& v);

/**
 * @brief 打印向量
 * @param v 要打印的向量
 * @param asColumnVector 是否按列向量打印（默认true）
 * @param precision 小数精度（默认4位）
 */
void print(const std::vector<double>& v, int precision = 4, bool asColumnVector = true);

/**
 * @brief 打印稠密矩阵
 * @param matrix 要打印的矩阵（二维向量）
 * @param precision 小数精度（默认4位）
 */
void print(const std::vector<std::vector<double>>& matrix, int precision = 4);

/**
 * @brief 创建全零向量
 * @param n 向量维数
 * @return 全零向量
 */
std::vector<double> zeros(size_t n);

/**
 * @brief 创建全一向量
 * @param n 向量维数
 * @return 全一向量
 */
std::vector<double> ones(size_t n);

/**
 * @brief 创建随机向量
 * @param n 向量维数
 * @param min 最小值
 * @param max 最大值
 * @return 随机向量
 */
std::vector<double> random(size_t n, double min = 0.0, double max = 1.0);

/**
 * @brief 向量是否近似相等
 * @param a 第一个向量
 * @param b 第二个向量
 * @param tol 容差
 * @return 相等返回true，否则返回false
 */
bool areApproxEqual(const std::vector<double>& a, const std::vector<double>& b, double tol = 1e-12);

/**
 * @brief 向量是否绝对相等
 * @param a 第一个向量
 * @param b 第二个向量
 * @return 相等返回true，否则返回false
 */
inline bool areEqual(const std::vector<double>& a, const std::vector<double>& b) {
	return areApproxEqual(a, b, 0.0);
}

} // namespace VectorOps