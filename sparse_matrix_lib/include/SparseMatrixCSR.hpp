#pragma once
#include <vector>
#include <cstddef>
#include <iostream>

/**
 * @class SparseMatrixCSR
 * @brief 行压缩稀疏矩阵（CSR格式）实现
 * @details 介绍什么是CSR：
 * CSR格式使用三个一维数组存储稀疏矩阵：
 * - values: 非零元素值数组，存储非零元素的数值
 * - col_indices: 列索引数组，存储每个非零元素对应的列索引
 * - row_ptrs: 行指针数组，存储每个非零元素所在的行的起始位置
 * 
 * 示例：
 *     | 1 1 0 4 |
 * A = | 0 0 0 5 |
 *     | 0 0 0 0 |
 *     | 1 0 4 0 |
 * 
 * 对应的CSR格式为：
 * values      = [1, 1, 4, 5, 1, 4];
 * col_indices = [0, 1, 3, 3, 0, 2];
 * row_ptrs    = [0, 3, 4, 4, 6];
 * 
 */
class SparseMatrixCSR {
private:
	size_t rows;      ///< 行数
	size_t cols;      ///< 列数
	std::vector<double> values;      ///< 非零元素值
	std::vector<size_t> col_indices; ///< 列索引
	std::vector<size_t> row_ptrs;    ///< 行指针

public:
	// ============================================================= //
	// ======================= 构造与析构函数 ======================= //
	// ============================================================= //

	/**
	 * @brief 构造函数：默认构造函数
	 */
	SparseMatrixCSR() = default;

	/**
	 * @brief 构造函数：构造一个 n x n 全零矩阵
	 * @param n 矩阵边长
	 */
	explicit SparseMatrixCSR(const size_t n) : SparseMatrixCSR(n, n) {};

	/**
	 * @brief 构造函数：构造一个 r x c 全零矩阵
	 * @param r 行数
	 * @param c 列数
	 */
	explicit SparseMatrixCSR(const size_t r, const size_t c);

	/**
	 * @brief 构造函数：稠密矩阵转稀疏矩阵
	 * @param dense 稠密矩阵
	 */
	explicit SparseMatrixCSR(const std::vector<std::vector<double>>& dense);

	/**
	 * @brief 构造函数：从COO格式矩阵构造CSR格式矩阵
	 * @param rows 矩阵行数
	 * @param cols 矩阵列数
	 * @param row_indices 非零元素行索引
	 * @param col_indices 非零元素列索引
	 * @param values 非零元素值
	 */
	explicit SparseMatrixCSR(const size_t rows, const size_t cols,
		const std::vector<size_t>& row_indices,
		const std::vector<size_t>& col_indices,
		const std::vector<double>& values);

	/**
	 * @brief 构造函数：直接从CSR格式构造CSR格式矩阵
	 * @param rows 矩阵行数
	 * @param cols 矩阵列数
	 * @param values 非零元素值
	 * @param col_indices 非零元素列索引
	 * @param row_ptrs 非零元素所在行的起始位置
	 */
	explicit SparseMatrixCSR(const size_t rows, const size_t cols,
		const std::vector<double>& values,
		const std::vector<size_t>& col_indices,
		const std::vector<size_t>& row_ptrs);

	/**
	 * @brief 构造函数：从文件构造矩阵
	 * @param filename 文件名
	 */
	explicit SparseMatrixCSR(const std::string& filename);

	/**
	 * @brief 析构函数
	 */
	~SparseMatrixCSR() = default;

	// ========================================================= //
	// ======================= 运算符重载 ======================= //
	// ========================================================= //

	/**
	 * @brief 运算符重载：拷贝赋值
	 * @param other 待拷贝的矩阵
	 * @return 赋值后的矩阵
	 */
	SparseMatrixCSR& operator=(const SparseMatrixCSR& other) = default;

	/**
	 * @brief 运算符重载：移动赋值
	 * @param other 待移动的矩阵
	 * @return 移动后的矩阵
	 */
	SparseMatrixCSR& operator=(SparseMatrixCSR&& other) noexcept = default;

	/**
	 * @brief 运算符重载：判断严格相等
	 * @param other 待比较的矩阵
	 * @return 是否严格相等
	 */
	bool operator==(const SparseMatrixCSR& other) const { return this->isApproxEqualto(other, 0); };

	/**
	 * @brief 运算符重载：判断严格不等
	 * @param other 待比较的矩阵
	 * @return 是否严格不等
	 */
	bool operator!=(const SparseMatrixCSR& other) const { return this->isNotApproxEqualto(other, 0); };

	/**
	 * @brief 运算符重载：矩阵累加
	 * @param other 加数矩阵
	 * @return 矩阵累加结果
	 */
	SparseMatrixCSR& operator+=(const SparseMatrixCSR& other);

	/**
	 * @brief 运算符重载：矩阵加法
	 * @param other 加数矩阵
	 * @return 矩阵相加结果
	 */
	SparseMatrixCSR operator+(const SparseMatrixCSR& other) const;

	/**
	 * @brief 运算符重载：矩阵累减
	 * @param other 减数矩阵
	 * @return 矩阵累减结果
	 */
	SparseMatrixCSR& operator-=(const SparseMatrixCSR& other);

	/**
	 * @brief 运算符重载：矩阵减法
	 * @param other 减数矩阵
	 * @return 矩阵相减结果
	 */
	SparseMatrixCSR operator-(const SparseMatrixCSR& other) const;

	/**
	 * @brief 运算符重载：矩阵乘法（稠密矩阵）
	 * @param other 乘数矩阵（稠密矩阵）
	 * @return 矩阵相乘结果
	 */
	std::vector<std::vector<double>> operator*(const std::vector<std::vector<double>>& other) const;

	/**
	 * @brief 运算符重载：矩阵乘向量
	 * @param other 向量
	 * @return 矩阵乘向量结果
	 */
	std::vector<double> operator*(const std::vector<double>& other) const;

	/**
	 * @brief 运算符重载：矩阵数乘
	 * @param scalar 数值
	 * @return 矩阵数乘结果
	 */
	SparseMatrixCSR operator*(const double scalar) const noexcept;

	/**
	 * @brief 运算符重载：矩阵数乘（左侧）
	 * @param scalar 数值
	 * @return 矩阵数乘结果
	 */
	friend SparseMatrixCSR operator*(const double scalar, const SparseMatrixCSR& mat) noexcept { return mat * scalar; };

	/**
	 * @brief 运算符重载：矩阵数除
	 * @param scalar 数值
	 * @return 矩阵数除结果
	 */
	SparseMatrixCSR operator/(const double scalar) const;

	/**
	 * @brief 运算符重载：单目负号
	 * @return 逐元素相反数的矩阵
	 */
	SparseMatrixCSR operator-() const;

	/**
	 * @brief 运算符重载：只读访问矩阵元素
	 * @param row 行索引
	 * @param col 列索引
	 * @return 矩阵元素值
	 */
	double operator()(const size_t row, const size_t col) const;

	/**
	 * @brief 运算符重载：流输出（左移）
	 * @param os 输出流
	 * @param mat 矩阵
	 * @return 输出流
	 */
	friend std::ostream& operator<<(std::ostream& os, const SparseMatrixCSR& mat);

	// ======================================================= //
	// ======================= 功能函数 ======================= //
	// ======================================================= //

	/**
	 * @brief 深拷贝
	 * @param other 待拷贝的矩阵
	 */
	SparseMatrixCSR(const SparseMatrixCSR& other) = default;

	/**
	 * @brief 移动构造函数
	 * @param other 待移动的矩阵
	 */
	SparseMatrixCSR(SparseMatrixCSR&& other) noexcept = default;

	/**
	 * @brief 获取矩阵行数
	 * @return 矩阵行数
	 */
	size_t getRows() const noexcept;

	/**
	 * @brief 获取矩阵列数
	 * @return 矩阵列数
	 */
	size_t getCols() const noexcept;

	/**
	 * @brief 获取矩阵大小（行数和列数）
	 * @return 矩阵大小
	 */
	std::pair<size_t, size_t> getSize() const noexcept;

	/**
	 * @brief 获取非零元素数目
	 * @return 非零元素数目
	 */
	size_t getNNZ() const noexcept;

	/**
	 * @brief 更新/插入/置零元素
	 * @param row 行索引
	 * @param col 列索引
	 * @param value 元素值
	 */
	void setValue(const size_t row, const size_t col, const double value);

	/**
	 * @brief 更新/插入/置零元素
	 * @param elements 元素三元组数组（行索引、列索引、值）
	 */
	void setValue(const std::tuple<size_t, size_t, double> element);

	/**
	 * @brief 批量更新/插入/置零元素
	 * @param row_indices 行索引数组
	 * @param col_indices 列索引数组
	 * @param values 元素值数组
	 */
	void setValue(const std::vector<size_t>& row_indices,
		const std::vector<size_t>& col_indices,
		const std::vector<double>& values);

	/**
	 * @brief 批量更新/插入/置零元素
	 * @param elements 元素三元组数组（行索引、列索引、值）
	 */
	void setValue(std::vector<std::tuple<size_t, size_t, double>> elements);

	/**
	 * @brief 判断近似相等
	 * @param other 待比较的矩阵
	 * @param epsilon 容许误差
	 * @return 是否近似相等
	 */
	bool isApproxEqualto(const SparseMatrixCSR& other, const double tolerance) const;

	/**
	 * @brief 判断近似不等
	 * @param other 待比较的矩阵
	 * @param epsilon 容许误差
	 * @return 是否近似不等
	 */
	bool isNotApproxEqualto(const SparseMatrixCSR& other, const double tolerance) const {
		return !this->isApproxEqualto(other, tolerance);
	};

	/**
	 * @brief 判断非常不等
	 * @param other 待比较的矩阵
	 * @param epsilon 容许误差
	 * @return 是否近似不等
	 */
	bool isFarAwayFrom(const SparseMatrixCSR& other, const double tolerance) const {
		return this->isNotApproxEqualto(other, tolerance);
	};

	/**
	 * @brief 计算矩阵的范数
	 * @param norm_type 范数类型
	 * @return 矩阵的范数
	 */
	double norm(const std::string& norm_type) const;

	/**
	 * @brief 矩阵转置
	 * @return 转置后的矩阵
	 */
	SparseMatrixCSR transpose() const noexcept;

	/**
	 * @brief 矩阵转置（简写）
	 * @return 转置后的矩阵
	 */
	SparseMatrixCSR t() const noexcept { return this->transpose(); };

	/**
	 * @brief 转为稠密矩阵
	 * @return 稠密矩阵
	 */
	std::vector<std::vector<double>> toDense() const noexcept;

	/**
	 * @brief 获取对角线元素（以向量形式返回）
	 * @return 对角线元素构成的向量
	 */
	std::vector<double> getDiagonalVector() const noexcept;

	/**
	 * @brief 获取对角线元素（以CSR格式的对角矩阵形式返回）
	 * @return 对角线元素构成的CSR格式的对角矩阵
	 */
	SparseMatrixCSR getDiagonalMatrix() const noexcept;

	/**
	 * @brief 获取上三角矩阵（包含对角线）
	 * @return 上三角矩阵
	 */
	SparseMatrixCSR getUpperTriangularMatrix() const;

	/**
	 * @brief 获取下三角矩阵（包含对角线）
	 * @return 下三角矩阵
	 */
	SparseMatrixCSR getLowerTriangularMatrix() const;

	/**
	 * @brief 判断矩阵是否为方阵
	 * @return 是否为方阵
	 */
	bool isSquare() const noexcept;

	/**
	 * @brief 判断矩阵是否为对称矩阵
	 * @return 是否为对称矩阵
	 */
	bool isSymmetric(const double tolerance = 0) const;

	/**
	 * @brief 判断矩阵是否为上三角矩阵
	 * @return 是否为上三角矩阵
	 */
	bool isUpperTriangular() const noexcept;

	/**
	 * @brief 判断矩阵是否为下三角矩阵
	 * @return 是否为下三角矩阵
	 */
	bool isLowerTriangular() const noexcept;

	/**
	 * @brief 判断矩阵是否为对角矩阵
	 * @return 是否为对角矩阵
	 */
	bool isDiagonal() const noexcept;

	/**
	 * @brief 矩阵保存到文件
	 * @param filename 文件名
	 */
	void saveToFile(const std::string& filename) const;

	/**
	 * @brief 矩阵打印
	 */
	void print() const noexcept { std::cout << *this <<std::endl; };

	/**
	 * @brief 矩阵稠密打印
	 * @param precision 精度，默认4位
	 */
	void printDense(const int precision = 4) const noexcept;

	/**
	 * @brief 求解线性方程组 Ax = b（GMRES求解器）
	 * @param b 右端项向量
	 * @return 解向量
	 */
	std::vector<double> solve(const std::vector<double>& b) const;
};