#include "VectorOperations.hpp"
#include "SpecificMatrixGenerators.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <algorithm>

using namespace std;

namespace MatrixGenerators {

namespace Poisson2D {

// ======================= 工具函数实现 ======================= //

vector<double> flattenMatrix(const vector<vector<double>>& matrix_2d) {
	if (matrix_2d.empty()) return {};
	
	size_t M = matrix_2d.size();
	size_t N = matrix_2d[0].size();
	vector<double> flattened(M * N);
	
	for (size_t i = 0; i < M; ++i) {
		if (matrix_2d[i].size() != N) {
			throw invalid_argument("Input matrix has inconsistent row sizes");
		}
		for (size_t j = 0; j < N; ++j) {
			flattened[i * N + j] = matrix_2d[i][j];
		}
	}
	return flattened;
}

vector<vector<double>> unflattenVector(const vector<double>& vec_1d, size_t M, size_t N) {
	if (vec_1d.size() != M * N) {
		throw invalid_argument("Vector size does not match M * N");
	}
	
	vector<vector<double>> matrix_2d(M, vector<double>(N));
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			matrix_2d[i][j] = vec_1d[i * N + j];
		}
	}
	return matrix_2d;
}

// 二维索引到一维索引的映射
inline size_t index2Dto1D(size_t i, size_t j, size_t N) {
	return i * N + j;
}

// ======================= 矩阵生成函数实现 ======================= //

SparseMatrixCSR generatePoissonMatrix2D(const PoissonConfig2D& config) {
	size_t M = config.M - 1; // 网格数减1为实际内部点数
	size_t N = config.N - 1;
	size_t total_points = M * N; // 总未知格点数
	
	double hx = config.hx();
	double hy = config.hy();
	double hx2 = hx * hx;
	double hy2 = hy * hy;
	
	// 五点差分格式的系数(-Δ)
	double center_coeff = 2.0/hx2 + 2.0/hy2;
	double x_coeff = -1.0/hx2;
	double y_coeff = -1.0/hy2;
	
	vector<double> values;
	vector<size_t> col_indices;
	vector<size_t> row_ptrs(total_points + 1, 0);
	
	// 按行拉直：预先计算每行的非零元素个数
	for (size_t i = 0; i < M; ++i) { // i = 0:M-1
		for (size_t j = 0; j < N; ++j) { // j = 0:N-1
			size_t row_idx = index2Dto1D(i, j, N);
			size_t nnz_in_row = 0;
			
			// 检查-x邻居
			if (i > 0) nnz_in_row++;
			// 检查-y邻居
			if (j > 0) nnz_in_row++;
			// 中心点
			nnz_in_row++;
			// 检查+y邻居  
			if (j < N - 1) nnz_in_row++;
			// 检查+x邻居
			if (i < M - 1) nnz_in_row++;
			
			row_ptrs[row_idx + 1] = row_ptrs[row_idx] + nnz_in_row;
		}
	}
	
	// 调整values和col_indices的大小
	size_t total_nnz = row_ptrs[total_points];
	values.resize(total_nnz);
	col_indices.resize(total_nnz);
	
	// 填充矩阵元素
	for (size_t i = 0; i < M; ++i) { // i = 0:M-1
		for (size_t j = 0; j < N; ++j) { // j = 0:N-1
			size_t row_idx = index2Dto1D(i, j, N);
			size_t elem_ptr = row_ptrs[row_idx];
			
			/// 顺序是重要的！否则会导致CSR的同行元素未形成正确顺序
			// -x邻居 (i-1, j)
			if (i > 0) {
				size_t bottom_idx = index2Dto1D(i-1, j, N);
				values[elem_ptr] = y_coeff;
				col_indices[elem_ptr] = bottom_idx;
				elem_ptr++;
			}
			// -y邻居 (i, j-1)
			if (j > 0) {
				size_t left_idx = index2Dto1D(i, j-1, N);
				values[elem_ptr] = x_coeff;
				col_indices[elem_ptr] = left_idx;
				elem_ptr++;
			}
			// 中心元素
			values[elem_ptr] = center_coeff;
			col_indices[elem_ptr] = row_idx;
			elem_ptr++;
			// +y邻居 (i, j+1)
			if (j < N - 1) {
				size_t right_idx = index2Dto1D(i, j+1, N);
				values[elem_ptr] = x_coeff;
				col_indices[elem_ptr] = right_idx;
				elem_ptr++;
			}
			// +x邻居 (i+1, j)
			if (i < M - 1) {
				size_t top_idx = index2Dto1D(i+1, j, N);
				values[elem_ptr] = y_coeff;
				col_indices[elem_ptr] = top_idx;
				elem_ptr++;
			}
		}
	}
	
	return SparseMatrixCSR(total_points, total_points, values, col_indices, row_ptrs);
}

// ======================= 右端项生成函数实现 ======================= //

vector<double> generatePoissonRHS2D(
	const PoissonConfig2D& config,
	function<double(double, double)> f,
	function<double(double, double)> g) {
	
	size_t M = config.M - 1;
	size_t N = config.N - 1;
	double hx = config.hx();
	double hy = config.hy();
	double hx2 = hx * hx;
	double hy2 = hy * hy;
	
	vector<double> rhs(M * N, 0.0);
	
	// 计算网格点坐标并填充右端项
	for (size_t i = 0; i < M; ++i) { // i = 0:M-1
		double x = config.x_min + (i + 1) * hx; // 网格中心
		for (size_t j = 0; j < N; ++j) { // j = 0:N-1
			double y = config.y_min + (j + 1) * hy; // 网格中心
			size_t idx = index2Dto1D(i, j, N);

			// -x邻居 (i-1, j)
			if (i == 0) {
				rhs[idx] += g(x - hx, y) / hx2;
			}
			// -y邻居 (i, j-1)
			if (j == 0) {
				rhs[idx] += g(x, y - hy) / hy2;
			}
			// 中心元素：右端项为 f(x,y)
			rhs[idx] += f(x, y);
			// +y邻居 (i, j+1)
			if (j == N - 1) {
				rhs[idx] += g(x, y + hy) / hy2;
			}
			// +x邻居 (i+1, j)
			if (i == M - 1) {
				rhs[idx] += g(x + hx, y) / hx2;
			}
		}
	}
	
	return rhs;
}

vector<double> generatePoissonRHS2D(
	const PoissonConfig2D& config,
	const vector<vector<double>>& f_discrete,
	const vector<vector<double>>& g_discrete) {
	
	size_t M = config.M;
	size_t N = config.N;
	
	// 验证输入尺寸
	if (f_discrete.size() != M || f_discrete[0].size() != N) {
		throw invalid_argument("f_discrete dimensions do not match config M, N");
	}
	if (g_discrete.size() != M || g_discrete[0].size() != N) {
		throw invalid_argument("g_discrete dimensions do not match config M, N");
	}
	
	vector<double> rhs(M * N, 0.0);
	
	for (size_t i = 0; i < M; ++i) {
		for (size_t j = 0; j < N; ++j) {
			size_t idx = index2Dto1D(i, j, N);
			
			// 内部点使用f_discrete
			rhs[idx] = f_discrete[i][j];
			
			// 边界点使用g_discrete
			bool is_boundary = (i == 0) || (i == M - 1) || (j == 0) || (j == N - 1);
			if (is_boundary) {
				rhs[idx] = g_discrete[i][j];
			}
		}
	}
	
	return rhs;
}

// ======================= 完整问题生成函数实现 ======================= //

PoissonProblem2D generatePoisson2D(
	const PoissonConfig2D& config,
	function<double(double, double)> f,
	function<double(double, double)> g) {
	
	PoissonProblem2D problem;
	problem.config = config;
	
	// 生成矩阵和右端项
	problem.A = generatePoissonMatrix2D(config);
	problem.b = generatePoissonRHS2D(config, f, g);
	
	// 求解线性方程组
	problem.solution_1d = problem.A.solve(problem.b);
	problem.solution_2d = unflattenVector(problem.solution_1d, config.M - 1, config.N - 1);
	
	return problem;
}

PoissonProblem2D generatePoisson2D(
	const PoissonConfig2D& config,
	const vector<vector<double>>& f_discrete,
	const vector<vector<double>>& g_discrete) {
	
	PoissonProblem2D problem;
	problem.config = config;
	
	// 生成矩阵和右端项
	problem.A = generatePoissonMatrix2D(config);
	problem.b = generatePoissonRHS2D(config, f_discrete, g_discrete);
	
	return problem;
}

// ======================= PoissonProblem2D 成员函数实现 ======================= //

void PoissonProblem2D::saveToFiles(const string& basename) const {
	// 保存矩阵 A
	A.saveToFile(basename + "_matrix.coo");

	// 保存右端项 b
	VectorOps::writeVec(b, basename + "_rhs.vec");

	// 保存解 solution_2d
	VectorOps::writeMat(solution_2d, basename + "_solution.mat");

	// 保存网格
	vector<vector<double>> X(config.M - 1, vector<double>(config.N - 1));
	vector<vector<double>> Y(config.M - 1, vector<double>(config.N - 1));
	for (size_t i = 0; i < config.M - 1; ++i) {
		for (size_t j = 0; j < config.N - 1; ++j) {
			X[i][j] = config.x_min + (i + 1) * config.hx();
			Y[i][j] = config.y_min + (j + 1) * config.hy();
		}
	}
	VectorOps::writeMat(X, basename + "_X.mat");
	VectorOps::writeMat(Y, basename + "_Y.mat");
}

void PoissonProblem2D::printInfo() const {
	cout << "----------- 2D Poisson Problem Info -----------" << endl;
	cout << "Grid: " << config.M << " x " << config.N << " cells" << endl;
	cout << "Domain: [" << config.x_min << ", " << config.x_max << "] x [" 
		<< config.y_min << ", " << config.y_max << "]" << endl;
	cout << "Step sizes: hx = " << config.hx() << ", hy = " << config.hy() << endl;
	cout << "Matrix size: " << A.getRows() << " x " << A.getCols() << endl;
	cout << "Matrix NNZ: " << A.getNNZ() << endl;
	cout << "RHS vector size: " << b.size() << endl;
	cout << "Solution computed: " << (solution_2d.empty() ? "No" : "Yes") << endl;
	
	if (!solution_2d.empty()) {
		// 使用 std::min_element 和 std::max_element
		auto flat_solution = flattenMatrix(solution_2d);
		double min_val = *min_element(flat_solution.begin(), flat_solution.end());
		double max_val = *max_element(flat_solution.begin(), flat_solution.end());
		cout << "Solution range: [" << min_val << ", " << max_val << "]" << endl;
	}
	cout << "-----------------------------------------------" << endl;
}

} // namespace Poisson2D

} // namespace MatrixGenerators