#include "SparseMatrixSolvers.hpp"
#include "SparseMatrixCSR.hpp"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include "VectorOperations.hpp"

using namespace std;
using namespace VectorOps;

namespace SparseMatrixSolvers {

// GMRES求解器实现
vector<double> solveGMRES(const SparseMatrixCSR& A,
							const vector<double>& b,
							const vector<double>& x0,
							size_t restart,
							size_t max_restarts,
							double tol) {
	// 基本检查
	if (!A.isSquare()) {
		throw invalid_argument("GMRES需要输入为方阵");
	}
	if (A.getRows() != b.size()) {
		throw invalid_argument("矩阵维度必须与向量大小匹配");
	}
	
	const size_t n = A.getRows();
	
	// 初始化解
	vector<double> x;
	if (x0.empty()) {
		x = vector<double>(n, 0.0);
	} else {
		if (x0.size() != n) {
			throw invalid_argument("初始解维度不匹配");
		}
		x = x0;
	}
	
	// 计算初始残差
	vector<double> r = b - A * x;
	double r_norm = norm(r);
	double b_norm = norm(b);
	
	// 如果初始残差已经足够小，直接返回
	if (r_norm < tol * max(1.0, b_norm)) {
		cout << "GMRES: 初始残差已经足够小，直接返回" << endl;
		return x;
	}
	
	// 主重启循环
	for (size_t restart_count = 0; restart_count < max_restarts; ++restart_count) {
		// Arnoldi过程初始化
		vector<vector<double>> V(restart + 1, vector<double>(n, 0.0));
		vector<double> g(restart + 1, 0.0);
		vector<vector<double>> H(restart + 1, vector<double>(restart, 0.0));
		
		// 第一个基向量 v0 = r / ||r||
		V[0] = r;
		for (double& val : V[0]) {
			val /= r_norm;
		}
		g[0] = r_norm;
		
		// Givens旋转存储
		vector<double> cs(restart, 0.0);  // cosine
		vector<double> sn(restart, 0.0);  // sine
		
		// Arnoldi迭代
		size_t j = 0;
		for (; j < restart; ++j) {
			// w = A * v_j
			vector<double> w = A * V[j];
			
			// 修正的Gram-Schmidt正交化
			for (size_t i = 0; i <= j; ++i) {
				H[i][j] = dot(V[i], w);
				for (size_t k = 0; k < n; ++k) {
					w[k] -= H[i][j] * V[i][k];
				}
			}
			
			H[j+1][j] = norm(w);
			
			// 检查breakdown
			if (abs(H[j+1][j]) < 1e-14) {
				// cout << "GMRES: Arnoldi breakdown at iteration " << j << endl;
				break;
			}
			
			// 新的基向量 v_{j+1} = w / H[j+1][j]
			V[j+1] = w;
			for (double& val : V[j+1]) {
				val /= H[j+1][j];
			}
			
			// 应用之前的Givens旋转到H的第j列
			for (size_t i = 0; i < j; ++i) {
				double temp = cs[i] * H[i][j] + sn[i] * H[i+1][j];
				H[i+1][j] = -sn[i] * H[i][j] + cs[i] * H[i+1][j];
				H[i][j] = temp;
			}
			
			// 计算新的Givens旋转
			double h1 = H[j][j];
			double h2 = H[j+1][j];
			double nu = sqrt(h1 * h1 + h2 * h2);
			
			if (abs(nu) < 1e-14) {
				nu = 1e-14;
			}
			
			cs[j] = h1 / nu;
			sn[j] = h2 / nu;
			
			// 应用Givens旋转到H和g
			H[j][j] = nu;
			H[j+1][j] = 0.0;
			
			// 更新g
			double temp = cs[j] * g[j] + sn[j] * g[j+1];
			g[j+1] = -sn[j] * g[j] + cs[j] * g[j+1];
			g[j] = temp;
			
			// 检查收敛
			double residual = abs(g[j+1]);
			if (residual < tol * max(1.0, b_norm)) {
				j++;  // 当前步已收敛
				break;
			}
		}
		
		// 解的索引范围
		size_t k = j;  // 实际使用的Krylov子空间维数
		
		// 回代求解上三角系统 H(1:k,1:k) * y = g(1:k)
		vector<double> y(k, 0.0);
		for (int i = k-1; i >= 0; --i) {
			y[i] = g[i];
			for (size_t m = i+1; m < k; ++m) {
				y[i] -= H[i][m] * y[m];
			}
			y[i] /= H[i][i];
		}
		
		// 更新解 x = x + V(1:k) * y
		for (size_t i = 0; i < k; ++i) {
			for (size_t m = 0; m < n; ++m) {
				x[m] += V[i][m] * y[i];
			}
		}
		
		// 计算新的残差
		r = b - A * x;
		r_norm = norm(r);
		
		// 输出重启信息
		cout << "GMRES restart " << restart_count + 1 
				<< ", residual norm: " << r_norm << endl;
		
		// 检查最终收敛
		if (r_norm < tol * max(1.0, b_norm)) {
			cout << "GMRES converged after " << (restart_count + 1) << " restarts" << endl;
			return x;
		}
		
		// 如果Arnoldi提前终止且残差没有改善，可能无法继续收敛
		if (j < restart && r_norm >= 0.9 * norm(b - A * vector<double>(n, 0.0))) {
			cout << "GMRES: Stagnation detected, stopping early" << endl;
			break;
		}
	}
	
	cout << "GMRES reached maximum restarts (" << max_restarts 
			<< ") without convergence. Final residual: " << r_norm << endl;
	return x;
}

} // namespace SparseMatrixSolvers