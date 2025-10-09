#include "VectorOperations.hpp"
#include <string>
#include <random>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;

namespace VectorOps {

// ======================= 基本运算实现 ======================= //

double norm(const vector<double>& v, const string& normType) {
	if (normType == "2") {
		double sum = 0.0;
		for (double val : v) {
			sum += val * val;
		}
		return sqrt(sum);
	} else if (normType == "1") {
		double sum = 0.0;
		for (double val : v) {
			sum += abs(val);
		}
		return sum;
	} else if (normType == "inf") {
		if (v.empty()) return 0.0;
		double max_val = abs(v[0]);
		for (size_t i = 1; i < v.size(); ++i) {
			max_val = max(max_val, abs(v[i]));
		}
		return max_val;
	} else {
		throw invalid_argument("不支持的范数类型: " + normType);
	}
}

double dot(const vector<double>& a, const vector<double>& b) {
	if (a.size() != b.size()) {
		throw invalid_argument("计算点积时，向量大小必须相同");
	}
	
	double result = 0.0;
	for (size_t i = 0; i < a.size(); ++i) {
		result += a[i] * b[i];
	}
	return result;
}

// ======================= 算术运算符实现 ======================= //

vector<double> operator+(const vector<double>& a, const vector<double>& b) {
	if (a.size() != b.size()) {
		throw invalid_argument("向量大小必须相同");
	}
	
	vector<double> result(a.size());
	for (size_t i = 0; i < a.size(); ++i) {
		result[i] = a[i] + b[i];
	}
	return result;
}

vector<double> operator-(const vector<double>& a, const vector<double>& b) {
	if (a.size() != b.size()) {
		throw invalid_argument("向量大小必须相同");
	}
	
	vector<double> result(a.size());
	for (size_t i = 0; i < a.size(); ++i) {
		result[i] = a[i] - b[i];
	}
	return result;
}

vector<double> operator*(const vector<double>& v, double scalar) {
	vector<double> result(v.size());
	for (size_t i = 0; i < v.size(); ++i) {
		result[i] = v[i] * scalar;
	}
	return result;
}

vector<double> operator*(double scalar, const vector<double>& v) {
	return v * scalar;  // 复用上面的实现
}

vector<double> operator/(const vector<double>& v, double scalar) {
	if (abs(scalar) < 1e-14) {
		throw invalid_argument("不能除以零");
	}
	
	vector<double> result(v.size());
	for (size_t i = 0; i < v.size(); ++i) {
		result[i] = v[i] / scalar;
	}
	return result;
}

vector<double> operator-(const vector<double>& v) {
	vector<double> result(v.size());
	for (size_t i = 0; i < v.size(); ++i) {
		result[i] = -v[i];
	}
	return result;
}

// ======================= 复合赋值运算符实现 ======================= //

vector<double>& operator+=(vector<double>& a, const vector<double>& b) {
	if (a.size() != b.size()) {
		throw invalid_argument("向量大小必须相同");
	}
	
	for (size_t i = 0; i < a.size(); ++i) {
		a[i] += b[i];
	}
	return a;
}

vector<double>& operator-=(vector<double>& a, const vector<double>& b) {
	if (a.size() != b.size()) {
		throw invalid_argument("向量大小必须相同");
	}
	
	for (size_t i = 0; i < a.size(); ++i) {
		a[i] -= b[i];
	}
	return a;
}

vector<double>& operator*=(vector<double>& v, double scalar) {
	for (size_t i = 0; i < v.size(); ++i) {
		v[i] *= scalar;
	}
	return v;
}

vector<double>& operator/=(vector<double>& v, double scalar) {
	if (abs(scalar) < 1e-14) {
		throw invalid_argument("不能除以零");
	}
	
	for (size_t i = 0; i < v.size(); ++i) {
		v[i] /= scalar;
	}
	return v;
}

// ======================= 输入输出和工具函数实现 ======================= //

vector<double> readVec(const string& filename) {
	// 一行一个数据
	ifstream fin(filename);
	if (!fin) {
		throw invalid_argument("无法打开文件: " + filename);
	}
	vector<double> result;
	double value;
	while (fin >> value) {
		result.push_back(value);
	}
	fin.close();
	return result;
}

ostream& operator<<(ostream& os, const vector<double>& v) {
	if (v.empty()) {
		os << "[]";
		return os;
	}
	
	// 计算最大显示宽度
	size_t max_width = 0;
	for (double val : v) {
		ostringstream oss;
		oss << fixed << setprecision(4);  // 默认精度
		oss << val;
		string s = oss.str();
		if (s.length() > max_width) {
			max_width = s.length();
		}
	}
	
	// 输出向量（默认按行向量格式）
	os << "[";
	for (size_t i = 0; i < v.size(); ++i) {
		ostringstream oss;
		oss << fixed << setprecision(4);
		oss << v[i];
		os << setw(max_width) << oss.str();
		
		if (i < v.size() - 1) {
			os << " ";
		}
	}
	os << "]";
	
	// 恢复默认格式
	return os << defaultfloat << setprecision(6);
}

void print(const vector<double>& v, int precision, bool asColumnVector) {
	// 先告知向量长度
	cout << "向量长度: " << v.size() << endl;
	
	if (v.empty()) {
		cout << "（空向量）" << endl;
		return;
	}

	// 计算最大显示宽度
	size_t max_width = 0;
	for (double val : v) {
		ostringstream oss;
		oss << fixed << setprecision(precision);
		oss << val;
		string s = oss.str();
		if (s.length() > max_width) {
			max_width = s.length();
		}
	}

	if (asColumnVector) {
		// 按列向量打印（每个元素占一行）
		for (size_t i = 0; i < v.size(); ++i) {
			cout << fixed << setprecision(precision) 
					<< setw(max_width) << v[i] << endl;
		}
	} else {
		// 按行向量打印（所有元素在一行）
		cout << "[";
		for (size_t i = 0; i < v.size(); ++i) {
			ostringstream oss;
			oss << fixed << setprecision(precision);
			oss << v[i];
			cout << setw(max_width) << oss.str();
			
			if (i < v.size() - 1) {
				cout << " ";
			}
		}
		cout << "]" << endl;
	}

	// 恢复默认格式
	cout << defaultfloat << setprecision(6);
}

void print(const vector<vector<double>>& matrix, int precision) {
	// 检查矩阵是否为空
	if (matrix.empty()) {
		cout << "矩阵大小: ( 0 x 0 )" << endl;
		cout << "(空矩阵)" << endl;
		return;
	}
	
	size_t rows = matrix.size();
	size_t cols = matrix[0].size();
	
	// 检查矩阵的合法性：所有行的列数是否相同
	for (size_t i = 1; i < rows; ++i) {
		if (matrix[i].size() != cols) {
			cout << "错误：矩阵行大小不一致！" << endl;
			return;
		}
	}
	
	// 处理0列的情况
	if (cols == 0) {
		cout << "矩阵大小: ( " << rows << " x 0 )" << endl;
		cout << "(空矩阵)" << endl;
		return;
	}
	
	cout << "矩阵大小: ( " << rows << " x " << cols << " )" << endl;

	size_t max_width = 0;

	// 遍历所有元素，计算它们格式化后的最大显示宽度
	for (const auto& row : matrix) {
		for (double val : row) {
			ostringstream oss;
			oss << fixed << setprecision(precision);
			oss << val;
			string s = oss.str();
			if (s.length() > max_width) {
				max_width = s.length();
			}
		}
	}

	// 打印矩阵内容
	for (const auto& row : matrix) {
		for (double val : row) {
			ostringstream os;
			os << fixed << setprecision(precision);
			os << val;
			cout << setw(max_width) << os.str() << " ";
		}
		cout << endl;
	}
	
	// 恢复默认格式
	cout << defaultfloat << setprecision(6);
}

vector<double> zeros(size_t n) {
	return vector<double>(n, 0.0);
}

vector<double> ones(size_t n) {
	return vector<double>(n, 1.0);
}

vector<double> random(size_t n, double min, double max) {
	vector<double> result(n);
	random_device rd;
	mt19937 gen(rd());
	uniform_real_distribution<double> dis(min, max);
	
	for (size_t i = 0; i < n; ++i) {
		result[i] = dis(gen);
	}
	return result;
}

bool areApproxEqual(const vector<double>& a, const vector<double>& b, double tol) {
	if (a.size() != b.size()) {
		return false;
	}
	
	for (size_t i = 0; i < a.size(); ++i) {
		if (abs(a[i] - b[i]) > tol) {
			return false;
		}
	}
	return true;
}

} // namespace VectorOps