#include <iostream>
#include <vector>
#include <fstream>
#include "SparseMatrixCSR.hpp"

using namespace std;

int main() {
	cout << "当前运行的是 TestCore 测试程序。" << endl;

	cout << "========================================" << endl;
	cout << "测试1a - 创建全零 CSR 矩阵(0行)" << endl;
	const size_t testRows_1a = 0;
	const size_t testCols_1a = 4;
	try {
		SparseMatrixCSR matrix_1a(testRows_1a, testCols_1a);
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}
	
	cout << "........................................" << endl;
	cout << "测试1b - 创建全零 CSR 矩阵并流式输出" << endl;
	const size_t testRows_1b = 3;
	const size_t testCols_1b = 4;
	SparseMatrixCSR matrix_1b(testRows_1b, testCols_1b);
	matrix_1b.print();

	cout << "----------------------------------------" << endl;
	cout << "测试2 - 输出矩阵的行数、列数、大小、非零元素个数" << endl;
	const size_t testRows_2 = 3;
	const size_t testCols_2 = 4;
	SparseMatrixCSR matrix_2(testRows_2, testCols_2);
	cout << "行数：" << matrix_2.getRows() << endl;
	cout << "列数：" << matrix_2.getCols() << endl;
	cout << "大小：(" << matrix_2.getSize().first << "x" << matrix_2.getSize().second << ")" << endl;
	cout << "非零元素个数：" << matrix_2.getNNZ() << endl;

	cout << "----------------------------------------" << endl;
	cout << "测试3a - 稠密矩阵转CSR稀疏矩阵" << endl;
	const vector<vector<double>> denseMatrix_3a = {
		{1, 1, 0, 4},
		{0, 0, 0, 5},
		{0, 0, 0, 0},
		{1, 0, 4, 0}
	};
	SparseMatrixCSR csrMatrix_3a(denseMatrix_3a);
	cout << "稠密矩阵转换得到的CSR矩阵: " << endl;
	csrMatrix_3a.print();

	cout << "........................................" << endl;
	cout << "测试3b - 稠密矩阵转CSR稀疏矩阵(0行)" << endl;
	const vector<vector<double>> denseMatrix_3b = {};
	try {
		SparseMatrixCSR csrMatrix_3b(denseMatrix_3b);
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}

	cout << "........................................." << endl;
	cout << "测试3c - 稠密矩阵转CSR稀疏矩阵(0列)" << endl;
	const vector<vector<double>> denseMatrix_3c = { {}, {}, {} };
	try {
		SparseMatrixCSR csrMatrix_3c(denseMatrix_3c);
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}

	cout << ".........................................." << endl;
	cout << "测试3d - 稠密矩阵转CSR稀疏矩阵(非方阵)" << endl;
	const vector<vector<double>> denseMatrix_3d = {
		{1, 1, 0, 4},
		{0, 0, 0, 5},
		{0, 0, 0, 0},
		{1, 0, 4}
	};
	try {
		SparseMatrixCSR csrMatrix_3d(denseMatrix_3d);
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}

	cout << "----------------------------------------" << endl;
	cout << "测试4a - COO稀疏矩阵转CSR稀疏矩阵" << endl;
	const size_t rows_4a = 4;
	const size_t cols_4a = 4;
	const vector<size_t> row_indices_4a = { 0, 0, 0, 1, 3, 3 };
	const vector<size_t> col_indices_4a = { 0, 1, 3, 3, 0, 2 };
	const vector<double> values_4a = { 1, 1, 4, 5, 1, 4 };
	SparseMatrixCSR csrMatrix_4a(rows_4a, cols_4a, row_indices_4a, col_indices_4a, values_4a);
	csrMatrix_4a.print();

	csrMatrix_4a.saveToFile("./sparse_matrix_lib/test/data/matrix.coo");

	cout << "........................................." << endl;
	cout << "测试4b - COO稀疏矩阵转CSR稀疏矩阵(列优先顺序)" << endl;
	const size_t rows_4b = 4;
	const size_t cols_4b = 4;
	const vector<size_t> row_indices_4b = { 0, 3, 0, 3, 0, 1 };
	const vector<size_t> col_indices_4b = { 0, 0, 1, 2, 3, 3 };
	const vector<double> values_4b = { 1, 1, 1, 4, 4, 5 };
	SparseMatrixCSR csrMatrix_4b(rows_4b, cols_4b, row_indices_4b, col_indices_4b, values_4b);
	csrMatrix_4b.print();

	cout << "........................................." << endl;
	cout << "测试4c - COO稀疏矩阵转CSR稀疏矩阵(打乱顺序)" << endl;
	const size_t rows_4c = 4;
	const size_t cols_4c = 4;
	const vector<size_t> row_indices_4c = { 0, 0, 1, 0, 3, 3 };
	const vector<size_t> col_indices_4c = { 3, 1, 3, 0, 2, 0 };
	const vector<double> values_4c = { 4, 1, 5, 1, 4, 1 };
	SparseMatrixCSR csrMatrix_4c(rows_4c, cols_4c, row_indices_4c, col_indices_4c, values_4c);
	csrMatrix_4c.print();

	cout << "........................................." << endl;
	cout << "测试4d - COO稀疏矩阵转CSR稀疏矩阵(0行)" << endl;
	const size_t rows_4d = 0;
	const size_t cols_4d = 4;
	const vector<size_t> row_indices_4d = {};
	const vector<size_t> col_indices_4d = {};
	const vector<double> values_4d = {};
	try {
		SparseMatrixCSR csrMatrix_4d(rows_4d, cols_4d, row_indices_4d, col_indices_4d, values_4d);
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}

	cout << "........................................." << endl;
	cout << "测试4e - COO稀疏矩阵转CSR稀疏矩阵(长度不一致)" << endl;
	const size_t rows_4e = 4;
	const size_t cols_4e = 4;
	const vector<size_t> row_indices_4e = { 0, 0, 0, 1, 3, 3};
	const vector<size_t> col_indices_4e = { 0, 1, 3, 3, 0, 2};
	const vector<double> values_4e = { 1, 1, 4, 5, 1 };
	try {
		SparseMatrixCSR csrMatrix_4e(rows_4e, cols_4e, row_indices_4e, col_indices_4e, values_4e);
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}
	
	cout << "........................................." << endl;
	cout << "测试4f - COO稀疏矩阵转CSR稀疏矩阵(索引越界)" << endl;
	const size_t rows_4f = 4;
	const size_t cols_4f = 4;
	const vector<size_t> row_indices_4f = { 0, 0, 0, 1, 3, 3, 4 };
	const vector<size_t> col_indices_4f = { 0, 1, 3, 3, 0, 2, 4 };
	const vector<double> values_4f = { 1, 1, 4, 5, 1, 4, 12.456 };
	try {
		SparseMatrixCSR csrMatrix_4f(rows_4f, cols_4f, row_indices_4f, col_indices_4f, values_4f);
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}

	cout << "-----------------------------------------" << endl;
	cout << "测试5a - 从.coo文件读入CSR稀疏矩阵" << endl;
	SparseMatrixCSR csrMatrix_5a("./sparse_matrix_lib/test/data/matrix.coo");
	csrMatrix_5a.print();

	cout << "........................................." << endl;
	cout << "测试5b - 从.mtx文件读入CSR稀疏矩阵(未来实现)" << endl;
	try {
		SparseMatrixCSR csrMatrix_5b("./sparse_matrix_lib/test/data/matrix.mtx");
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}

	cout << "........................................." << endl;
	cout << "测试5c - 从.wtf文件读入CSR稀疏矩阵(不支持)" << endl;
	try {
		SparseMatrixCSR csrMatrix_5c("./sparse_matrix_lib/test/data/matrix.wtf");
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}

	cout << "........................................." << endl;
	cout << "测试5d - 从无后缀名的文件读入CSR稀疏矩阵(不支持)" << endl;
	try {
		SparseMatrixCSR csrMatrix_5d("./sparse_matrix_lib/test/data/matrix");
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}

	cout << "........................................." << endl;
	cout << "测试5e - 从有错误的.coo文件读入CSR稀疏矩阵" << endl;
	try {
		SparseMatrixCSR csrMatrix_5e("./sparse_matrix_lib/test/data/matrix_err.coo");
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}

	cout << "........................................." << endl;
	cout << "测试5f - 从不存在的文件读入CSR稀疏矩阵" << endl;
	try {
		SparseMatrixCSR csrMatrix_5f("./sparse_matrix_lib/test/data/NonExist.coo");
	}
	catch(const std::exception& e) {
		std::cerr << e.what() << '\n';
	}

	cout << "-----------------------------------------" << endl;
	cout << "测试6 - 保存到.coo文件" << endl;
	const size_t rows_6 = 4;
	const size_t cols_6 = 4;
	const vector<size_t> row_indices_6 = { 0, 0, 0, 1, 3, 3 };
	const vector<size_t> col_indices_6 = { 0, 1, 3, 3, 0, 2 };
	const vector<double> values_6 = { 1, 1, 4, 5, 1, 4 };
	SparseMatrixCSR csrMatrix_6(rows_6, cols_6, row_indices_6, col_indices_6, values_6);
	csrMatrix_6.saveToFile("./sparse_matrix_lib/test/data/matrix_OUTPUT.coo");
	// 输出文件到流以检查
	std::ifstream inputFile("./sparse_matrix_lib/test/data/matrix_OUTPUT.coo"); // 打开刚才保存的 COO 文件
	if (!inputFile.is_open()) {
		std::cerr << "无法打开文件用于读取: ./sparse_matrix_lib/test/data/matrix_OUTPUT.coo" << std::endl;
	} else {
		std::cout << "文件内容（COO 格式，每行：行 列 值）：" << std::endl;
		std::string line;
		while (std::getline(inputFile, line)) {
			std::cout << line << std::endl; // 直接输出每一行
		}
		inputFile.close();  // 关闭文件（其实析构时也会自动关，但显式关闭是好习惯）
		std::cout << "文件内容已打印到控制台。" << std::endl;
	}
	
	return 0;
}