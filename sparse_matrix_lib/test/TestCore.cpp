#include "SparseMatrixCSR.hpp"
#include "VectorOperations.hpp"
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
using namespace VectorOps;

int main() {
	cout << "当前运行的是 TestCore 测试程序。" << endl;

	cout << "=========================================" << endl;
	cout << "测试1a - 创建全零 CSR 矩阵(0行)" << endl;
	const size_t testRows_1a = 0;
	const size_t testCols_1a = 4;
	try { SparseMatrixCSR matrix_1a(testRows_1a, testCols_1a); } catch(const exception& e) { cerr << e.what() << '\n'; }
	
	cout << "........................................." << endl;
	cout << "测试1b - 创建全零 CSR 矩阵并流式输出" << endl;
	const size_t testRows_1b = 3;
	const size_t testCols_1b = 4;
	SparseMatrixCSR matrix_1b(testRows_1b, testCols_1b);
	matrix_1b.print();

	cout << "-----------------------------------------" << endl;
	cout << "测试2 - 输出矩阵的行数、列数、大小、非零元素个数" << endl;
	const size_t testRows_2 = 3;
	const size_t testCols_2 = 4;
	SparseMatrixCSR matrix_2(testRows_2, testCols_2);
	cout << "行数：" << matrix_2.getRows() << endl;
	cout << "列数：" << matrix_2.getCols() << endl;
	cout << "大小：(" << matrix_2.getSize().first << "x" << matrix_2.getSize().second << ")" << endl;
	cout << "非零元素个数：" << matrix_2.getNNZ() << endl;

	cout << "-----------------------------------------" << endl;
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

	cout << "........................................." << endl;
	cout << "测试3b - 稠密矩阵转CSR稀疏矩阵(0行)" << endl;
	const vector<vector<double>> denseMatrix_3b = {};
	try { SparseMatrixCSR csrMatrix_3b(denseMatrix_3b); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "........................................." << endl;
	cout << "测试3c - 稠密矩阵转CSR稀疏矩阵(0列)" << endl;
	const vector<vector<double>> denseMatrix_3c = { {}, {}, {} };
	try { SparseMatrixCSR csrMatrix_3c(denseMatrix_3c); } catch(const exception& e) { cerr << e.what() << '\n';	}

	cout << "........................................." << endl;
	cout << "测试3d - 稠密矩阵转CSR稀疏矩阵(非方阵)" << endl;
	const vector<vector<double>> denseMatrix_3d = {
		{1, 1, 0, 4},
		{0, 0, 0, 5},
		{0, 0, 0, 0},
		{1, 0, 4}
	};
	try { SparseMatrixCSR csrMatrix_3d(denseMatrix_3d); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "-----------------------------------------" << endl;
	cout << "测试4a - COO稀疏矩阵转CSR稀疏矩阵(1-based)" << endl;
	const size_t rows_4a = 4;
	const size_t cols_4a = 4;
	const vector<size_t> row_indices_4a = { 1, 1, 1, 2, 4, 4 };
	const vector<size_t> col_indices_4a = { 1, 2, 4, 4, 1, 3 };
	const vector<double> values_4a = { 1, 1, 4, 5, 1, 4 };
	SparseMatrixCSR csrMatrix_4a(rows_4a, cols_4a, row_indices_4a, col_indices_4a, values_4a);
	csrMatrix_4a.print();

	cout << "........................................." << endl;
	cout << "测试4b - COO稀疏矩阵转CSR稀疏矩阵(1-based)(列优先顺序)" << endl;
	const size_t rows_4b = 4;
	const size_t cols_4b = 4;
	const vector<size_t> row_indices_4b = { 1, 4, 1, 4, 1, 2 };
	const vector<size_t> col_indices_4b = { 1, 1, 2, 3, 4, 4 };
	const vector<double> values_4b = { 1, 1, 1, 4, 4, 5 };
	SparseMatrixCSR csrMatrix_4b(rows_4b, cols_4b, row_indices_4b, col_indices_4b, values_4b);
	csrMatrix_4b.print();

	cout << "........................................." << endl;
	cout << "测试4c - COO稀疏矩阵转CSR稀疏矩阵(1-based)(打乱顺序)" << endl;
	const size_t rows_4c = 4;
	const size_t cols_4c = 4;
	const vector<size_t> row_indices_4c = { 1, 1, 2, 1, 4, 4 };
	const vector<size_t> col_indices_4c = { 4, 2, 4, 1, 3, 1 };
	const vector<double> values_4c = { 4, 1, 5, 1, 4, 1 };
	SparseMatrixCSR csrMatrix_4c(rows_4c, cols_4c, row_indices_4c, col_indices_4c, values_4c);
	csrMatrix_4c.print();

	cout << "........................................." << endl;
	cout << "测试4d - COO稀疏矩阵转CSR稀疏矩阵(1-based)(0行)" << endl;
	const size_t rows_4d = 0;
	const size_t cols_4d = 4;
	const vector<size_t> row_indices_4d = {};
	const vector<size_t> col_indices_4d = {};
	const vector<double> values_4d = {};
	try { SparseMatrixCSR csrMatrix_4d(rows_4d, cols_4d, row_indices_4d, col_indices_4d, values_4d); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "........................................." << endl;
	cout << "测试4e - COO稀疏矩阵转CSR稀疏矩阵(1-based)(长度不一致)" << endl;
	const size_t rows_4e = 4;
	const size_t cols_4e = 4;
	const vector<size_t> row_indices_4e = { 1, 1, 1, 2, 4, 4};
	const vector<size_t> col_indices_4e = { 1, 2, 4, 4, 1, 3};
	const vector<double> values_4e = { 1, 1, 4, 5, 1 };
	try { SparseMatrixCSR csrMatrix_4e(rows_4e, cols_4e, row_indices_4e, col_indices_4e, values_4e); } catch(const exception& e) { cerr << e.what() << '\n'; }
	
	cout << "........................................." << endl;
	cout << "测试4f - COO稀疏矩阵转CSR稀疏矩阵(1-based)(索引越界)" << endl;
	const size_t rows_4f = 4;
	const size_t cols_4f = 4;
	const vector<size_t> row_indices_4f = { 1, 1, 1, 2, 4, 4, 5 };
	const vector<size_t> col_indices_4f = { 1, 2, 4, 4, 1, 3, 5 };
	const vector<double> values_4f = { 1, 1, 4, 5, 1, 4, 12.456 };
	try { SparseMatrixCSR csrMatrix_4f(rows_4f, cols_4f, row_indices_4f, col_indices_4f, values_4f); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "-----------------------------------------" << endl;
	cout << "测试4.5a - CSR格式数据生成CSR稀疏矩阵(0-based)" << endl;
	const size_t rows__4a = 4;
	const size_t cols__4a = 4;
	const vector<double> values__4a =      { 1, 1, 4, 5, 1, 4 };
	const vector<size_t> col_indices__4a = { 0, 1, 3, 3, 0, 2 };
	const vector<size_t> row_ptrs__4a =    { 0, 3, 4, 4, 6 };
	SparseMatrixCSR csrMatrix__4a(rows__4a, cols__4a, values__4a, col_indices__4a, row_ptrs__4a);
	csrMatrix__4a.print();
	csrMatrix__4a.printDense();

	cout << "........................................." << endl;
	cout << "测试4.5b - CSR格式数据生成CSR稀疏矩阵(0-based)（失败案例）" << endl;
	cout << "行数为0：" << endl;
	size_t rows__4b = 0;
	size_t cols__4b = 4;
	vector<double> values__4b =      { 1, 1, 4, 5, 1, 4 };
	vector<size_t> col_indices__4b = { 0, 1, 3, 3, 0, 2 };
	vector<size_t> row_ptrs__4b =    { 0, 3, 4, 4, 6 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "列数为0：" << endl;
	rows__4b = 4;
	cols__4b = 0;
	values__4b =      { 1, 1, 4, 5, 1, 4 };
	col_indices__4b = { 0, 1, 3, 3, 0, 2 };
	row_ptrs__4b =    { 0, 3, 4, 4, 6 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "列索引和值长度不一致：" << endl;
	rows__4b = 4;
	cols__4b = 4;
	values__4b =      { 1, 1, 4, 5, 1, 4 };
	col_indices__4b = { 0, 1, 3, 3, 0, 2, 4 };
	row_ptrs__4b =    { 0, 3, 4, 4, 6 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "行指针数量不等于行数+1：" << endl;
	rows__4b = 4;
	cols__4b = 4;
	values__4b =      { 1, 1, 4, 5, 1, 4 };
	col_indices__4b = { 0, 1, 3, 3, 0, 2 };
	row_ptrs__4b =    { 0, 3, 4, 4, 6, 6 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "列索引超过矩阵维度：" << endl;
	rows__4b = 4;
	cols__4b = 4;
	values__4b =      { 1, 1, 4, 5, 1, 4 };
	col_indices__4b = { 0, 1, 3, 3, 0, 4 };
	row_ptrs__4b =    { 0, 3, 4, 4, 6 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "行指针非递增：" << endl;
	rows__4b = 4;
	cols__4b = 4;
	values__4b =      { 1, 1, 4, 5, 1, 4 };
	col_indices__4b = { 0, 1, 3, 3, 0, 2 };
	row_ptrs__4b =    { 0, 3, 4, 3, 6 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "值显式包含零元素：" << endl;
	rows__4b = 4;
	cols__4b = 4;
	values__4b =      { 1, 1, 4, 5, 0, 4 };
	col_indices__4b = { 0, 1, 3, 3, 0, 2 };
	row_ptrs__4b =    { 0, 3, 4, 4, 6 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "某行列索引非递增：" << endl;
	rows__4b = 4;
	cols__4b = 4;
	values__4b =      { 1, 1, 4, 5, 1, 4 };
	col_indices__4b = { 0, 3, 1, 3, 0, 2 };
	row_ptrs__4b =    { 0, 3, 4, 4, 6 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "重复元素：" << endl;
	rows__4b = 4;
	cols__4b = 4;
	values__4b =      { 1, 1, 4, 5, 1, 4 };
	col_indices__4b = { 0, 0, 3, 3, 0, 2 };
	row_ptrs__4b =    { 0, 3, 4, 4, 6 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "行指针的第一个元素不为0：" << endl;
	rows__4b = 4;
	cols__4b = 4;
	values__4b =      { 1, 1, 4, 5, 1, 4 };
	col_indices__4b = { 0, 1, 3, 3, 0, 2 };
	row_ptrs__4b =    { 3, 3, 4, 4, 6 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "行指针的最后一个元素不等于非零元素个数：" << endl;
	rows__4b = 4;
	cols__4b = 4;
	values__4b =      { 1, 1, 4, 5, 1, 4 };
	col_indices__4b = { 0, 1, 3, 3, 0, 2 };
	row_ptrs__4b =    { 0, 3, 4, 4, 7 };
	try { SparseMatrixCSR csrMatrix__4b(rows__4b, cols__4b, values__4b, col_indices__4b, row_ptrs__4b); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "-----------------------------------------" << endl;
	cout << "测试5a - 从.coo文件读入CSR稀疏矩阵(1-based)" << endl;
	SparseMatrixCSR csrMatrix_5a("./sparse_matrix_lib/test/data/matrix.coo");
	csrMatrix_5a.print();

	cout << "........................................." << endl;
	cout << "测试5b - 从.mtx文件读入CSR稀疏矩阵(未来实现)" << endl;
	try { SparseMatrixCSR csrMatrix_5b("./sparse_matrix_lib/test/data/matrix.mtx"); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "........................................." << endl;
	cout << "测试5c - 从.wtf文件读入CSR稀疏矩阵(不支持)" << endl;
	try { SparseMatrixCSR csrMatrix_5c("./sparse_matrix_lib/test/data/matrix.wtf"); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "........................................." << endl;
	cout << "测试5d - 从无后缀名的文件读入CSR稀疏矩阵(不支持)" << endl;
	try { SparseMatrixCSR csrMatrix_5d("./sparse_matrix_lib/test/data/matrix"); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "........................................." << endl;
	cout << "测试5e - 从有错误的.coo文件读入CSR稀疏矩阵" << endl;
	try { SparseMatrixCSR csrMatrix_5e("./sparse_matrix_lib/test/data/matrix_err.coo"); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "........................................." << endl;
	cout << "测试5f - 从不存在的文件读入CSR稀疏矩阵" << endl;
	try { SparseMatrixCSR csrMatrix_5f("./sparse_matrix_lib/test/data/NonExist.coo"); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "-----------------------------------------" << endl;
	cout << "测试6 - 保存到.coo文件(1-based)" << endl;
	const size_t rows_6 = 4;
	const size_t cols_6 = 4;
	const vector<size_t> row_indices_6 = { 1, 1, 1, 2, 4, 4 };
	const vector<size_t> col_indices_6 = { 1, 2, 4, 4, 1, 3 };
	const vector<double> values_6 = { 1, 1, 4, 5, 1, 4 };
	SparseMatrixCSR csrMatrix_6(rows_6, cols_6, row_indices_6, col_indices_6, values_6);
	csrMatrix_6.saveToFile("./sparse_matrix_lib/test/data/matrix_OUTPUT.coo");
	// 输出文件到流以检查
	ifstream inputFile("./sparse_matrix_lib/test/data/matrix_OUTPUT.coo"); // 打开刚才保存的 COO 文件
	if (!inputFile.is_open()) {
		cerr << "无法打开文件用于读取: ./sparse_matrix_lib/test/data/matrix_OUTPUT.coo" << endl;
	} else {
		cout << "文件内容（COO 格式，行 列 值，1-based 索引）：" << endl;
		string line;
		while (getline(inputFile, line)) {
			cout << line << endl; // 直接输出每一行
		}
		inputFile.close();  // 关闭文件（其实析构时也会自动关，但显式关闭是好习惯）
		cout << "文件内容已打印。" << endl;
	}
	
	cout << "-----------------------------------------" << endl;
	cout << "测试7a - 插入单元素到CSR稀疏矩阵(1-based)" << endl;
	const size_t rows_7a = 4;
	const size_t cols_7a = 4;
	SparseMatrixCSR csrMatrix_7a(rows_7a, cols_7a);
	cout << "初始矩阵：" << endl;
	csrMatrix_7a.print();
	cout << "新增元素(1,1,1.22), (1,2,23.4)：" << endl;
	csrMatrix_7a.setValue(1, 1, 1.22);
	csrMatrix_7a.setValue(1, 2, 23.4);
	csrMatrix_7a.print();
	cout << "更新元素(1,1,1)：" << endl;
	csrMatrix_7a.setValue(1, 1, 1);
	csrMatrix_7a.print();
	cout << "新增元素(3,3,-1)：" << endl;
	csrMatrix_7a.setValue(3, 3, -1);
	csrMatrix_7a.print();
	cout << "删除元素(1,2,0)：" << endl;
	csrMatrix_7a.setValue(1, 2, 0);
	csrMatrix_7a.print();

	cout << "........................................." << endl;
	cout << "测试7b - 其他方法插入到CSR稀疏矩阵(1-based)" << endl;
	const size_t rows_7b = 4;
	const size_t cols_7b = 4;
	SparseMatrixCSR csrMatrix_7b(rows_7b, cols_7b);
	cout << "初始矩阵：" << endl;
	csrMatrix_7b.print();
	cout << "设置元素tuple{(1,1,1)}：" << endl;
	tuple<size_t, size_t, double> t1(1, 1, 1);
	csrMatrix_7b.setValue(t1);
	csrMatrix_7b.print();
	cout << "设置元素(vector{3,2,1,2,4}, vector{4,3,3,2,1}, vector{-1,-2,-3,-4,-5})：" << endl;
	vector<size_t> row_indices_7b = { 3, 2, 1, 2, 4 };
	vector<size_t> col_indices_7b = { 4, 3, 3, 2, 1 };
	vector<double> values_7b = { -1, -2, -3, -4, -5 };
	csrMatrix_7b.setValue(row_indices_7b, col_indices_7b, values_7b);
	csrMatrix_7b.print();
	cout << "新增元素vector{tuple{1,1,0}, tuple{2,3,0}, tuple{2,2,0}}：" << endl;
	vector<tuple<size_t, size_t, double>> tuples_7b = { make_tuple(1, 1, 0), make_tuple(2, 3, 0), make_tuple(2, 2, 0) };
	csrMatrix_7b.setValue(tuples_7b);
	csrMatrix_7b.print();

	cout << "........................................." << endl;
	cout << "测试7c - 错误插入到CSR稀疏矩阵(1-based)" << endl;
	const size_t rows_7c = 4;
	const size_t cols_7c = 4;
	SparseMatrixCSR csrMatrix_7c(rows_7c, cols_7c);
	cout << "初始矩阵：" << endl;
	csrMatrix_7c.print();
	cout << "设置错误元素组({1,1,1}, {3,2,1}, {3,3})：" << endl;
	vector<size_t> row_indices_7c = { 1, 1, 1 };
	vector<size_t> col_indices_7c = { 3, 2, 1 };
	vector<double> values_7c = { 3, 3 };
	try { csrMatrix_7c.setValue(row_indices_7c, col_indices_7c, values_7c); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "设置元素(5,1,1)：" << endl;
	try { csrMatrix_7c.setValue(5, 1, 1); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "设置元素(1,5,1)：" << endl;
	try { csrMatrix_7c.setValue(1, 5, 1); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "设置元素(0,1,3)：" << endl;
	try { csrMatrix_7c.setValue(0, 1, 3); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "设置元素(1,0,2)：" << endl;
	try { csrMatrix_7c.setValue(1, 0, 2); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "设置元素{(3,2,1), (3,2,-1)}：" << endl;
	tuple<size_t, size_t, double> t2(3, 2, 1);
	tuple<size_t, size_t, double> t3(3, 2, -1);
	vector<tuple<size_t, size_t, double>> tuples_7c = { t2, t3 };
	try { csrMatrix_7c.setValue(tuples_7c); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "当前矩阵：" << endl;
	csrMatrix_7c.print();

	cout << "-----------------------------------------" << endl;
	cout << "测试8 - 测试矩阵运算符重载功能" << endl;
	const size_t rows_8 = 4;
	const size_t cols_8 = 4;
	vector<vector<double>> values_8A = { {1, 2, 3, 4}, {5, 6, 7, 8}, {0, 0, 0, 0}, {0, 0, 0, 0} };
	vector<vector<double>> values_8B = { {1, 2, 3, 4}, {0, 0, 0, 0}, {0, 0, 0, 0}, {13, 14, 15, 16} };
	SparseMatrixCSR csrMatrix_8A(values_8A);
	SparseMatrixCSR csrMatrix_8B(values_8B);
	SparseMatrixCSR csrMatrix_8temp;
	vector<double> vec_8temp;
	vector<vector<double>> denseMatrix_8temp;
	cout << "矩阵A：" << endl;
	csrMatrix_8A.print();
	cout << "矩阵B：" << endl;
	csrMatrix_8B.print();
	cout << "矩阵括号索引(1-based)" << endl;
	cout << "矩阵A(1,1)：" << csrMatrix_8A(1, 1) << endl;
	cout << "矩阵A(2,3)：" << csrMatrix_8A(2, 3) << endl;
	cout << "矩阵A(3,2)：" << csrMatrix_8A(3, 2) << endl;
	cout << "矩阵A(4,4)：" << csrMatrix_8A(4, 4) << endl;
	try { cout << "矩阵A(5,1)：" << csrMatrix_8A(5, 1) << endl; } catch(const exception& e) { cerr << e.what() << '\n'; }
	try { cout << "矩阵A(1,5)：" << csrMatrix_8A(1, 5) << endl; } catch(const exception& e) { cerr << e.what() << '\n'; }
	try { cout << "矩阵A(0,1)：" << csrMatrix_8A(0, 1) << endl; } catch(const exception& e) { cerr << e.what() << '\n'; }
	try { cout << "矩阵A(1,0)：" << csrMatrix_8A(1, 0) << endl; } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "-----------------------------------------" << endl;
	cout << "测试9 - CSR稀疏矩阵转稠密矩阵以及稠密打印" << endl;
	vector<vector<double>> values_9 = { {3.1, -4, 1.5, -9.2, 0.3}, {-6.5, -3.5, 89.7, 9.32, 0}, {-38.4, 62.6, -433.8, 3.279, 0}, {5, 0, -2.88, 419716, -9} };
	SparseMatrixCSR csrMatrix_9(values_9);
	cout << "CSR稀疏矩阵稀疏打印：" << endl;
	csrMatrix_9.print();
	cout << "CSR稀疏矩阵转稠密矩阵并显式打印：" << endl;
	vector<vector<double>> denseMatrix_9 = csrMatrix_9.toDense();
	print(denseMatrix_9);
	cout << "CSR稀疏矩阵稠密打印：" << endl;
	csrMatrix_9.printDense();
	cout << "CSR稀疏矩阵稠密打印（presicion=2）：" << endl;
	csrMatrix_9.printDense(2);
	cout << "CSR稀疏矩阵稠密打印（presicion=8）：" << endl;
	csrMatrix_9.printDense(8);

	cout << "-----------------------------------------" << endl;
	cout << "测试10 - CSR稀疏矩阵转置" << endl;
	vector<vector<double>> values_10 = { {3.1, -4, 1.5, -9.2, 0.3}, {-6.5, -3.5, 89.7, 9.32, 0}, {-38.4, 62.6, -433.8, 3.279, 0}, {5, 0, -2.88, 419716, -9} };
	SparseMatrixCSR csrMatrix_10(values_10);
	cout << "CSR稀疏矩阵：" << endl;
	csrMatrix_10.printDense();
	cout << "CSR稀疏矩阵的转置：" << endl;
	SparseMatrixCSR csrMatrix_10T = csrMatrix_10.transpose();
	csrMatrix_10T.printDense();
	cout << "CSR稀疏矩阵的转置（使用简写）：" << endl;
	SparseMatrixCSR csrMatrix_10t = csrMatrix_10.t();
	csrMatrix_10T.printDense();

	cout << "-----------------------------------------" << endl;
	cout << "测试11 - CSR稀疏矩阵范数" << endl;
	vector<vector<double>> values_11 = { {3.1, -4, 1.5, -9.2, 0.3}, {-6.5, -3.5, 89.7, 9.32, 0}, {-38.4, 62.6, -433.8, 3.279, 0}, {5, 0, -2.88, 419716, -9} };
	SparseMatrixCSR csrMatrix_11(values_11);
	cout << "CSR稀疏矩阵：" << endl;
	csrMatrix_11.printDense();
	cout << "CSR稀疏矩阵的无穷范数：" << endl;
	double norm_inf_11 = csrMatrix_11.norm("inf");
	cout << norm_inf_11 << endl;
	cout << "CSR稀疏矩阵的1范数：" << endl;
	double norm_1_11 = csrMatrix_11.norm("1");
	cout << norm_1_11 << endl;
	cout << "CSR稀疏矩阵的frobenius范数：" << endl;
	double norm_fro_11 = csrMatrix_11.norm("fro");
	cout << norm_fro_11 << endl;
	cout << "CSR稀疏矩阵的按元素max范数：" << endl;
	double norm_max_11 = csrMatrix_11.norm("max");
	cout << norm_max_11 << endl;
	cout << "CSR稀疏矩阵的范数（填写错误）：" << endl;
	try { double norm_none_11 = csrMatrix_11.norm("这是一个不支持的范数类型"); } catch(const exception& e) { cerr << e.what() << '\n'; }
	
	cout << "-----------------------------------------" << endl;
	cout << "测试12 - 矩阵判断" << endl;
	SparseMatrixCSR csrMatrix_12a(3, 3);
	csrMatrix_12a.printDense(0);
	cout << "全零方阵：isSquare() = " << csrMatrix_12a.isSquare() << ", isSymmetric() = " << csrMatrix_12a.isSymmetric() << ", isUpperTriangular() = " << csrMatrix_12a.isUpperTriangular() << ", isLowerTriangular() = " << csrMatrix_12a.isLowerTriangular() << ", isDiagonal() = " << csrMatrix_12a.isDiagonal() << endl;
	SparseMatrixCSR csrMatrix_12b(4, 3);
	csrMatrix_12b.printDense(0);
	cout << "全零非方阵：isSquare() = " << csrMatrix_12b.isSquare() << ", isSymmetric() = " << csrMatrix_12b.isSymmetric() << ", isUpperTriangular() = " << csrMatrix_12b.isUpperTriangular() << ", isLowerTriangular() = " << csrMatrix_12b.isLowerTriangular() << ", isDiagonal() = " << csrMatrix_12b.isDiagonal() << endl;
	vector<vector<double>> values_12c = { {1, 2, 3, 4}, {5, 6, 7, 8}, {0, 0, 0, 0}, {0, 0, 0, 0} };
	SparseMatrixCSR csrMatrix_12c(values_12c);
	csrMatrix_12c.printDense(0);
	cout << "一般非对称方阵：isSquare() = " << csrMatrix_12c.isSquare() << ", isSymmetric() = " << csrMatrix_12c.isSymmetric() << ", isUpperTriangular() = " << csrMatrix_12c.isUpperTriangular() << ", isLowerTriangular() = " << csrMatrix_12c.isLowerTriangular() << ", isDiagonal() = " << csrMatrix_12c.isDiagonal() << endl;
	vector<vector<double>> values_12d = { {1, 2, 3, 4}, {2, 1, 0, 0}, {3, 0, 1, 0}, {4, 0, 0, 1}, {-9, -9, -8, -7}};
	SparseMatrixCSR csrMatrix_12d(values_12d);
	csrMatrix_12d.printDense(0);
	cout << "一般非方阵：isSquare() = " << csrMatrix_12d.isSquare() << ", isSymmetric() = " << csrMatrix_12d.isSymmetric() << ", isUpperTriangular() = " << csrMatrix_12d.isUpperTriangular() << ", isLowerTriangular() = " << csrMatrix_12d.isLowerTriangular() << ", isDiagonal() = " << csrMatrix_12d.isDiagonal() << endl;
	vector<vector<double>> values_12e = { {1, 2, 3, 4}, {2, 1, 0, 0}, {3, 0, 1, 0}, {4, 0, 0, 1} };
	SparseMatrixCSR csrMatrix_12e(values_12e);
	csrMatrix_12e.printDense(0);
	cout << "对称方阵：isSquare() = " << csrMatrix_12e.isSquare() << ", isSymmetric() = " << csrMatrix_12e.isSymmetric() << ", isUpperTriangular() = " << csrMatrix_12e.isUpperTriangular() << ", isLowerTriangular() = " << csrMatrix_12e.isLowerTriangular() << ", isDiagonal() = " << csrMatrix_12e.isDiagonal() << endl;
	vector<vector<double>> values_12f = { {1, 2, 3, 4}, {0, 2, 3, 4}, {0, 0, 3, 4}, {0, 0, 0, 4} };
	SparseMatrixCSR csrMatrix_12f(values_12f);
	csrMatrix_12f.printDense(0);
	cout << "一般上三角矩阵：isSquare() = " << csrMatrix_12f.isSquare() << ", isSymmetric() = " << csrMatrix_12f.isSymmetric() << ", isUpperTriangular() = " << csrMatrix_12f.isUpperTriangular() << ", isLowerTriangular() = " << csrMatrix_12f.isLowerTriangular() << ", isDiagonal() = " << csrMatrix_12f.isDiagonal() << endl;
	vector<vector<double>> values_12g = { {1, 0, 0, 0}, {2, 1, 0, 0}, {3, 2, 1, 0}, {4, 3, 2, 1} };
	SparseMatrixCSR csrMatrix_12g(values_12g);
	csrMatrix_12g.printDense(0);
	cout << "一般下三角矩阵：isSquare() = " << csrMatrix_12g.isSquare() << ", isSymmetric() = " << csrMatrix_12g.isSymmetric() << ", isUpperTriangular() = " << csrMatrix_12g.isUpperTriangular() << ", isLowerTriangular() = " << csrMatrix_12g.isLowerTriangular() << ", isDiagonal() = " << csrMatrix_12g.isDiagonal() << endl;
	vector<vector<double>> values_12h = { {1, 0, 0, 0}, {0, 2, 0, 0}, {0, 0, 3, 0}, {0, 0, 0, 4} };
	SparseMatrixCSR csrMatrix_12h(values_12h);
	csrMatrix_12h.printDense(0);
	cout << "对角矩阵：isSquare() = " << csrMatrix_12h.isSquare() << ", isSymmetric() = " << csrMatrix_12h.isSymmetric() << ", isUpperTriangular() = " << csrMatrix_12h.isUpperTriangular() << ", isLowerTriangular() = " << csrMatrix_12h.isLowerTriangular() << ", isDiagonal() = " << csrMatrix_12h.isDiagonal() << endl;
	vector<vector<double>> values_12i = { {1, 2, 3, 4+1e-5}, {2, 2, 0, 0}, {3, 0, 3, 0}, {4, 0, 0, 4}};
	SparseMatrixCSR csrMatrix_12i(values_12i);
	csrMatrix_12i.printDense(5);
	cout << "近似对称矩阵：isSquare() = " << csrMatrix_12i.isSquare() << ", isSymmetric() = " << csrMatrix_12i.isSymmetric() << ", isSymmetric(1e-5) = " << csrMatrix_12i.isSymmetric(1e-5) << endl;

	cout << "-----------------------------------------" << endl;
	cout << "测试13 - 矩阵获取" << endl;
	vector<vector<double>> values_13a = { {3.1, -4, 1.5, -9.2, 0.3}, {-6.5, -3.5, 89.7, 9.32, 0}, {-38.4, 62.6, -433.8, 3.279, 0}, {5, 0, -2.88, 419716, -9} };
	SparseMatrixCSR csrMatrix_13a(values_13a);
	cout << "CSR稀疏矩阵：" << endl;
	csrMatrix_13a.printDense();
	cout << "获取矩阵的对角元，生成为向量：" << endl;
	vector<double> DiagVec_13a = csrMatrix_13a.getDiagonalVector();
	print(DiagVec_13a);
	cout << endl;
	cout << "获取矩阵的对角元，生成为矩阵：" << endl;
	SparseMatrixCSR DiagMat_13a = csrMatrix_13a.getDiagonalMatrix();
	DiagMat_13a.print();
	vector<vector<double>> values_13A = { {3.1, -4, 1.5, -9.2}, {-6.5, -3.5, 89.7, 9.32}, {-38.4, 62.6, -433.8, 3.279}, {5, 0, -2.88, 419716} };
	SparseMatrixCSR csrMatrix_13A(values_13A);
	cout << "获取矩阵的上三角部分：" << endl;
	SparseMatrixCSR UpperMat_13A = csrMatrix_13A.getUpperTriangularMatrix();
	UpperMat_13A.printDense();
	cout << "获取矩阵的上三角部分（非方阵）" << endl;
	try { SparseMatrixCSR UpperMat_13a = csrMatrix_13a.getUpperTriangularMatrix(); } catch(const exception& e) { cerr << e.what() << '\n'; }
	cout << "获取矩阵的下三角部分：" << endl;
	SparseMatrixCSR LowerMat_13A = csrMatrix_13A.getLowerTriangularMatrix();
	LowerMat_13A.printDense();
	cout << "获取矩阵的下三角部分（非方阵）" << endl;
	try { SparseMatrixCSR LowerMat_13a = csrMatrix_13a.getLowerTriangularMatrix(); } catch(const exception& e) { cerr << e.what() << '\n'; }

	cout << "=========================================" << endl;
	cout << "TestCore 测试结束。" << endl;
	return 0;
}