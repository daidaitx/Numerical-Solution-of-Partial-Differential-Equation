#include <iostream>
#include "SparseMatrixCSR.hpp"

using namespace std;

int main() {
	cout << "当前运行的是 TestCore 测试程序。" << endl;

	cout << "测试1 - 创建全零 CSR 矩阵并流式输出" << endl;
	const size_t testRows = 3;
	const size_t testCols = 4;
	SparseMatrixCSR matrix(testRows, testCols);
	matrix.print();

	cout << "测试2 - 输出矩阵的行数、列数、大小、非零元素个数" << endl;
	cout << "行数：" << matrix.getRows() << endl;
	cout << "列数：" << matrix.getCols() << endl;
	cout << "大小：(" << matrix.getSize().first << "x" << matrix.getSize().second << ")" << endl;
	cout << "非零元素个数：" << matrix.getNNZ() << endl;

	
	return 0;
}