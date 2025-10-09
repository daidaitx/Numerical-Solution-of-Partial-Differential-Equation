#include "SparseMatrixCSR.hpp"
#include "SparseMatrixSolvers.hpp"
#include "VectorOperations.hpp"
#include <iostream>
#include <vector>
#include <fstream>

using namespace std;
using namespace VectorOps;

int main() {
	cout << "当前运行的是 TestSpMV 测试程序。" << endl;

	cout << "=========================================" << endl;
	cout << "测试1 - 判断矩阵严格相等、严格不等、近似相等、近似不等" << endl;
	vector<vector<double>> matrix_1a = {{1, 0, 3}, {0, 5, 0}, {0, 8, 9}};
	vector<vector<double>> matrix_1b = {{1, 0, 3}, {0, 5, 0}, {0, 8, -9}};
	vector<vector<double>> matrix_1c = {{1, 0, 3}, {0, 5, 0}, {0, 8, 9+1e-5}};
	SparseMatrixCSR csrMatrix_1a(matrix_1a);
	SparseMatrixCSR csrMatrix_1b(matrix_1b);
	SparseMatrixCSR csrMatrix_1c(matrix_1c);
	cout << "矩阵1a：" << endl;
	csrMatrix_1a.printDense(5);
	cout << "矩阵1b：" << endl;
	csrMatrix_1b.printDense(5);
	cout << "矩阵1c：" << endl;
	csrMatrix_1c.printDense(5);
	cout << "矩阵1a == 矩阵1a：" << (csrMatrix_1a == csrMatrix_1a) << endl;
	cout << "矩阵1a == 矩阵1b：" << (csrMatrix_1a == csrMatrix_1b) << endl;
	cout << "矩阵1a == 矩阵1c：" << (csrMatrix_1a == csrMatrix_1c) << endl;
	cout << "矩阵1a != 矩阵1a：" << (csrMatrix_1a != csrMatrix_1a) << endl;
	cout << "矩阵1a != 矩阵1b：" << (csrMatrix_1a != csrMatrix_1b) << endl;
	cout << "矩阵1a != 矩阵1c：" << (csrMatrix_1a != csrMatrix_1c) << endl;
	cout << "矩阵1a ≈= 矩阵1a：" << (csrMatrix_1a.isApproxEqualto(csrMatrix_1a, 1e-5)) << endl;
	cout << "矩阵1a ≈= 矩阵1b：" << (csrMatrix_1a.isApproxEqualto(csrMatrix_1b, 1e-5)) << endl;
	cout << "矩阵1a ≈= 矩阵1c：" << (csrMatrix_1a.isApproxEqualto(csrMatrix_1c, 1e-5)) << endl;
	cout << "矩阵1a !≈ 矩阵1a：" << (csrMatrix_1a.isFarAwayFrom(csrMatrix_1a, 1e-5)) << endl;
	cout << "矩阵1a !≈ 矩阵1b：" << (csrMatrix_1a.isFarAwayFrom(csrMatrix_1b, 1e-5)) << endl;
	cout << "矩阵1a !≈ 矩阵1c：" << (csrMatrix_1a.isFarAwayFrom(csrMatrix_1c, 1e-5)) << endl;

	cout << "-----------------------------------------" << endl;
	cout << "测试2 - 矩阵累加、相加、累减、相减、单目负号" << endl;
	vector<vector<double>> matrix_2a = {{1, 0, 3}, {0, 5, 0}, {0, 8, 9}};
	vector<vector<double>> matrix_2b = {{1, 4, -3}, {1, 0, 0}, {0, -9, -9}};
	vector<vector<double>> matrix_2c = {{-1, 0, 1}, {0, 0, 0}, {-1, 0, 0}};
	SparseMatrixCSR csrMatrix_2a(matrix_2a);
	SparseMatrixCSR csrMatrix_2b(matrix_2b);
	SparseMatrixCSR csrMatrix_2c(matrix_2c);
	SparseMatrixCSR csrMatrix_2d(3, 4);
	cout << "矩阵2a：" << endl;
	csrMatrix_2a.printDense(0);
	cout << "矩阵2b：" << endl;
	csrMatrix_2b.printDense(0);
	cout << "矩阵2c：" << endl;
	csrMatrix_2c.printDense(0);
	cout << "矩阵2d：" << endl;
	csrMatrix_2d.printDense(0);
	cout << "矩阵2a + 矩阵2b：" << endl;
	(csrMatrix_2a + csrMatrix_2b).printDense(0);
	(csrMatrix_2a + csrMatrix_2b).print();
	cout << "矩阵2a - 矩阵2b：" << endl;
	(csrMatrix_2a - csrMatrix_2b).printDense(0);
	(csrMatrix_2a - csrMatrix_2b).print();
	cout << "-矩阵2a：" << endl;
	(-csrMatrix_2a).printDense(0);
	cout << "矩阵2c += 矩阵2a：" << endl;
	csrMatrix_2c += csrMatrix_2a;
	csrMatrix_2c.printDense(0);
	csrMatrix_2c.print();
	cout << "矩阵2c -= 矩阵2b：" << endl;
	csrMatrix_2c -= csrMatrix_2b;
	csrMatrix_2c.printDense(0);
	csrMatrix_2c.print();
	cout << "错误维度：矩阵2a += 矩阵2d：" << endl;
	try { csrMatrix_2a += csrMatrix_2d;	} catch (const exception& e) { cout << e.what() << endl; }
	cout << "错误维度：矩阵2a -= 矩阵2d：" << endl;
	try { csrMatrix_2a -= csrMatrix_2d;	} catch (const exception& e) { cout << e.what() << endl; }
	cout << "错误维度：矩阵2a + 矩阵2d：" << endl;
	try { (csrMatrix_2a + csrMatrix_2d).printDense(0); } catch (const exception& e) { cout << e.what() << endl; }
	cout << "错误维度：矩阵2a - 矩阵2d：" << endl;
	try { (csrMatrix_2a - csrMatrix_2d).printDense(0); } catch (const exception& e) { cout << e.what() << endl; }

	cout << "-----------------------------------------" << endl;
	cout << "测试3 - 矩阵数乘、数除" << endl;
	vector<vector<double>> matrix_3 = {{1, 0, 3}, {0, 5, 0}, {0, 8, 9}};
	SparseMatrixCSR csrMatrix_3(matrix_3);
	cout << "矩阵3：" << endl;
	csrMatrix_3.printDense(0);
	csrMatrix_3.print();
	cout << "矩阵3 * 2.0：" << endl;
	(csrMatrix_3 * 2.0).printDense();
	(csrMatrix_3 * 2.0).print();
	cout << "2.0 * 矩阵3：" << endl;
	(2.0 * csrMatrix_3).printDense();
	(2.0 * csrMatrix_3).print();
	cout << "矩阵3 * 0.0：" << endl;
	(csrMatrix_3 * 0.0).printDense();
	(csrMatrix_3 * 0.0).print();
	cout << "0.0 * 矩阵3：" << endl;
	(0.0 * csrMatrix_3).printDense();
	(0.0 * csrMatrix_3).print();
	cout << "矩阵3 / 2.0：" << endl;
	(csrMatrix_3 / 2.0).printDense();
	(csrMatrix_3 / 2.0).print();
	cout << "矩阵3 / 0.0：" << endl;
	try { (csrMatrix_3 / 0.0).printDense(); } catch (const exception& e) { cout << e.what() << endl; }

	cout << "-----------------------------------------" << endl;
	cout << "测试4 - CSR格式稀疏矩阵乘以向量、稠密矩阵" << endl;
	vector<vector<double>> matrix_4 = {{1, 0, 3}, {0, 5, 0}, {0, 8, 9}};
	vector<double> vector_4 = {1, 2, 3};
	vector<double> wrong_vector_4 = {1, 2, 3, 4};
	SparseMatrixCSR csrMatrix_4(matrix_4);
	vector<vector<double>> denseMatrix_4 = {{1, 0, 0}, {2, 1, 0}, {3, -1, 1}};
	vector<vector<double>> wrong1_denseMatrix_4 = {{1, 0, 0}, {2, 1, 0}, {3, -1, 1}, {0, 0, 0}};
	vector<vector<double>> wrong2_denseMatrix_4 = {{1, 0, 0}, {2, 1, 0}, {3, -1}};
	cout << "矩阵4：" << endl;
	csrMatrix_4.printDense(0);
	cout << "向量4：" << endl;
	print(vector_4, 0);
	cout << "错误向量4：" << endl;
	print(wrong_vector_4, 0);
	cout << "稠密矩阵4：" << endl;
	print(denseMatrix_4, 0);
	cout << "错误1稠密矩阵4：" << endl;
	print(wrong1_denseMatrix_4, 0);
	cout << "错误2稠密矩阵4：" << endl;
	try { print(wrong2_denseMatrix_4, 0); } catch (const exception& e) { cout << e.what() << endl; }
	cout << "矩阵4 * 向量4：" << endl;
	vector<double> result_4_1 = csrMatrix_4 * vector_4;
	print(result_4_1);
	cout << "矩阵4 * 稠密矩阵4：" << endl;
	vector<vector<double>> result_4_2 = csrMatrix_4 * denseMatrix_4;
	print(result_4_2, 4);
	cout << "矩阵4 * 错误向量4：" << endl;
	try { csrMatrix_4 * wrong_vector_4; } catch (const exception& e) { cout << e.what() << endl; }
	cout << "矩阵4 * 错误1稠密矩阵4：" << endl;
	try { csrMatrix_4 * wrong1_denseMatrix_4; } catch (const exception& e) { cout << e.what() << endl; }
	cout << "矩阵4 * 错误2稠密矩阵4：" << endl;
	try { csrMatrix_4 * wrong2_denseMatrix_4; } catch (const exception& e) { cout << e.what() << endl; }

	cout << "-----------------------------------------" << endl;
	cout << "测试5 - CSR格式稀疏矩阵GMRES法求解线性方程组" << endl;
	vector<vector<double>> matrix_5 = {{1, 0, 3}, {0, 5, 0}, {0, 8, 9}};
	vector<double> vector_5 = {1, 2, 3};
	vector<double> wrong_vector_5 = {1, 2, 3, 4};
	SparseMatrixCSR csrMatrix_5(matrix_5);
	cout << "矩阵5：" << endl;
	csrMatrix_5.printDense(0);
	cout << "向量5：" << endl;
	print(vector_5);
	vector<double> x = csrMatrix_5.solve(vector_5);
	cout << "解向量：" << endl;
	print(x);
	// 验证残差
	vector<double> b = csrMatrix_5 * x;
	vector<double> r = vector_5 - b;
	cout << "残差：" << norm(r) << endl;
	cout << "不匹配的大小：" << endl;
	try { csrMatrix_5.solve(wrong_vector_5); } catch (const exception& e) { cout << e.what() << endl; }

	cout << "=========================================" << endl;
	cout << "TestSpMV 测试结束。" << endl;
	return 0;
}