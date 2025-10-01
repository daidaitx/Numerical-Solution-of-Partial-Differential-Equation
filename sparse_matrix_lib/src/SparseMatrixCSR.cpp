#include "SparseMatrixCSR.hpp"
#include <stdexcept>

// 构造函数：构造一个 r x c 全零矩阵
SparseMatrixCSR::SparseMatrixCSR(size_t r, size_t c) : 
	rows(r), cols(c),
	values(), col_indices(),
	row_ptrs(r + 1, 0)  // 初始化行指针为全0
{
	if (r == 0 || c == 0) 
		throw std::invalid_argument("矩阵行数和列数必须均大于0。");
}

// getter函数
size_t SparseMatrixCSR::getRows() const noexcept {
	return rows;
}

size_t SparseMatrixCSR::getCols() const noexcept {
	return cols;
}

std::pair<size_t, size_t> SparseMatrixCSR::getSize() const noexcept {
	return {rows, cols};
}

size_t SparseMatrixCSR::getNNZ() const noexcept {
	return values.size();
}

// 流输出矩阵
std::ostream& operator<<(std::ostream& os, const SparseMatrixCSR& mat) {
    size_t rows = mat.rows;
    size_t cols = mat.cols;
    const auto& values = mat.values;
    const auto& col_indices = mat.col_indices;
    const auto& row_ptrs = mat.row_ptrs;

    os << "CSR稀疏矩阵 (" << rows << " x " << cols << "), NNZ = " << values.size() << "\n";
    os << "格式: (row, col) -> value\n";

    for (size_t i = 0; i < rows; ++i) {
        size_t start = row_ptrs[i];
        // size_t end = (i + 1 < rows) ? row_ptrs[i + 1] : values.size();
		size_t end = row_ptrs[i + 1];
        for (size_t j = start; j < end; ++j) {
            os << "  (" << i << ", " << col_indices[j] << ") -> " << values[j] << "\n";
		}
    }

    return os;
}