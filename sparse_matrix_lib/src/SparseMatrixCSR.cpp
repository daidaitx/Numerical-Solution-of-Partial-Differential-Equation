#include "SparseMatrixCSR.hpp"
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cctype>
#include <stdexcept>

// 构造函数：构造一个 r x c 全零矩阵
SparseMatrixCSR::SparseMatrixCSR(size_t r, size_t c) : 
	rows(r), cols(c),
	values(), col_indices(),
	row_ptrs(r + 1, 0)  // 初始化行指针为全0
{
	if (r == 0 || c == 0) {
		throw std::invalid_argument("矩阵行数和列数必须均大于0。");
	}
}

// 构造函数：从稠密矩阵（二维 vector）构造 CSR 稀疏矩阵
SparseMatrixCSR::SparseMatrixCSR(const std::vector<std::vector<double>>& dense) {
	// 1. 检查矩阵是否为空
	if (dense.empty() || dense[0].empty()) {
		throw std::invalid_argument("矩阵行数和列数必须均大于0。");
	}

	// 2. 获取行数和列数
	rows = dense.size();
	cols = dense[0].size();

	// 3. 检查是否是“矩形”矩阵（所有行的列数相同）
	for (const auto& row : dense) {
		if (row.size() != cols) {
			throw std::invalid_argument("输入的稠密矩阵必须是矩形矩阵，每一行的列数必须相同。");
		}
	}

	// 4. 初始化 CSR 数据结构
	values.clear();
	col_indices.clear();
	row_ptrs.clear();

	row_ptrs.resize(rows + 1, 0);  // row_ptrs[i] 表示第 i 行在 values 中的起始位置

	// 5. 遍历每一行、每一列，收集非零元素
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			double val = dense[i][j];
			if (val != 0.0) {
				values.push_back(val);
				col_indices.push_back(j);
			}
		}
		row_ptrs[i + 1] = values.size();  // 下一行的起始位置 = 当前 values 的大小
	}
}

// 构造函数：从 COO 格式的矩阵构造 CSR 稀疏矩阵
SparseMatrixCSR::SparseMatrixCSR(
    const size_t rows, const size_t cols,
    const std::vector<size_t>& coo_row_indices,
    const std::vector<size_t>& coo_col_indices,
    const std::vector<double>& coo_values)
    : rows(rows), cols(cols)
{
    // ===========================
    // 1. 输入合法性检查
    // ===========================
    if (rows == 0 || cols == 0) {
        throw std::invalid_argument("错误：矩阵的行数和列数必须均大于 0。");
    }

    const size_t nnz = coo_values.size();
    if (coo_row_indices.size() != nnz || coo_col_indices.size() != nnz) {
        throw std::invalid_argument("错误：COO 格式错误，行、列、值的数量必须一致。");
    }

    for (size_t i = 0; i < nnz; ++i) {
        if (coo_row_indices[i] >= rows || coo_col_indices[i] >= cols) {
            throw std::out_of_range("错误：COO 格式错误，行或列索引超出矩阵维度。");
        }
    }

    // ===========================
    // 2. 构造 COO 三元组 (row, col, value) 并排序（先行后列）
    // ===========================
    std::vector<std::tuple<size_t, size_t, double>> triples(nnz);
    for (size_t i = 0; i < nnz; ++i) {
        triples[i] = std::make_tuple(coo_row_indices[i], coo_col_indices[i], coo_values[i]);
    }

    // 按行优先，然后按列排序 => 保证每行内部 col 有序
    std::sort(triples.begin(), triples.end(),
        [](const std::tuple<size_t, size_t, double>& a,
           const std::tuple<size_t, size_t, double>& b) {
            if (std::get<0>(a) != std::get<0>(b))
                return std::get<0>(a) < std::get<0>(b);  // 按行号升序
            return std::get<1>(a) < std::get<1>(b);      // 同一行，按列号升序
        });

    // ===========================
    // 3. 统计每行非零元个数，构建 row_ptrs（前缀和）
    // ===========================
    std::vector<size_t> row_nnz_counts(rows, 0);
    for (const auto& t : triples) {
        size_t row = std::get<0>(t);
        row_nnz_counts[row]++;
    }

    row_ptrs.resize(rows + 1, 0);
    for (size_t r = 0; r < rows; ++r) {
        row_ptrs[r + 1] = row_ptrs[r] + row_nnz_counts[r];
    }

    // ===========================
    // 4. 填充 values 和 col_indices
    // ===========================
    const size_t nnz_sorted = triples.size();  // == nnz
    this->values.resize(nnz_sorted);
    this->col_indices.resize(nnz_sorted);

    // 每行当前写入的位置
    std::vector<size_t> row_write_positions(rows, 0);

    for (size_t i = 0; i < nnz_sorted; ++i) {
        size_t row = std::get<0>(triples[i]);
        size_t col = std::get<1>(triples[i]);
        double value = std::get<2>(triples[i]);

        size_t start = row_ptrs[row];  // 该行的起始索引
        size_t pos = start + row_write_positions[row];

        this->values[pos] = value;
        this->col_indices[pos] = col;

        row_write_positions[row]++;  // 该行写入位置 +1
    }
}

// 构造函数：从文件构造 CSR 稀疏矩阵
SparseMatrixCSR::SparseMatrixCSR(const std::string& filename) {
	// --- 1. 获取文件后缀名 ---
    size_t dot_pos = filename.find_last_of('.');
    if (dot_pos == std::string::npos || dot_pos == filename.length() - 1) {
        throw std::runtime_error("错误：文件没有后缀名，请提供 .coo 或 .mtx 后缀以检测格式");
    }

    // 提取后缀（例如 ".coo"）
    std::string ext = filename.substr(dot_pos); // 如 ".coo"
    // 转小写（可选，使判断更鲁棒，如 .COO 也能识别）
    for (char &c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));

    // --- 2. 根据后缀名选择处理逻辑 ---
    if (ext == ".coo") {
        // === COO 格式：每行是 row col value ===

        std::ifstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("无法打开 .coo 文件: " + filename);
        }

        std::vector<size_t> coo_rows;
        std::vector<size_t> coo_cols;
        std::vector<double> coo_vals;

        size_t max_row = 0, max_col = 0;

        std::string line;
        while (std::getline(file, line)) {
            // 跳过空行和注释（以 # 或 % 开头）
            if (line.empty() || line[0] == '#' || line[0] == '%') {
                continue;
            }

            std::istringstream iss(line);
            size_t row, col;
            double val;

            if (!(iss >> row >> col >> val)) {
                throw std::runtime_error("文件格式错误，无法解析 COO 行: " + line);
            }

            coo_rows.push_back(row);
            coo_cols.push_back(col);
            coo_vals.push_back(val);

            if (row > max_row) max_row = row;
            if (col > max_col) max_col = col;
        }

        // 推断矩阵的行数和列数（假设行号、列号从 0 开始）
        const size_t num_rows = max_row + 1;
        const size_t num_cols = max_col + 1;

        // 调用已有的 COO -> CSR 构造函数
        *this = SparseMatrixCSR(num_rows, num_cols, coo_rows, coo_cols, coo_vals);

    } else if (ext == ".mtx") {
        // === Matrix Market 格式：未来支持 ===
        throw std::runtime_error("未来支持：Matrix Market 格式 (.mtx)，当前暂未实现");

    } else {
        // === 其他不支持的格式 ===
        throw std::runtime_error("错误：不支持的文件格式 '" + ext + "'，请使用 .coo 或 .mtx 后缀");
    }
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

// 输出矩阵到流
std::ostream& operator<<(std::ostream& os, const SparseMatrixCSR& mat) {
	size_t rows = mat.rows;
	size_t cols = mat.cols;
	const auto& values = mat.values;
	const auto& col_indices = mat.col_indices;
	const auto& row_ptrs = mat.row_ptrs;

	os << "CSR稀疏矩阵 (" << rows << " x " << cols << "), NNZ = " << values.size() << "\n";
	os << "格式: (row, col) -> value\n";

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = row_ptrs[i]; j < row_ptrs[i + 1]; ++j) {
			os << "  (" << i+1 << ", " << col_indices[j]+1 << ") -> " << values[j] << "\n";
		}
	}

	return os;
}

// 输出矩阵到.coo文件
void SparseMatrixCSR::saveToFile(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("无法打开文件用于写入: " + filename);
    }

    // 遍历每一行
    for (size_t r = 0; r < rows; ++r) {
        size_t start = row_ptrs[r];          // 当前行起始位置
        size_t end = row_ptrs[r + 1];        // 当前行结束位置

        // 遍历该行的所有非零元素
        for (size_t pos = start; pos < end; ++pos) {
            size_t col = col_indices[pos];   // 列号
            double val = values[pos];        // 值

            // 写入：行 列 值
            file << r << " " << col << " " << val << "\r\n";
        }
    }
	std::cout << "矩阵已保存到文件: " << filename << std::endl;
    // 文件会在 file 析构时自动关闭
}