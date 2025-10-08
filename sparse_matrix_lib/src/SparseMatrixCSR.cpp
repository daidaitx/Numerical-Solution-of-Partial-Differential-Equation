#include "SparseMatrixCSR.hpp"
#include <iostream>
#include <vector>
#include <tuple>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <cctype>
#include <stdexcept>
#include <iomanip>
#include <string>
#include <cmath>
#include <unordered_map>

using namespace std;

// 构造函数：构造一个 r x c 全零矩阵
SparseMatrixCSR::SparseMatrixCSR(size_t r, size_t c) : 
	rows(r), cols(c),
	values(), col_indices(),
	row_ptrs(r + 1, 0)  // 初始化行指针为全0
{
	if (r == 0 || c == 0) {
		throw invalid_argument("矩阵行数和列数必须均大于0。");
	}
}

// 构造函数：从稠密矩阵（二维 vector）构造 CSR 稀疏矩阵
SparseMatrixCSR::SparseMatrixCSR(const vector<vector<double>>& dense) {
	// 1. 检查矩阵是否为空
	if (dense.empty() || dense[0].empty()) {
		throw invalid_argument("矩阵行数和列数必须均大于0。");
	}

	// 2. 获取行数和列数
	rows = dense.size();
	cols = dense[0].size();

	// 3. 检查是否是“矩形”矩阵（所有行的列数相同）
	for (const auto& row : dense) {
		if (row.size() != cols) {
			throw invalid_argument("输入的稠密矩阵必须是矩形矩阵，每一行的列数必须相同。");
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
	const vector<size_t>& coo_row_indices,
	const vector<size_t>& coo_col_indices,
	const vector<double>& coo_values)
	: rows(rows), cols(cols)
{
	// ===========================
	// 1. 输入合法性检查
	// ===========================
	if (rows == 0 || cols == 0) {
		throw invalid_argument("错误：矩阵的行数和列数必须均大于 0。");
	}

	const size_t nnz = coo_values.size();
	if (coo_row_indices.size() != nnz || coo_col_indices.size() != nnz) {
		throw invalid_argument("错误：COO 格式错误，行、列、值的数量必须一致。");
	}

	for (size_t i = 0; i < nnz; ++i) {
		if (coo_row_indices[i] >= rows || coo_col_indices[i] >= cols) {
			throw out_of_range("错误：COO 格式错误，行或列索引超出矩阵维度。");
		}
	}

	// ===========================
	// 2. 构造 COO 三元组 (row, col, value) 并排序（先行后列）
	// ===========================
	vector<tuple<size_t, size_t, double>> triples(nnz);
	for (size_t i = 0; i < nnz; ++i) {
		triples[i] = make_tuple(coo_row_indices[i], coo_col_indices[i], coo_values[i]);
	}

	// 按行优先，然后按列排序 => 保证每行内部 col 有序
	sort(triples.begin(), triples.end(),
		[](const tuple<size_t, size_t, double>& a,
		const tuple<size_t, size_t, double>& b) {
			if (get<0>(a) != get<0>(b))
				return get<0>(a) < get<0>(b);  // 按行号升序
			return get<1>(a) < get<1>(b);      // 同一行，按列号升序
		});

	// ===========================
	// 3. 统计每行非零元个数，构建 row_ptrs（前缀和）
	// ===========================
	vector<size_t> row_nnz_counts(rows, 0);
	for (const auto& t : triples) {
		size_t row = get<0>(t);
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
	vector<size_t> row_write_positions(rows, 0);

	for (size_t i = 0; i < nnz_sorted; ++i) {
		size_t row = get<0>(triples[i]);
		size_t col = get<1>(triples[i]);
		double value = get<2>(triples[i]);

		size_t start = row_ptrs[row];  // 该行的起始索引
		size_t pos = start + row_write_positions[row];

		this->values[pos] = value;
		this->col_indices[pos] = col;

		row_write_positions[row]++;  // 该行写入位置 +1
	}
}

// 构造函数：CSR 格式数据构造 CSR 稀疏矩阵
SparseMatrixCSR::SparseMatrixCSR(const size_t rows, const size_t cols,
		const vector<double>& values,
		const vector<size_t>& col_indices,
		const vector<size_t>& row_ptrs) : 
		rows(rows), cols(cols),	values(values), col_indices(col_indices), row_ptrs(row_ptrs) {
	// 检查输入参数
	if (rows == 0) {
		throw invalid_argument("错误：矩阵的行数必须大于 0。");
	}
	if (cols == 0) {
		throw invalid_argument("错误：矩阵的列数必须大于 0。");
	}
	if (values.size() != col_indices.size()) {
		throw invalid_argument("错误：CSR 格式错误，列指标、值的数量必须一致。");
	}
	if (row_ptrs.size() != rows + 1) {
		throw invalid_argument("错误：CSR 格式错误，行指针的数量必须为 rows + 1。");
	}
	for (size_t i = 0; i < values.size(); ++i) {
		if (col_indices[i] >= cols) {
			throw out_of_range("错误：CSR 格式错误，列索引超出矩阵维度。");
		}
	}
	for (size_t i = 0; i < rows; ++i) {
		if (row_ptrs[i] > row_ptrs[i + 1]) {
			throw invalid_argument("错误：CSR 格式错误，行指针必须递增。");
		}
	}
	// 不能包含显式零元素
	for (double val : values) {
		if (val == 0.0) {
			throw invalid_argument("错误：CSR 格式错误，不能包含显式零元素。");
		}
	}
	// 行指针的第一个元素必须为 0
	if (row_ptrs[0] != 0) {
		throw invalid_argument("错误：CSR 格式错误，行指针的第一个元素必须为 0。");
	}
	// 行指针的最后一个元素必须等于 values.size()
	if (row_ptrs[rows] != values.size()) {
		throw invalid_argument("错误：CSR 格式错误，行指针的最后一个元素必须等于非零元素个数。");
	}
	// 检查输入必须无重复、良排序：每一行内的列索引必须严格递增
	for (size_t i = 0; i < rows; ++i) {
		const size_t start = row_ptrs[i];
		const size_t end = row_ptrs[i + 1];
		// 当前行没有非零元素，跳过
		if (start == end) {
			continue;
		}
		// 检查当前行的列索引是否严格递增
		for (size_t j = start + 1; j < end; ++j) {
			if (col_indices[j - 1] >= col_indices[j]) {
				throw invalid_argument("错误：CSR 格式错误，每行内的列索引必须严格递增。");
			}
		}
	}
}

// 构造函数：从文件构造 CSR 稀疏矩阵
SparseMatrixCSR::SparseMatrixCSR(const string& filename) {
	// --- 1. 获取文件后缀名 ---
	size_t dot_pos = filename.find_last_of('.');
	if (dot_pos == string::npos || dot_pos == filename.length() - 1) {
		throw runtime_error("错误：文件没有后缀名，请提供 .coo 或 .mtx 后缀以检测格式");
	}

	// 提取后缀（例如 ".coo"）
	string ext = filename.substr(dot_pos); // 如 ".coo"
	// 转小写（可选，使判断更鲁棒，如 .COO 也能识别）
	for (char &c : ext) c = static_cast<char>(tolower(static_cast<unsigned char>(c)));

	// --- 2. 根据后缀名选择处理逻辑 ---
	if (ext == ".coo") {
		// === COO 格式：每行是 row col value ===

		ifstream file(filename);
		if (!file.is_open()) {
			throw runtime_error("无法打开 .coo 文件: " + filename);
		}

		vector<size_t> coo_rows;
		vector<size_t> coo_cols;
		vector<double> coo_vals;

		size_t max_row = 0, max_col = 0;

		string line;
		while (getline(file, line)) {
			// 跳过空行和注释（以 # 或 % 开头）
			if (line.empty() || line[0] == '#' || line[0] == '%') {
				continue;
			}

			istringstream iss(line);
			size_t row, col;
			double val;

			if (!(iss >> row >> col >> val)) {
				throw runtime_error("文件格式错误，无法解析 COO 行: " + line);
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
		throw runtime_error("未来支持：Matrix Market 格式 (.mtx)，当前暂未实现");

	} else {
		// === 其他不支持的格式 ===
		throw runtime_error("错误：不支持的文件格式 '" + ext + "'，请使用 .coo 或 .mtx 后缀");
	}
}

// getter函数
size_t SparseMatrixCSR::getRows() const noexcept {
	return rows;
}

size_t SparseMatrixCSR::getCols() const noexcept {
	return cols;
}

pair<size_t, size_t> SparseMatrixCSR::getSize() const noexcept {
	return {rows, cols};
}

size_t SparseMatrixCSR::getNNZ() const noexcept {
	return values.size();
}

// 输出矩阵到流
ostream& operator<<(ostream& os, const SparseMatrixCSR& mat) {
	size_t rows = mat.rows;
	size_t cols = mat.cols;
	const auto& values = mat.values;
	const auto& col_indices = mat.col_indices;
	const auto& row_ptrs = mat.row_ptrs;

	os << "CSR稀疏矩阵 (" << rows << " x " << cols << "), NNZ = " << values.size() << "\n";
	os << "格式: (row, col) -> value\n";

	if (values.empty()) {
		os << "（全零矩阵）\n";
		return os;
	}
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = row_ptrs[i]; j < row_ptrs[i + 1]; ++j) {
			os << "  (" << i+1 << ", " << col_indices[j]+1 << ") -> " << values[j] << "\n";
		}
	}

	return os;
}

// 输出矩阵到.coo文件
void SparseMatrixCSR::saveToFile(const string& filename) const {
	ofstream file(filename);
	if (!file.is_open()) {
		throw runtime_error("无法打开文件用于写入: " + filename);
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
			file << r << " " << col << " " << val << "\n";
		}
	}
	cout << "矩阵已保存到文件: " << filename << endl;
	// 文件会在 file 析构时自动关闭
}

// 设置CSR稀疏矩阵元素
void SparseMatrixCSR::setValue(vector<tuple<size_t, size_t, double>> elements) {
	// 1. 根据行列索引对elements升序排列
	sort(elements.begin(), elements.end(),
		[](const tuple<size_t, size_t, double>& a, const tuple<size_t, size_t, double>& b) {
			if (get<0>(a) != get<0>(b)) {
				return get<0>(a) < get<0>(b);  // 按行号升序
			} else if (get<1>(a) != get<1>(b)) {
				return get<1>(a) < get<1>(b);      // 同一行，按列号升序
			} else {
				throw invalid_argument("新插入的元素包含位置完全相同的元素。");
			}
		});
	// 2. 对每一行，插入列指标、值，更新行指针（同时检查输入合法性），然后排序
	if (get<0>(elements[0]) == 0) {
		throw invalid_argument("插入元素的坐标（行）请使用 1-based 记法，不应出现0。");
	}
	if (get<0>(elements[elements.size()-1]) > getRows()) {
		throw invalid_argument("插入元素的坐标（行）超出矩阵维度。");
	}
	size_t elem_id = 0;
	size_t insert_ptr = 0;
	int update_row_ptrs = 0;
	for (size_t row_id = 0; row_id < getRows(); ++row_id) {
		// 更新行指针
		row_ptrs[row_id+1] += update_row_ptrs;
		if (elem_id == elements.size()) {
			continue;
		}
		if (get<1>(elements[elem_id]) == 0) {
			throw invalid_argument("插入元素的坐标（列）请使用 1-based 记法，不应出现0。");
		}
		if (get<1>(elements[elem_id]) > getCols()) {
			throw invalid_argument("插入元素的坐标（列）超出矩阵维度。");
		}
		while (get<0>(elements[elem_id]) == row_id + 1) {
			/**
			 * 现在考虑第 row_id 行：
			 * 要插入的元素是 elements[elem_id]，
			 * 待插入的区间是 col_indices 的 row_ptrs[row_id] 到 row_ptrs[row_id+1) 之间。
			 */
			// 插入/更新/删除元素
			insert_ptr = row_ptrs[row_id];
			while (insert_ptr < row_ptrs[row_id+1] && col_indices[insert_ptr] < get<1>(elements[elem_id])-1) {
				++insert_ptr;
			}
			if (insert_ptr == row_ptrs[row_id+1] || col_indices[insert_ptr] > get<1>(elements[elem_id])-1) {
				// 插入新元素
				if (get<2>(elements[elem_id]) != 0.0) {
					col_indices.insert(col_indices.begin() + insert_ptr, get<1>(elements[elem_id])-1);
					values.insert(values.begin() + insert_ptr, get<2>(elements[elem_id]));
					++row_ptrs[row_id+1];
					++update_row_ptrs;
				}
			} else {
				if (get<2>(elements[elem_id]) != 0.0) {
					// 更新旧元素
					values[insert_ptr] = get<2>(elements[elem_id]);
				} else {
					// 删除旧元素
					col_indices.erase(col_indices.begin() + insert_ptr);
					values.erase(values.begin() + insert_ptr);
					--row_ptrs[row_id+1];
					--update_row_ptrs;
				}
			}
			++elem_id;
			if (elem_id == elements.size()) {
				break;
			}
		}
	}
}

void SparseMatrixCSR::setValue(const tuple<size_t, size_t, double> element) {
	vector<tuple<size_t, size_t, double>> elements = {element};
	setValue(elements);
}

void SparseMatrixCSR::setValue(const vector<size_t>& row_indices,const vector<size_t>& col_indices,const vector<double>& values) {
	if (row_indices.size() != values.size() || col_indices.size() != values.size()) {
		throw invalid_argument("行、列、值的数量必须一致。");
	}
	// 包装为elements
	vector<tuple<size_t, size_t, double>> elements;
	for (size_t i = 0; i < row_indices.size(); ++i) {
		elements.push_back(make_tuple(row_indices[i], col_indices[i], values[i]));
	}
	setValue(elements);
}

void SparseMatrixCSR::setValue(const size_t row, const size_t col, const double value) {
	setValue(make_tuple(row, col, value));
}

// 运算符重载：括号运算
double SparseMatrixCSR::operator()(const size_t row, const size_t col) const {
	if(row == 0) {
		throw invalid_argument("坐标（行）请使用 1-based 记法，不应出现0。");
	}
	if(row > rows) {
		throw invalid_argument("坐标（行）超出矩阵维度。");
	}
	if(col == 0) {
		throw invalid_argument("坐标（列）请使用 1-based 记法，不应出现0。");
	}
	if(col > cols) {
		throw invalid_argument("坐标（列）超出矩阵维度。");
	}
	for(size_t i = row_ptrs[row-1]; i < row_ptrs[row]; ++i) {
		if(col_indices[i] == col-1) {
			return values[i];
		}
	}
	return 0.0;
}

// 转为稠密矩阵
vector<vector<double>> SparseMatrixCSR::toDense() const noexcept {
	vector<vector<double>> dense(rows, vector<double>(cols, 0.0));
	for (size_t i = 0; i < rows; ++i) {
		for (size_t k = row_ptrs[i]; k < row_ptrs[i + 1]; ++k) {
			size_t col = col_indices[k];
			dense[i][col] = values[k];
		}
	}
	return dense;
}

// 矩阵稠密打印
void SparseMatrixCSR::printDense(const int precision) const noexcept {
	size_t max_width = 0;

	// 遍历非零元素 values[]，计算它们格式化后的最大显示宽度
	for (double val : values) {
		ostringstream oss;
		oss << fixed << setprecision(precision);  // 控制小数位数，比如 4 位
		oss << val;
		string s = oss.str();
		if (s.length() > max_width) {
			max_width = s.length();
		}
	}

	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			double val = 0.0;
			for (size_t k = row_ptrs[i]; k < row_ptrs[i + 1]; ++k) {
				if (col_indices[k] == j) {
					val = values[k];
					break;
				}
			}
			ostringstream oss;
			oss << fixed << setprecision(precision);
			oss << val;
			cout << setw(max_width) << oss.str() << " ";
		}
		cout << endl;
	}
}

// 矩阵转置
SparseMatrixCSR SparseMatrixCSR::transpose() const noexcept {
	const size_t num_rows = this->cols;  // 转置后行数 = 原列数
	const size_t num_cols = this->rows;  // 转置后列数 = 原行数
	
	// 如果矩阵为空，直接返回空矩阵
	if (this->values.empty()) {
		return SparseMatrixCSR(num_rows, num_cols);
	}

	// Step 1: 统计原矩阵中每一列（即转置后每一行）有多少非零元素
	vector<size_t> col_counts(num_rows, 0);  // 大小为转置后行数 = 原列数
	for (size_t k = 0; k < this->values.size(); ++k) {
		size_t col = this->col_indices[k];
		col_counts[col]++;
	}

	// Step 2: 计算转置后的 row_ptrs（即列的前缀和）
	vector<size_t> transposed_row_ptrs(num_rows + 1, 0);
	for (size_t i = 0; i < num_rows; ++i) {
		transposed_row_ptrs[i + 1] = transposed_row_ptrs[i] + col_counts[i];
	}

	// Step 3: 填充转置后的 values, col_indices
	vector<double> transposed_values(this->values.size());
	vector<size_t> transposed_col_indices(this->col_indices.size());

	// 辅助数组：记录每一列（转置行）当前写入位置
	vector<size_t> current_pos(num_rows, 0);

	for (size_t i = 0; i < this->rows; ++i) {
		for (size_t k = this->row_ptrs[i]; k < this->row_ptrs[i + 1]; ++k) {
			size_t col = this->col_indices[k];       // 原列号 → 转置后行号
			size_t pos = transposed_row_ptrs[col] + current_pos[col];
			transposed_values[pos] = this->values[k];
			transposed_col_indices[pos] = i;  // 原行号 → 转置后列号
			current_pos[col]++;
		}
	}

	// Step 4: 构造并返回新的 CSR（转置后矩阵）
	return SparseMatrixCSR(num_rows, num_cols, 
						transposed_values, 
						transposed_col_indices, 
						transposed_row_ptrs);
}

// 矩阵范数
double SparseMatrixCSR::norm(const string& norm_type) const {
	if (norm_type == "row" || norm_type == "inf" || norm_type == "infinity" || norm_type == "INF" || norm_type == "INFINITY") {
		// 行和范数（每行绝对值之和的最大值）→ ∞-范数
		double max_row_sum = 0.0;
		for (size_t i = 0; i < rows; ++i) {
			double row_sum = 0.0;
			for (size_t k = row_ptrs[i]; k < row_ptrs[i + 1]; ++k) {
				row_sum += fabs(values[k]); // 累加绝对值
			}
			max_row_sum = max(max_row_sum, row_sum);
		}
		return max_row_sum;
	} else if (norm_type == "col" || norm_type == "column" || norm_type == "1" || norm_type == "one" || norm_type == "ONE") {
		// 列和范数（每列绝对值之和的最大值）→ 1-范数
		vector<double> col_sums(cols, 0.0); // 每列的绝对值之和
		for (size_t i = 0; i < rows; ++i) {
			for (size_t k = row_ptrs[i]; k < row_ptrs[i + 1]; ++k) {
				size_t col = col_indices[k];
				col_sums[col] += fabs(values[k]);
			}
		}
		double max_col_sum = 0.0;
		for (double sum : col_sums) {
			max_col_sum = max(max_col_sum, sum);
		}
		return max_col_sum;
	} else if (norm_type == "fro" || norm_type == "Fro" || norm_type == "FRO" || norm_type == "frobenius" || norm_type == "Frobenius" || norm_type == "FROBENIUS") {
		// Frobenius 范数：sqrt( sum(|A_ij|^2) )
		double sum_sq = 0.0;
		for (double val : values) {
			sum_sq += val * val; // 已经是非负，不需要 fabs
		}
		return sqrt(sum_sq);
	} else if (norm_type == "max" || norm_type == "maximum" || norm_type == "MAX" || norm_type == "MAXIMUM") {
		// 按元素最大范数：max( |A_ij| )
		double max_val = 0.0;
		for (double val : values) {
			max_val = max(max_val, fabs(val));
		}
		return max_val;
	} else {
		throw invalid_argument("不支持的范数类型：" + norm_type);
	}
}

// 判断矩阵是否为方阵
bool SparseMatrixCSR::isSquare() const noexcept {
	return this->rows == this->cols;
}

// 获取矩阵的对角元素，输出为向量
vector<double> SparseMatrixCSR::getDiagonalVector() const noexcept {
	const size_t n = min(this->rows, this->cols);
	vector<double> diagonal(n, 0.0);
	
	for (size_t i = 0; i < n; ++i) {
		// 在当前行中查找对角线元素 (i, i)
		for (size_t k = this->row_ptrs[i]; k < this->row_ptrs[i + 1]; ++k) {
			if (this->col_indices[k] == i) {
				diagonal[i] = this->values[k];
				break;
			}
		}
	}
	
	return diagonal;
}

// 获取矩阵的对角元素，输出为CSR格式稀疏对角矩阵
SparseMatrixCSR SparseMatrixCSR::getDiagonalMatrix() const noexcept {
	const size_t n = min(this->rows, this->cols);
	vector<double> diag_values;
	vector<size_t> diag_col_indices;
	vector<size_t> diag_row_ptrs(n + 1, 0);
	
	// 收集对角线元素
	for (size_t i = 0; i < n; ++i) {
		for (size_t k = this->row_ptrs[i]; k < this->row_ptrs[i + 1]; ++k) {
			if (this->col_indices[k] == i) {
				diag_values.push_back(this->values[k]);
				diag_col_indices.push_back(i);
				break;
			}
		}
		diag_row_ptrs[i + 1] = diag_values.size();
	}
	
	return SparseMatrixCSR(n, n, diag_values, diag_col_indices, diag_row_ptrs);
}

// 获取矩阵的上三角部分
SparseMatrixCSR SparseMatrixCSR::getUpperTriangularMatrix() const {
	if (!this->isSquare()) {
		throw invalid_argument("要获取矩阵的上三角部分，需要矩阵是一个方阵");
	}
	
	vector<double> upper_values;
	vector<size_t> upper_col_indices;
	vector<size_t> upper_row_ptrs(this->rows + 1, 0);
	
	for (size_t i = 0; i < this->rows; ++i) {
		for (size_t k = this->row_ptrs[i]; k < this->row_ptrs[i + 1]; ++k) {
			size_t j = this->col_indices[k];
			// 只保留对角线及上三角元素 (j >= i)
			if (j >= i) {
				upper_values.push_back(this->values[k]);
				upper_col_indices.push_back(j);
			}
		}
		upper_row_ptrs[i + 1] = upper_values.size();
	}
	
	return SparseMatrixCSR(this->rows, this->cols, upper_values, upper_col_indices, upper_row_ptrs);
}

// 获取矩阵的下三角部分
SparseMatrixCSR SparseMatrixCSR::getLowerTriangularMatrix() const {
	if (!this->isSquare()) {
		throw invalid_argument("要获取矩阵的下三角部分，需要矩阵是一个方阵");
	}
	
	vector<double> lower_values;
	vector<size_t> lower_col_indices;
	vector<size_t> lower_row_ptrs(this->rows + 1, 0);
	
	for (size_t i = 0; i < this->rows; ++i) {
		for (size_t k = this->row_ptrs[i]; k < this->row_ptrs[i + 1]; ++k) {
			size_t j = this->col_indices[k];
			// 只保留对角线及下三角元素 (j <= i)
			if (j <= i) {
				lower_values.push_back(this->values[k]);
				lower_col_indices.push_back(j);
			}
		}
		lower_row_ptrs[i + 1] = lower_values.size();
	}
	
	return SparseMatrixCSR(this->rows, this->cols, lower_values, lower_col_indices, lower_row_ptrs);
}

// 判断矩阵是否对称
bool SparseMatrixCSR::isSymmetric(const double tolerance) const {
	if (!this->isSquare()) {
		return false;
	}
	// 构建转置矩阵进行比较
	SparseMatrixCSR transposed = this->transpose();
	// 调用 isApproxEqualto 进行比较
	return isApproxEqualto(transposed, tolerance);
}

// 判断矩阵是否是上三角矩阵
bool SparseMatrixCSR::isUpperTriangular() const noexcept {
	if (!this->isSquare()) {
		return false;
	}
	
	for (size_t i = 0; i < this->rows; ++i) {
		for (size_t k = this->row_ptrs[i]; k < this->row_ptrs[i + 1]; ++k) {
			size_t j = this->col_indices[k];
			// 如果发现下三角区域（j < i）有非零元素，则不是上三角矩阵
			if (j < i) {
				return false;
			}
		}
	}
	
	return true;
}

// 判断矩阵是否是下三角矩阵
bool SparseMatrixCSR::isLowerTriangular() const noexcept {
	if (!this->isSquare()) {
		return false;
	}
	
	for (size_t i = 0; i < this->rows; ++i) {
		for (size_t k = this->row_ptrs[i]; k < this->row_ptrs[i + 1]; ++k) {
			size_t j = this->col_indices[k];
			// 如果发现上三角区域（j > i）有非零元素，则不是下三角矩阵
			if (j > i) {
				return false;
			}
		}
	}
	
	return true;
}

// 判断矩阵是否是对角矩阵
bool SparseMatrixCSR::isDiagonal() const noexcept {
	if (!this->isSquare()) {
		return false;
	}
	
	for (size_t i = 0; i < this->rows; ++i) {
		for (size_t k = this->row_ptrs[i]; k < this->row_ptrs[i + 1]; ++k) {
			size_t j = this->col_indices[k];
			// 如果发现非对角线元素（j != i），则不是对角矩阵
			if (j != i) {
				return false;
			}
		}
	}
	
	return true;
}

// 判断两个矩阵是否近似相等
bool SparseMatrixCSR::isApproxEqualto(const SparseMatrixCSR& other, const double tolerance) const {
	SparseMatrixCSR DiffMat = *this - other;
	return DiffMat.norm("max") <= tolerance * (other.norm("max") + this->norm("max") + 1e-10);
}

// 矩阵加法(+=)
SparseMatrixCSR& SparseMatrixCSR::operator+=(const SparseMatrixCSR& other) {
	// 检查维度是否匹配
	if (this->rows != other.rows || this->cols != other.cols) {
		throw invalid_argument("矩阵相加要求两矩阵大小相同");
	}
	
	// 临时存储结果
	vector<double> result_values;
	vector<size_t> result_col_indices;
	vector<size_t> result_row_ptrs(this->rows + 1, 0);
	
	// 逐行合并
	for (size_t i = 0; i < this->rows; ++i) {
		size_t pos1 = this->row_ptrs[i];  // 当前矩阵当前行的起始位置
		size_t pos2 = other.row_ptrs[i];  // 另一个矩阵当前行的起始位置
		size_t end1 = this->row_ptrs[i + 1];  // 当前矩阵当前行的结束位置
		size_t end2 = other.row_ptrs[i + 1];  // 另一个矩阵当前行的结束位置
		
		// 归并两个有序列表（按列索引）
		while (pos1 < end1 && pos2 < end2) {
			size_t col1 = this->col_indices[pos1];
			size_t col2 = other.col_indices[pos2];
			
			if (col1 < col2) {
				// 只存在于当前矩阵
				result_values.push_back(this->values[pos1]);
				result_col_indices.push_back(col1);
				pos1++;
			} else if (col1 > col2) {
				// 只存在于另一个矩阵
				result_values.push_back(other.values[pos2]);
				result_col_indices.push_back(col2);
				pos2++;
			} else {
				// 两个矩阵都有，相加
				double sum = this->values[pos1] + other.values[pos2];
				if (abs(sum) > 1e-12) {
					result_values.push_back(sum);
					result_col_indices.push_back(col1);
				}
				pos1++;
				pos2++;
			}
		}
		
		// 处理剩余元素（当前矩阵）
		while (pos1 < end1) {
			result_values.push_back(this->values[pos1]);
			result_col_indices.push_back(this->col_indices[pos1]);
			pos1++;
		}
		
		// 处理剩余元素（另一个矩阵）
		while (pos2 < end2) {
			result_values.push_back(other.values[pos2]);
			result_col_indices.push_back(other.col_indices[pos2]);
			pos2++;
		}
		
		result_row_ptrs[i + 1] = result_values.size();
	}
	
	// 更新当前矩阵
	this->values = move(result_values);
	this->col_indices = move(result_col_indices);
	this->row_ptrs = move(result_row_ptrs);
	
	return *this;
}

// 矩阵加法(+)
SparseMatrixCSR SparseMatrixCSR::operator+(const SparseMatrixCSR& other) const {
	SparseMatrixCSR result = *this;  // 拷贝当前矩阵
	result += other;                 // 使用 += 实现加法
	return result;
}

// 矩阵减法(-=)
SparseMatrixCSR& SparseMatrixCSR::operator-=(const SparseMatrixCSR& other) {
	// 检查维度是否匹配
	if (this->rows != other.rows || this->cols != other.cols) {
		throw invalid_argument("矩阵相减要求两矩阵大小相同");
	}
	
	// 临时存储结果
	vector<double> result_values;
	vector<size_t> result_col_indices;
	vector<size_t> result_row_ptrs(this->rows + 1, 0);
	
	// 逐行合并
	for (size_t i = 0; i < this->rows; ++i) {
		size_t pos1 = this->row_ptrs[i];  // 当前矩阵当前行的起始位置
		size_t pos2 = other.row_ptrs[i];  // 另一个矩阵当前行的起始位置
		size_t end1 = this->row_ptrs[i + 1];  // 当前矩阵当前行的结束位置
		size_t end2 = other.row_ptrs[i + 1];  // 另一个矩阵当前行的结束位置
		
		// 归并两个有序列表（按列索引）
		while (pos1 < end1 && pos2 < end2) {
			size_t col1 = this->col_indices[pos1];
			size_t col2 = other.col_indices[pos2];
			
			if (col1 < col2) {
				// 只存在于当前矩阵
				result_values.push_back(this->values[pos1]);
				result_col_indices.push_back(col1);
				pos1++;
			} else if (col1 > col2) {
				// 只存在于另一个矩阵，取负
				result_values.push_back(-other.values[pos2]);
				result_col_indices.push_back(col2);
				pos2++;
			} else {
				// 两个矩阵都有，相减
				double diff = this->values[pos1] - other.values[pos2];
				if (abs(diff) > 1e-12) {
					result_values.push_back(diff);
					result_col_indices.push_back(col1);
				}
				pos1++;
				pos2++;
			}
		}
		
		// 处理剩余元素（当前矩阵）
		while (pos1 < end1) {
			result_values.push_back(this->values[pos1]);
			result_col_indices.push_back(this->col_indices[pos1]);
			pos1++;
		}
		
		// 处理剩余元素（另一个矩阵），取负
		while (pos2 < end2) {
			result_values.push_back(-other.values[pos2]);
			result_col_indices.push_back(other.col_indices[pos2]);
			pos2++;
		}
		
		result_row_ptrs[i + 1] = result_values.size();
	}
	
	// 更新当前矩阵
	this->values = move(result_values);
	this->col_indices = move(result_col_indices);
	this->row_ptrs = move(result_row_ptrs);
	
	return *this;
}

// 矩阵减法(-)
SparseMatrixCSR SparseMatrixCSR::operator-(const SparseMatrixCSR& other) const {
	SparseMatrixCSR result = *this;  // 拷贝当前矩阵
	result -= other;                 // 使用 -= 实现减法
	return result;
}

// 矩阵取负(单目-)
SparseMatrixCSR SparseMatrixCSR::operator-() const {
	SparseMatrixCSR result = *this;
	for (double& value : result.values) {
		value = -value;
	}
	return result;
}

// 矩阵数乘(*)
SparseMatrixCSR SparseMatrixCSR::operator*(const double scalar) const noexcept {
	// 如果标量为0，返回空矩阵
	if (abs(scalar) == 0.0) {
		return SparseMatrixCSR(this->rows, this->cols);
	}
	
	// 创建当前矩阵的副本
	SparseMatrixCSR result = *this;
	
	// 将所有非零元素乘以标量
	for (double& value : result.values) {
		value *= scalar;
	}
	
	return result;
}

// 矩阵数除(/)
SparseMatrixCSR SparseMatrixCSR::operator/(const double scalar) const {
	// 检查除数是否为0
	if (abs(scalar) == 0.0) {
		throw invalid_argument("不能除以0");
	}
	
	// 创建当前矩阵的副本
	SparseMatrixCSR result = *this;
	
	// 将所有非零元素除以标量
	for (double& value : result.values) {
		value /= scalar;
	}
	
	return result;
}

// CSR稀疏矩阵-向量乘法(A * vec)
vector<double> SparseMatrixCSR::operator*(const vector<double>& other) const {
	// 检查维度匹配：矩阵列数必须等于向量长度
	if (this->cols != other.size()) {
		throw invalid_argument("CSR矩阵列数必须等于向量长度");
	}
	
	vector<double> result(this->rows, 0.0);
	
	// 对于每一行，计算与向量的点积
	for (size_t i = 0; i < this->rows; ++i) {
		double sum = 0.0;
		// 遍历该行的所有非零元素
		for (size_t k = this->row_ptrs[i]; k < this->row_ptrs[i + 1]; ++k) {
			size_t j = this->col_indices[k];
			sum += this->values[k] * other[j];
		}
		result[i] = sum;
	}
	
	return result;
}

// CSR稀疏矩阵-稠密矩阵乘法(A * dense_mat)
vector<vector<double>> SparseMatrixCSR::operator*(const vector<vector<double>>& other) const {
	// 检查维度匹配
	if (this->cols != other.size()) {
		throw invalid_argument("CSR矩阵列数必须等于稠密矩阵行数");
	}
	
	if (other.empty()) {
		return {};
	}
	
	// 检查稠密矩阵的列一致性
	const size_t other_cols = other[0].size();
	for (size_t i = 1; i < other.size(); ++i) {
		if (other[i].size() != other_cols) {
			throw invalid_argument("稠密矩阵的列数必须一致");
		}
	}
	
	// 结果矩阵：this->rows × other_cols
	vector<vector<double>> result(this->rows, vector<double>(other_cols, 0.0));
	
	// 高性能实现：单次遍历稀疏矩阵
	for (size_t i = 0; i < this->rows; ++i) {
		// 对于当前行的每个非零元素
		for (size_t k = this->row_ptrs[i]; k < this->row_ptrs[i + 1]; ++k) {
			size_t j = this->col_indices[k];
			double value = this->values[k];
			
			// 更新结果矩阵的当前行
			for (size_t col = 0; col < other_cols; ++col) {
				result[i][col] += value * other[j][col];
			}
		}
	}
	
	return result;
}