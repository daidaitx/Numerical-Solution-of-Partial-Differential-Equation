#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
import os

def resolve_file(arg: str, data_dir: Path) -> Path:
    """在 data_dir 中按文件名查找文件。"""
    if '/' in arg or '\\' in arg or Path(arg).is_absolute():
        raise ValueError(f"参数必须是文件名而不是路径：{arg}")

    if not data_dir.exists() or not data_dir.is_dir():
        raise FileNotFoundError(f"数据目录不存在: {data_dir}")

    # 先尝试完全匹配文件名（包括扩展名）
    for f in data_dir.iterdir():
        if f.is_file() and f.name == arg:
            return f

    # 再尝试按 stem 匹配（去掉扩展名的名称）
    for f in data_dir.iterdir():
        if f.is_file() and f.stem == arg:
            return f

    raise FileNotFoundError(f"在 {data_dir} 中找不到名为 '{arg}' 的文件")

def read_coo_matrix(file_path):
    """读取COO格式的稀疏矩阵文件"""
    with open(file_path, 'r') as f:
        # 读取第一行：行数、列数、非零元数
        first_line = f.readline().strip().split()
        rows = int(first_line[0])
        cols = int(first_line[1])
        nnz = int(first_line[2])
        
        # 读取非零元素数据
        row_indices = []
        col_indices = []
        values = []
        
        for line in f:
            data = line.strip().split()
            if len(data) == 3:
                # 将1-based索引转换为0-based
                row_indices.append(int(data[0]) - 1)
                col_indices.append(int(data[1]) - 1)
                values.append(float(data[2]))
    
    return rows, cols, nnz, row_indices, col_indices, values

def plot_sparse_matrix(rows, cols, row_indices, col_indices, output_path, title=None):
    """绘制稀疏矩阵的结构图"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # 这样左上角对应矩阵的(0,0)元素，右下角对应(rows-1, cols-1)元素
    ax.scatter(col_indices, row_indices, color='#1f77b4', s=0.8, marker='s', alpha=0.7)
    
    # 设置图形属性
    ax.set_xlabel('Column Index')
    ax.set_ylabel('Row Index')
    
    if title:
        ax.set_title(title)
    else:
        # 移除输出文件名中的.png扩展名用于标题
        title_name = os.path.basename(output_path).replace('.png', '')
        ax.set_title(f'Sparse Matrix Structure: {title_name}\nSize: {rows} × {cols}, Non-zeros: {len(row_indices)}')
    
    # 设置坐标轴范围
    ax.set_xlim(-0.5, cols - 0.5)
    ax.set_ylim(-0.5, rows - 0.5)
    
    # 设置网格
    ax.grid(True, alpha=0.3)
    
    # 设置坐标轴刻度
    max_ticks = 20  # 最大刻度数量
    if rows <= max_ticks:
        ax.set_yticks(np.arange(0, rows))
    else:
        step = rows // max_ticks
        ax.set_yticks(np.arange(0, rows, step))
    
    if cols <= max_ticks:
        ax.set_xticks(np.arange(0, cols))
    else:
        step = cols // max_ticks
        ax.set_xticks(np.arange(0, cols, step))
    
    ax.invert_yaxis()  # 反转y轴，使得(0,0)在左上角

    # 设置等比例，确保矩阵显示为正方形
    ax.set_aspect('equal')
    
    plt.tight_layout()
    
    # 保存图片到指定路径
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"图片已保存到: {output_path}")
    
    # 显示图片
    plt.show()

def main():
    if len(sys.argv) < 4:
        print("用法: python spy_CSR_matrix.py project_name matrix_file output_image_name")
        print("参数说明:")
        print("  project_name      - 项目名称（如: polygon_dirichlet_poisson_solver）")
        print("  matrix_file       - 矩阵文件名（如: matrix_A 或 matrix_A.coo）")
        print("  output_image_name - 输出图片名称（如: spy_matrix_A_506 或 spy_matrix_A_506.png）")
        sys.exit(1)
    
    # 解析命令行参数
    project_name = sys.argv[1]
    matrix_file = sys.argv[2]
    output_image_name = sys.argv[3]
    
    # 构建项目数据目录路径
    script_dir = Path(__file__).parent
    project_data_dir = script_dir.parent / project_name / "data"
    
    # 检查数据目录是否存在
    if not project_data_dir.exists():
        print(f"错误: 数据目录不存在: {project_data_dir}")
        sys.exit(1)
    
    data_dir = project_data_dir
    
    # 解析矩阵文件路径
    matrix_path = resolve_file(matrix_file, data_dir)
    
    # 确保输出文件名有.png扩展名
    if not output_image_name.lower().endswith('.png'):
        output_image_name += '.png'
    
    # 构建输出文件路径
    output_path = data_dir / output_image_name
    
    try:
        # 读取矩阵数据
        rows, cols, nnz, row_indices, col_indices, values = read_coo_matrix(matrix_path)
        
        print(f"矩阵信息:")
        print(f"  大小: {rows} × {cols}")
        print(f"  非零元素: {nnz}")
        print(f"  稀疏度: {nnz/(rows*cols)*100:.6f}%")
        
        # 绘制稀疏矩阵
        plot_sparse_matrix(rows, cols, row_indices, col_indices, output_path)
        
    except FileNotFoundError as e:
        print(f"文件未找到错误: {e}")
        print(f"请确保以下文件存在于 {data_dir}:")
        print(f"  - {matrix_file} 或 {matrix_file}.coo")
        sys.exit(1)
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()