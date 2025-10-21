#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
from pathlib import Path

def read_vertices(filename):
    """读取多边形顶点数据"""
    vertices = []
    with open(filename, 'r') as f:
        for line in f:
            x, y = map(float, line.strip().split())
            vertices.append([x, y])
    return np.array(vertices)

def read_grid_points(data_dir):
    """读取网格点数据：固定读取 data/grid_X.txt 和 data/grid_Y.txt 文件"""
    x_file = os.path.join(data_dir, 'grid_X.txt')
    y_file = os.path.join(data_dir, 'grid_Y.txt')
    with open(x_file, 'r') as f:
        x_coords = list(map(float, f.read().split()))
    with open(y_file, 'r') as f:
        y_coords = list(map(float, f.read().split()))
    return np.column_stack((x_coords, y_coords))

def plot_polygon(vertices, grid_points, output_file=None, prefix=None):
    """绘制多边形和网格点"""
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # 绘制多边形（使用深青色线条和半透明填充）
    ax.fill(vertices[:, 0], vertices[:, 1], 
            color='#2b8cbe', alpha=0.3)
    ax.plot(vertices[:, 0], vertices[:, 1], 
            color='#045a8d', linewidth=2, linestyle='-')
    ax.plot([vertices[-1, 0], vertices[0, 0]], 
            [vertices[-1, 1], vertices[0, 1]], 
            color='#045a8d', linewidth=2)  # 闭合多边形
    
    # 标记顶点（使用橙色圆点）
    ax.scatter(vertices[:, 0], vertices[:, 1], 
            color='#fd8d3c', s=50, edgecolor='white', 
            linewidth=1, zorder=3)
    
    # 绘制网格点
    ax.scatter(grid_points[:, 0], grid_points[:, 1],
            color='#7f7f7f', s=10, alpha=0.5,
            marker='.', label='Grid Points')
    
    # 设置图形属性
    ax.set_title('Polygon Region with ' + (grid_points.shape[0]).__str__() + ' Grid Points', fontsize=14, pad=20)
    ax.set_xlabel('x', fontsize=12, labelpad=10)
    ax.set_ylabel('y', fontsize=12, labelpad=10)
    ax.grid(True, linestyle=':', color='gray', alpha=0.4)
    ax.set_aspect('equal')
    
    # 自动保存到项目数据文件夹
    if prefix:
        auto_save_path = os.path.join(data_dir, f'{prefix}_grid.png')
        plt.savefig(auto_save_path, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"已自动保存到: {auto_save_path}")
    
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
        print(f"已保存到: {output_file}")

    plt.show()

if __name__ == "__main__":
    # 检查命令行参数
    if len(sys.argv) < 2:
        print("用法: python plot_grid.py <project_name>")
        print("例如: python plot_grid.py polygon_dirichlet_poisson_solver")
        sys.exit(1)
    
    project_name = sys.argv[1]
    
    # 获取当前脚本所在目录
    script_dir = Path(__file__).parent
    
    # 构建项目数据目录路径：上级目录/项目名称/data
    project_data_dir = script_dir.parent / project_name / "data"
    
    # 检查数据目录是否存在
    if not project_data_dir.exists():
        print(f"错误: 数据目录不存在: {project_data_dir}")
        sys.exit(1)
    
    data_dir = str(project_data_dir)
    
    # 数据文件路径
    vertices_file = os.path.join(data_dir, 'polygon_vertices.txt')

    try:
        # 读取数据
        vertices = read_vertices(vertices_file)
        grid_points = read_grid_points(data_dir)

        # 输出文件名为 grid_<N>.png，其中 N 是格点数量
        n_points = grid_points.shape[0]
        output_file = os.path.join(data_dir, f'grid_{n_points}.png')

        # 绘制图形
        plot_polygon(vertices, grid_points, output_file=output_file)
        
    except FileNotFoundError as e:
        print(f"文件未找到错误: {e}")
        print(f"请确保以下文件存在:")
        print(f"  - {vertices_file}")
        print(f"  - {os.path.join(data_dir, 'grid_X.txt')}")
        print(f"  - {os.path.join(data_dir, 'grid_Y.txt')}")
        sys.exit(1)
    except Exception as e:
        print(f"错误: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)