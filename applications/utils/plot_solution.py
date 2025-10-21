import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import matplotlib.tri as tri

def segments_intersect(a, b, c, d, tol=1e-10):
    """检查线段ab和cd是否相交，考虑数值误差"""
    def ccw(A, B, C):
        return (C[1]-A[1]) * (B[0]-A[0]) - (B[1]-A[1]) * (C[0]-A[0])
    
    def point_on_segment(p, a, b):
        """检查点p是否在线段ab上（考虑容忍度）"""
        # 检查点p是否在a和b的边界框内
        if (min(a[0], b[0]) - tol <= p[0] <= max(a[0], b[0]) + tol and
            min(a[1], b[1]) - tol <= p[1] <= max(a[1], b[1]) + tol):
            # 检查点p是否在直线ab上
            cross = ccw(a, b, p)
            return abs(cross) <= tol
        return False
    
    # 检查端点是否在另一条线段上
    if point_on_segment(a, c, d) or point_on_segment(b, c, d) or \
       point_on_segment(c, a, b) or point_on_segment(d, a, b):
        return True
    
    # 使用改进的跨立试验，考虑数值误差
    o1 = ccw(a, c, d)
    o2 = ccw(b, c, d)
    o3 = ccw(c, a, b)
    o4 = ccw(d, a, b)
    
    # 如果两条线段严格相交
    if (o1 * o2 < -tol) and (o3 * o4 < -tol):
        return True
    
    # 如果两条线段几乎共线，检查它们是否重叠
    if abs(o1) <= tol and abs(o2) <= tol and abs(o3) <= tol and abs(o4) <= tol:
        # 检查线段是否在x轴或y轴上有重叠
        if (max(a[0], b[0]) >= min(c[0], d[0]) - tol and
            min(a[0], b[0]) <= max(c[0], d[0]) + tol and
            max(a[1], b[1]) >= min(c[1], d[1]) - tol and
            min(a[1], b[1]) <= max(c[1], d[1]) + tol):
            return True
    
    return False

def resolve_file(arg: str, data_dir: Path) -> Path:
    """仅在 data_dir 中按文件名查找文件。"""
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

def read_polygon_vertices(data_dir):
    """读取多边形顶点"""
    vertices_file = data_dir / "polygon_vertices.txt"
    vertices = []
    with open(vertices_file, 'r') as f:
        for line in f:
            x, y = map(float, line.strip().split())
            vertices.append([x, y])
    return np.array(vertices)

if len(sys.argv) < 6:
    print("用法: python plot_solution.py project_name X_file Y_file U_file output_name")
    print("参数说明:")
    print("  project_name - 项目名称（如: polygon_dirichlet_poisson_solver）")
    print("  X_file       - X坐标文件名（如: grid_X 或 grid_X.txt）")
    print("  Y_file       - Y坐标文件名（如: grid_Y 或 grid_Y.txt）")
    print("  U_file       - 函数值文件名（如: U_exact 或 U_exact.txt）")
    print("  output_name  - 输出图片名称（如: exact_solution_506 或 exact_solution_506.png）")
    sys.exit(1)

# 解析命令行参数
project_name = sys.argv[1]
arg_X = sys.argv[2]
arg_Y = sys.argv[3]
arg_U = sys.argv[4]
output_name = sys.argv[5]

# 构建项目数据目录路径
script_dir = Path(__file__).parent
project_data_dir = script_dir.parent / project_name / "data"

# 检查数据目录是否存在
if not project_data_dir.exists():
    print(f"错误: 数据目录不存在: {project_data_dir}")
    sys.exit(1)

data_dir = project_data_dir

# 解析文件路径
X_path = resolve_file(arg_X, data_dir)
Y_path = resolve_file(arg_Y, data_dir)
U_path = resolve_file(arg_U, data_dir)

# 确保输出文件名有.png扩展名
if not output_name.lower().endswith('.png'):
    output_name += '.png'

try:
    # 读取数据
    X = np.loadtxt(str(X_path))
    Y = np.loadtxt(str(Y_path))
    U = np.loadtxt(str(U_path))
    
    # 读取多边形顶点
    polygon_vertices = read_polygon_vertices(data_dir)
    
    # 创建三角剖分
    triangles = tri.Triangulation(X, Y)
    
    # 创建多边形边界线段
    polygon_edges = []
    n_vertices = len(polygon_vertices)
    for i in range(n_vertices):
        start = polygon_vertices[i]
        end = polygon_vertices[(i+1) % n_vertices]
        polygon_edges.append((start, end))
    
    # 创建掩码：如果三角形的任何边与多边形边界相交，则掩码掉该三角形
    mask = np.zeros(len(triangles.triangles), dtype=bool)
    for i, tri_indices in enumerate(triangles.triangles):
        # 获取三角形的三个顶点
        A = (X[tri_indices[0]], Y[tri_indices[0]])
        B = (X[tri_indices[1]], Y[tri_indices[1]])
        C = (X[tri_indices[2]], Y[tri_indices[2]])
        
        # 三角形的三条边
        triangle_edges = [(A, B), (B, C), (C, A)]
        
        # 检查三角形的每条边是否与多边形边界相交
        intersects_boundary = False
        for tri_edge in triangle_edges:
            for poly_edge in polygon_edges:
                if segments_intersect(tri_edge[0], tri_edge[1], poly_edge[0], poly_edge[1]):
                    intersects_boundary = True
                    break
            if intersects_boundary:
                break
        
        # 如果三角形与边界相交，则掩码掉
        if intersects_boundary:
            mask[i] = True
    
    # 应用掩码
    triangles.set_mask(mask)
    
    # 创建包含两个子图的图形
    fig = plt.figure(figsize=(16, 7))
    
    # 第一个子图：3D图
    ax1 = fig.add_subplot(121, projection='3d')
    
    # 使用掩码后的三角剖分绘制曲面
    surf = ax1.plot_trisurf(triangles, U, cmap='viridis', alpha=0.9, edgecolor='none')
    
    # 在多边形边界处绘制轮廓
    polygon_closed = np.vstack([polygon_vertices, polygon_vertices[0]])  # 闭合多边形
    ax1.plot(polygon_closed[:, 0], polygon_closed[:, 1], 
            np.zeros_like(polygon_closed[:, 0]), 
            'r-', linewidth=2, label='Polygon Boundary')
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('u(x,y)')

    # 移除输出文件名中的.png扩展名用于标题
    title_name = output_name.replace('.png', '') if output_name.endswith('.png') else output_name

    ax1.set_title(f'3D Solution: {title_name}')
    ax1.legend()
    
    # 添加颜色条
    fig.colorbar(surf, ax=ax1, shrink=0.5, aspect=5)
    
    # 第二个子图：俯视图
    ax2 = fig.add_subplot(122)
    
    # 绘制俯视图（等高线图）
    contour = ax2.tricontourf(triangles, U, cmap='viridis', levels=20)
    ax2.tricontour(triangles, U, colors='black', linewidths=0.5, levels=10)
    
    # 绘制多边形边界
    ax2.plot(polygon_closed[:, 0], polygon_closed[:, 1], 
            'r-', linewidth=2, label='Polygon Boundary')
    
    # 标记网格点
    ax2.scatter(X, Y, c='black', s=1, alpha=0.5, marker='.')
    
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_title(f'Top View (Contour): {title_name}')
    ax2.set_aspect('equal')

    # 将图例放在图像外部右上角
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    # 添加颜色条
    fig.colorbar(contour, ax=ax2, shrink=0.5, aspect=5)
    
    plt.tight_layout()
    
    # 使用指定的输出文件名保存图片到项目数据目录
    out_path = data_dir / output_name
    plt.savefig(str(out_path), dpi=300, bbox_inches='tight')
    print(f"图片已保存到: {out_path}")
    plt.show()

except FileNotFoundError as e:
    print(f"文件未找到错误: {e}")
    print(f"请确保以下文件存在于 {data_dir}:")
    print(f"  - polygon_vertices.txt")
    print(f"  - {arg_X} 或 {arg_X}.txt")
    print(f"  - {arg_Y} 或 {arg_Y}.txt")
    print(f"  - {arg_U} 或 {arg_U}.txt")
    sys.exit(1)
except Exception as e:
    print(f"错误: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)