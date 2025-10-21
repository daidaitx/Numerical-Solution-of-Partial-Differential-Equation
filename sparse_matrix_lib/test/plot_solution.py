import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# 请在当前文件所在目录下的 data 文件夹中放入数据文件，文件名格式为：data_name_X.mat, data_name_Y.mat, data_name_solution.mat

# 读取数据
if len(sys.argv) < 2:
    print("请在命令行传入数据名称，例如：python path/to/plot_solution.py poisson_2d_1")
    sys.exit(1)

data_name = sys.argv[1]
file_path = os.path.join(os.path.dirname(__file__), 'data') + os.sep
file_path_name = file_path + data_name
try:
	U = np.loadtxt(file_path_name + '_solution.mat')
	X = np.loadtxt(file_path_name + '_X.mat')
	Y = np.loadtxt(file_path_name + '_Y.mat')

	# 绘制3D图
	fig = plt.figure(figsize=(10, 8))
	ax = fig.add_subplot(111, projection='3d')
	surf = ax.plot_surface(X, Y, U, cmap='viridis', alpha=0.9)

	ax.set_xlabel('x')
	ax.set_ylabel('y')
	ax.set_zlabel('u(x,y)')
	ax.set_title('2D Poisson Solution')

	fig.colorbar(surf, ax=ax, shrink=0.5, aspect=5)

	plt.tight_layout()
	plt.savefig(file_path_name + '_solution.png', dpi=300)
	plt.show()

except Exception as e:
	print(f"绘制图形失败。错误信息：{e}")