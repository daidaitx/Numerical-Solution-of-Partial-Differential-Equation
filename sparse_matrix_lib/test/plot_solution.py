import sys
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 读取数据
# data_name 需要输入进来, 例如 'poisson_2d_1'
if len(sys.argv) < 2:
    print("请在命令行传入数据名称，例如：python script.py poisson_2d_1")
    sys.exit(1)

data_name = sys.argv[1]
file_path = './sparse_matrix_lib/test/data/'
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
	print(f"读取数据失败，请检查数据名称是否正确或文件是否存在。错误信息：{e}")