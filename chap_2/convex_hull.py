# 在 $\mathbb{R}^2$ 中随机生成 10 个点，用 scipy.spatial.ConvexHull 计算其凸包，并绘制：
#   - 原始点（散点图）
#   - 凸包的边界（多边形）
#   - 用不同颜色区分凸包顶点和内部点

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

# 生成随机数
np.random.seed(66)
points = np.random.rand(10, 2)
hull = ConvexHull(points) # 返回points arr中组成凸包的顶点索引（二维以逆时针方向返回）
hull_points = points[hull.vertices] # 这里返回的是 numpy 的多维数组

# 差集找非凸包顶点的点
vertex_indices = set(hull.vertices)
all_indices = set(range(len(points)))
interior_indices = list(all_indices - vertex_indices)
interior_points = points[interior_indices]

# 绘图
plt.figure(figsize=(8, 6))
# 画点
plt.scatter(interior_points[:, 0], interior_points[:, 1], color='blue', label='Interior points')
plt.scatter(hull_points[:, 0], hull_points[:, 1], color='red', label='Hull vertices')
# 画线
hull_points_closed = np.vstack([hull_points, hull_points[0]]) # Vertical Stack 拼接起始点到最后，让plot画一个封闭曲线
plt.plot(hull_points_closed[:, 0], hull_points_closed[:, 1],
         color='green', linewidth=2, label='Convex hull boundary')
# 补全坐标图及图例
plt.xlabel('x')
plt.ylabel('y')
plt.title('Convex Hull of 10 Random Points')
plt.legend()
plt.grid(True)
plt.axis('equal') # 保证坐标轴单位长度一致
plt.show()