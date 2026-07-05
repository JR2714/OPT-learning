# Deepseek 

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

# 1. 准备固定的网格数据
x = np.linspace(-1.2, 1.2, 500)
y = np.linspace(-1.2, 1.2, 500)
X, Y = np.meshgrid(x, y)

# 2. 创建画布和初始绘图
fig, ax = plt.subplots(figsize=(6, 6))
plt.subplots_adjust(bottom=0.2)             # 为滑块留出空间

p_init = 1.0
Z = np.abs(X)**p_init + np.abs(Y)**p_init
cont = ax.contour(X, Y, Z, levels=[1], colors='blue', linewidths=2)

ax.set_aspect('equal')
ax.set_xlim(-1.2, 1.2)
ax.set_ylim(-1.2, 1.2)
ax.grid(True, linestyle='--', alpha=0.5)
ax.set_title(f'$p = {p_init:.2f}$')

# 3. 创建滑块（范围从 0.1 到 5.0）
ax_slider = plt.axes([0.2, 0.05, 0.6, 0.03])  # [左, 下, 宽, 高]
slider = Slider(
    ax=ax_slider,
    label='$p$ 值',
    valmin=0.1,
    valmax=2.0,
    valinit=p_init,
    valstep=0.01
)

# 4. 定义滑块更新函数
def update(val):
    p = slider.val
    # 重新计算 Z
    Z = np.abs(X)**p + np.abs(Y)**p

    # 清空当前坐标轴，重新绘制（因为 contour 无法直接更新数据）
    ax.clear()
    ax.contour(X, Y, Z, levels=[1], colors='blue', linewidths=2)
    ax.set_aspect('equal')
    ax.set_xlim(-1.2, 1.2)
    ax.set_ylim(-1.2, 1.2)
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.set_title(f'$p = {p:.2f}$')

    # 刷新画布
    fig.canvas.draw_idle()

# 绑定事件
slider.on_changed(update)

plt.show()
