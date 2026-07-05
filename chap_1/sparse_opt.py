import numpy as np
import cvxpy as cp

m, n = 128, 256
A = np.random.randn(m, n)
u = np.random.randn(n) * (np.random.rand(n) < 0.1)
b = A @ u

# 如果直接用 l2 最小范数解：
x_l2 = A.T @ np.linalg.solve(A @ A.T, b)  # Ax=b 的最小 l2 范数解
print(f"l2 解的稀疏度: {np.sum(np.abs(x_l2) > 1e-6)} / {n}")
print(f"l2 解的误差: {np.linalg.norm(x_l2 - u) / np.linalg.norm(u):.2e}")

# 定义变量
x_l1 = cp.Variable(n)

# 目标：最小化 l1 范数
objective = cp.Minimize(cp.norm(x_l1, 1))

# 约束：Ax = b
constraints = [A @ x_l1 == b]

# 构建并求解
problem = cp.Problem(objective, constraints)
problem.solve(solver=cp.ECOS)  # ECOS 是专门处理 l1 问题的求解器

# 提取 NumPy 数组
x_l1_value = x_l1.value

print(f"l1 解的稀疏度: {np.sum(np.abs(x_l1_value) > 1e-6)} / {n}")
print(f"l1 解的误差: {np.linalg.norm(x_l1_value - u) / np.linalg.norm(u):.2e}")
