# 第一章 最优化简介

## 1.1 最优化问题一般描述

$$
\begin{aligned}
    \min_x \quad & f(x) \\
    \text{s.t.} \quad & x \in \chi
\end{aligned}
$$

$x$ 称决策变量，$f$ 称目标函数，$\chi$ 称可行域。

## 1.2 稀疏优化实例

问题描述：

$$
\begin{aligned}
    \min_{x\in \mathbb{R}^n} \quad & ||x||_0 \\
    \text{s.t.} \quad & Ax=b
\end{aligned}
$$

简化：1、替换目标函数 $\ell_0$ 范数为 $\ell_1$ 范数（如下)：

$$
\begin{aligned}
    \min_{x\in \mathbb{R}^n} \quad & ||x||_1 \\
    \text{s.t.} \quad & Ax=b
\end{aligned}
$$

注： $\ell_2$范数等值线光滑，其解不满足稀疏性。

2、通过转化简化1的优化问题为无约束优化问题：

$$
\min_{x\in \mathbb{R}^n} \quad \mu||x||_1 + \frac{1}{2}||Ax - b||_2^2
$$

> 代码见：`sparse_opt.py`

其中最小 $\ell _2$ 范数的求解正确性如下：

对 $x^*=A^T(AA^T)b$ ，首先有 $Ax^*=AA^T(AA^T)^{-1}b=b$ , $x^*$ 确实为方程解；又对任意 $x$ 满足 $Ax=b$ ，考虑 $v=x-x^*$ ，有 $(x^*, v)=(A^T(AA^T)^{-1}b, v)=((AA^T)^{-1}b, Av)=0$ ；最后考察 $||x||_2=||x^*-v||_2=||x^*||_2+||v||_2\geq ||x^*||_2$ ，这一意味着 $x^*$ 为唯一最小 $\ell_2$ 范数解。

## 1.3 低秩矩阵恢复实例

对于矩阵 $M$ ，该矩阵除下标集 $\Omega$ 外的元素都是未知的，且矩阵由于统计原因是低秩的，要恢复该矩阵，则问题描述为：

$$
\begin{aligned}
    \min_{X\in \mathbb{R}^{m\times n}} \quad & \text{rank}(X), \\
    \text{s.t.} \quad & X_{ij}=M_{ij},\quad (i,j)\in\Omega.
\end{aligned}
$$

同1.2类似，有两个简化：

1、用核范数取代秩：

$$
\begin{aligned}
    \min_{X\in \mathbb{R}^{m\times n}} \quad & ||X||_*, \\
    \text{s.t.} \quad & X_{ij}=M_{ij},\quad (i,j)\in\Omega.
\end{aligned}
$$

其中 $||X||_*=\sum_i \sigma_i(X)$ 。

2、 二次罚函数形式：

$$
\min_{X\in \mathbb{R}^{m\times n}} \quad \mu||X||_*+\frac{1}{2}\sum_{(i,j)\in\Omega}(X_{ij}-M_{ij})^2.
$$

## 1.4 深度学习实例

### 1.4.1 MLP

实际上是一个大型的非凸优化问题：

$$
\min_x \quad \sum_{i=1}^m ||h(a;x)-b_i||^2+\lambda r(x)
$$

其中 $a$ 为输入参数， $x$ 为参数， $b$ 为目标， $r$ 为正则项。

### 1.4.2 CNN

实际上在 MLP 的前面加了处理数据的卷基层和下采样层，本质仍是优化问题。

## 1.5 最优化基本概念

求解最优化问题分为三步：构造模型 $\rightarrow$ 确定类型 $\rightarrow$ 求解问题。

### 1.5.1 问题分类

I. 根据决策变量在可行域中的情况分为连续/离散优化问题；

II. 根据是否有约束条件分为无约束优化/约束优化问题，约束优化问题可以转化为无约束优化问题；

III. 根据目标或约束函数中是否设计随机变量分为随机/确定性优化问题，随机优化问题可能带来计算的减少；

IV. 根据目标函数和约束函数是否线性分为线性/非线性规划问题；

V. 根据目标函数和可行域的凸性分为凸/非凸优化问题，凸优化由于最优解为定为全局而带来算法设计和理论分析的便利，非凸优化可能可以转化为凸优化问题；

**Def 1.1（最优解）**  对可行点 $\overline{x}\in \chi$ ，

（1）若 $f(\overline{x})\leq f(x),\quad \forall x\in\chi$ ，则 $\overline{x}$ 为全局最优解；

（2）若存在 $\overline{x}$ 的邻域 $N_\epsilon(\overline{x})$ 使得 $f(\overline{x})\leq f(x),\quad \forall x\in N_\epsilon(\overline{x})$ ，则称 $\overline{x}$ 为局部最优解，特别地，当不等号在出该点外严格小时称为严格局部最优解。

### 1.5.2 优化算法

- 显式解：可用代数表达式表示
- 迭代：无法显示求解
- 常用简化技巧：Taylor展开，对偶问题，拆分变量，块坐标下降
- 收敛速度：Q-收敛速度，R-收敛速度
