# 第二章习题

## 0.1

设 $f: \mathbb{R}^n \to \mathbb{R}$ 是连续函数。证明：  
(1) $f$ 是下半连续的；  
(2) $f$ 的所有 $\alpha$-下水平集是闭集；    
(3) 举一个反例：存在下半连续但不连续的函数

证：(1)由连续定义：$\lim_{x\to x_0}f(x)=f(x_0),\,\forall x_0\in\mathbf{dom}f$, 从而有 $\lim_{x\to x_0}f(x)\geq f(x_0),\,\forall x_0\in\mathbf{dom}f$ 成立，因此为下半连续的；（2）由于 $f$ 是下半连续的，从而其 $\alpha$ - 下水平集为闭集；（3）反例：
$$
f(x)=\begin{cases}
x-1,\quad & x\leq 0\\
x, \quad & x>0
\end{cases}
$$

## 0.2

证明：椭球 $\mathcal{E} = \{x \in \mathbb{R}^n \mid (x - x_c)^T P^{-1}(x - x_c) \leqslant 1\}$（其中 $P \succ 0$）是凸集。

证：任意 $x_1,\,x_2\in \mathcal{E},\, \theta\in[0,1]$ ，考虑 $x=\theta x_1+(1-\theta)x_2$ ，同时记 $x_i-x_c=y_i,\,i=1,2$ ，有：
$$
\begin{aligned}
&(x - x_c)^T P^{-1}(x - x_c)\\
=& \left[\theta y_1+(1-\theta) y_2\right]^TP^{-1}\left[\theta y_1+(1-\theta) y_2\right]\\
=& \lVert \theta y_1+(1-\theta)y_2 \rVert_{P^{-1}}^2\\
\leq & (\lVert \theta y_1 \rVert_{P^{-1}}+\lVert (1-\theta)y_2 \rVert_{P^{-1}})^2\\
=& (\theta + (1-\theta))^2\\
=& 1.
\end{aligned}
$$
故 $x\in \mathcal{E}$ ，从而椭球是凸集。
