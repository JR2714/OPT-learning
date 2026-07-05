# 第二章 基础知识


## 2.1 范数

### 向量范数

**Def.** 2.1（范数）向量空间 $\mathbb{R}^n$ 到 $\mathbb{R}$ 上的非负函数 $||\cdot ||$ 为范数，若它满足：

- 正定性：$\forall v\in \mathbb{R}^n$ ，有 $||v||\geq 0$ ，且 $||v||=0$ ***iff*** $v=0$ ；
- 齐次性：$\forall v\in \mathbb{R}^n, \forall \alpha\in \mathbb{R}$ ，有 $||\alpha v||=|\alpha|||v||$；
- 三角不等式：$\forall v,w \in \mathbb{R}^n$ ，有 $||v+w||\leq ||v||+||w||$.

### 常用范数

$\ell_p$ 范数( $p\geq 1$ ):
$$
||v||_p=(\sum_{i=1}^n|v_i|^p)^{\frac{1}{p}}
$$

特别地，当 $p=1$，单位球为"菱形"，解容易被"吸到"坐标轴，产生稀疏性；$p=2$ 时，单位球为"圆形"，旋转不变。

正定阵 $A$ 诱导的范数：$||x||_A=\sqrt{x^TAx}$.

### Cauchy 不等式

**Prop.** 2.1（Cauchy 不等式）对 $a,b\in\mathbb{R}^n$，有
$$
|a^Tb|\leq ||a||_2||b||_2
$$
且等号成立 ***iff*** $a,b$ 线性相关。

---

### 矩阵范数

同样地将矩阵空间视作定义域到 $\mathbb{R}$ 上的非负函数，若其满足正定性、齐次性、三角不等式则它可以被视为一个矩阵范数。

### 常用范数

对矩阵 $A\in\mathbb{R}^m\times\mathbb{R}^n$，
1. $\ell_1$ 范数：
$$
||A||_1=\sum_{i=1}^m\sum_{j=1}^n|a_{ij}|.
$$

2. $\ell_2$/Frobenius 范数：
$$
||A||_F=\sqrt{Tr(AA^T)}=\sqrt{\sum_{i,j}a_{i,j}^2}
$$

注：

- F 范数具有正交不变形：对任意正交阵 $U\in\mathbb{R}^{m\times m}, V\in\mathbb{R}^{n\times n}$，$||UAV||_F=||A||_F$。

- 原因：第一个等号由 Frobenius 范数定义即可得到，第二、四个等号由 $U,V$ 为正交阵得，第三个等号由 $Tr$运算的可交换性得到。

3. 算子范数

对于给定的向量范数$||\cdot||_{(m)},||\cdot||_{(n)}$ ，诱导的矩阵范数为：
$$
||A||_{(m,n)}=\max_{x\in\mathbb{R}^n, ||x||_(n)=1}||Ax||_{(m)},
$$
当给定的向量范数为 $\ell_p$ 范数时，诱导的范数为矩阵的 $p$ 范数。

注：

- 矩阵的2范数是该矩阵的最大奇异值。
- 相容性：$$||Ax||_{(m)}\leq ||A||_{(m,n)}||x||_{(n)}$$

证明" $||A||_2=\sigma_{max}(A)$ ：

由 $||A||_2 = \max_{||x||_2=1}||Ax||_2$ ，从几何上看，2范数即是该矩阵对应的线性变换将单位球面拉伸为椭球面的最长轴的一半，奇异值基本符合。
具体来说，考虑 $ ||A||_2^2=\max_{||x||_2=1} x^TA^TAx $ ，此时 $A^TA$ 对称，考虑其 Rayleigh 商 ，因 $||x||=1$ ，故写为 $R(A^TA,x)=x^TA^TAx$ ，此时正好是等式右边，而由 Rayleigh 商性质，$x$ 取 $A^TA$ 最大特征值对应的特征向量时恰为 $A$ 的最大奇异值的平方，从而原等式值即为最大奇异值。


4. 核范数
$$
||A||_*=\sum_{i=1}^r\sigma_i.
$$
其中$\sigma_i$为$A$的所有非零奇异值，$r=rank(A)$。

补充：

- F 范数看所有元素的平方和 → 衡量矩阵"大小”
- 2 范数（算子范数）看最大拉伸倍数 → 衡量矩阵"作用强度"
- 核范数看奇异值之和 → 促进低秩性（类似 ℓ₁ 促进稀疏）

### 内积

对 $A,B\in \mathbb{R}^{m\times n}$ 定义 Frobineius 内积：
$$
<A,B>=Tr(AB^T)=\sum_{i=1}^m\sum_{j=1}^na_{ij}b_{ij}.
$$

### Cauchy不等式

对 $A,B\in \mathbb{R}^{m\times n}$ ,
$$
|<A,B>|\leq ||A||_F||B||_F.
$$

## 2.2 导数

### 梯度

**Def.** 2.2（梯度）对函数 $f: \mathbb{R}^n\to\mathbb{R}$ ，且 $f$ 在 $x$ 的一个邻域中有意义，若存在向量 $g\in\mathbb{R}^n$ 满足：
$$
\lim_{p\to 0}\frac{f(x+p)-f(x)-g^Tp}{||p||}=0,
$$
其中 $||\cdot||$ 为任意向量范数（合理性由有限维空间范数的等价性保证），则称 $f$ 在点 $x$ 处 $Fr\'{e}chet$ 可微，此时 $g$ 称为 $f$ 在 $x$ 处的梯度，记作 $\nabla f(x)$。若 $f$ 在区域 $D$ 中每一点 $\nabla f(x)$ 存在，则称 $f$ 在 $D$ 上可微。

注：

- 令 $p=\varepsilon e_i$，可得 $$ \nabla f(x)=\left[\frac{\partial f(x)}{\partial x_1},\frac{\partial f(x)}{\partial x_2},\cdots,\frac{\partial f(x)}{\partial x_n}\right]^T $$

### Hesse 矩阵

**Def.** 2.3（Hesse 矩阵）若函数 $f(x): \mathbb{R}^n\to\mathbb{R}$ 在点 $x$ 处的二阶偏导 $\frac{\partial^2 f(x)}{\partial x_i\partial x_j},\quad i,j=1,2,\cdots,n$ 都存在，则
$$
\nabla^2 f(x)=
\begin{bmatrix}
\frac{\partial^2 f(x)}{\partial x_1^2} & \frac{\partial^2 f(x)}{\partial x_1 \partial x_2} & \cdots & \frac{\partial^2 f(x)}{\partial x_1 \partial x_n} \\
\frac{\partial^2 f(x)}{\partial x_2 \partial x_1} & \frac{\partial^2 f(x)}{\partial x_2^2} & \cdots & \frac{\partial^2 f(x)}{\partial x_2 \partial x_n} \\
\vdots & \vdots & \ddots & \vdots \\
\frac{\partial^2 f(x)}{\partial x_n \partial x_1} & \frac{\partial^2 f(x)}{\partial x_n \partial x_2} & \cdots & \frac{\partial^2 f(x)}{\partial x_n^2}
\end{bmatrix}
$$
称为 $f$ 在 $x$ 处的 Hesse 矩阵。

### Jacobi 矩阵

$f: \mathbb{R}^n\to \mathbb{R}^m$，定义它的 Jacobi 矩阵 $J(x)\in\mathbb{R}^{m\times n}$ ：
$$
J(x) =
\begin{bmatrix}
\frac{\partial f_1}{\partial x_1} & \frac{\partial f_1}{\partial x_2} & \cdots & \frac{\partial f_1}{\partial x_n} \\
\frac{\partial f_2}{\partial x_1} & \frac{\partial f_2}{\partial x_2} & \cdots & \frac{\partial f_2}{\partial x_n} \\
\vdots & \vdots & \ddots & \vdots \\
\frac{\partial f_m}{\partial x_1} & \frac{\partial f_m}{\partial x_2} & \cdots & \frac{\partial f_m}{\partial x_n}
\end{bmatrix}
$$

### Taylor 展开

**Th.** 2.1 设 $f: \mathbb{R}^n\to\mathbb{R}$ 是连续可微的， $p\in\mathbb{R}^n$ ，则
$$
f(x+p)=f(x)+\nabla f(x+tp)^Tp,
$$
其中 $t\in(0,1)$ ，若 $f$ 二阶连续可微，则：
$$
\nabla f(x+p)=\nabla f(x)+\int_0^1\nabla ^2 f(x+tp)p\, \mathrm{d} t,
$$
$$
f(x+p)=f(x)+\nabla f(x)^Tp+\frac{1}{2}p^T\nabla^2 f(x+tp)p.
$$

### 梯度 Lipschitz 连续

**Def.** 2.4（梯度 Lipschitz 连续）给定可微函数 $f$ ，若存在 $L\ge 0$ ，对 $\forall x,y \in \mathbf{dom}f$ 有：
$$
||\nabla f(x) - \nabla f(y)||\leq L||x-y||,
$$
则 $f$ 为梯度 Lipschitz 连续的，相应地记 $L$ 为 Lipschitz 常数。

**LM.** 2.1 <font color=blue>**（二次上界）**</font>设定义域为凸集的可微函数 $f(x)$ 为梯度 L-Lipschitz 连续的，则函数 $f(x)$ 有二次上界：
$$
f(y)\leq f(x)+\nabla f(x)^T(y-x)+\frac{L}{2}||y-x||^2,\quad \forall x,y\in \mathbf{dom}f.
$$

**Co.** 2.1 设可微函数 $f(x)$ 定义域为 $\mathbb{R}^n$ 且存在一个全局极小点 $x^*$ ，若 $f(x)$ 的梯度 L-Lipschitz 连续，则对任意 $x$ 有：
$$
\frac{1}{2L}||\nabla f(x)||^2 \leq f(x)-f(x^*).
$$

### 矩阵变量函数的导数

1. **F 可微**

对函数 $f: \mathbb{R}^{m\times n}\to \mathbb{R}$ 若 对 $X\in\mathbb{R}^{m\times n}$ ，存在矩阵 $G\in\mathbb{R}^{m\times n}$ 满足：
$$
\lim_{V\to 0}\frac{f(X+V)-f(X)-\left<G,V\right>}{||V||}=0,
$$
其中 $||\cdot||$ 任意矩阵范数，则称矩阵变量函数 $f$ 在 $X$ 处 $Fr\'{e}chet$ 可微， $G$ 为 $f$ 在 $Fr\'{e}chet$ 可微意义下的梯度。可以记为：
$$
\nabla f(X) =
\begin{bmatrix}
\frac{\partial f}{\partial x_{11}} & \frac{\partial f}{\partial x_{12}} & \cdots & \frac{\partial f}{\partial x_{1n}} \\
\frac{\partial f}{\partial x_{21}} & \frac{\partial f}{\partial x_{22}} & \cdots & \frac{\partial f}{\partial x_{2n}} \\
\vdots & \vdots & \ddots & \vdots \\
\frac{\partial f}{\partial x_{m1}} & \frac{\partial f}{\partial x_{m2}} & \cdots & \frac{\partial f}{\partial x_{mn}}
\end{bmatrix}
$$

2. **G 可微**

**Def.** 2.5（$G\hat{a}teaux$ 可微）设 $f(X)$ 为矩阵函数，若存在矩阵 $G\in\mathbb{R}^{m\times n}$ 对任意方向 $V\in\mathbb{R}^{m\times n}$ 满足：
$$
\lim_{t\to 0}\frac{f(X+tV)-f(X)-t\left<G,V\right>}{t}=0,
$$
则称 $f$ 关于 $X$ 是 $G\hat{a}teaux$ 可微的。其中 $G$ 称为 $f$ 在 $X$ 处在 $G\hat{a}teaux$ 可微意义下的梯度。

注：

- 性质：F 可微 $\Rightarrow$ G 可微。
- 由于 G 可微的定义计算上较为简便，计算多为使用 G 可微。

### 自动微分

自动微分使用计算机计算导数，分为前项模式和后向模式；前向模式依次计算每个中间变量的取值及其对父变量的偏导数值，后向模式则先利用前项模式计算各个节点的值，在根据计算图逆向计算函数关于中间变量的偏导数。后向模式较前向模式计算复杂度更低。

