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

- F 范数具有正交不变性：对任意正交阵 $U\in\mathbb{R}^{m\times m}, V\in\mathbb{R}^{n\times n}$，$||UAV||_F=||A||_F$。

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

> 一个例子：`backprop_manual.py`

## 2.3 广义实值函数

**Def 2.6**（广义实值函数）令 $\overline{\mathbb{R}}\overset{\mathrm{def}}{=}\mathbb{R}\cup\{\pm\infty\}$ 为广义实数空间，则映射 $f: \mathbb{R}^n\to\overline{\mathbb{R}}$ 称为广义实值函数。

### 适当函数

**Def 2.7**（适当函数）广义实值函数 $f$ ，非空集合 $\chi$ 。若存在 $x\in\chi$ 使得 $f(x)<+\infty$ ，且对任意 $x\in \chi$ 有 $f(x) > -\infty$ ，则称函数 $f$ 关于集合 $\chi$ 是适的。

### 闭函数

**Def 2.8**（$\alpha$ - 下水平集）对于广义实值函数 $f: \mathbb{R}^n\to\overline{\mathbb{R}}$，
$$
C_{\alpha} = \{x|f(x)\leq\alpha\}
$$
为 $f$ 的 $\alpha$ - 下水平集。

**Def 2.9**（上方图）对于广义实值函数 $f: \mathbb{R}^n\to\overline{\mathbb{R}}$，
$$
\mathbf{epi}f=\{(x,t)\in\mathbb{R}^{n+1}|f(x)\leq t\}
$$
称为 $f$ 的上方图。

**Def 2.10**（闭函数）设 $f: \mathbb{R}^n\to\overline{\mathbb{R}}$，若 **epi** $f$ 为闭集，则 $f$ 为闭函数。

**Def 2.11**（下半连续函数）设广义实值函数 $f: \mathbb{R}^n\to\overline{\mathbb{R}}$，若对任意 $x\in\mathbb{R}^n$ ，有
$$
\liminf_{y\to x} f(y)\geq f(x),\\
i.e.\, \underline{\lim}_{y\to x} f(y)\geq f(x),
$$
则称 $f(x)$ 为下半连续函数。

**Th 2.2**（闭函数与下半连续函数的等价性）设 $f: \mathbb{R}^n\to\overline{\mathbb{R}}$，则下列命题等价：
1. $f(x)$ 的任意 $\alpha$ - 下水平集都是闭集；
2. $f(x)$ 是下半连续的；
3. $f(x)$ 是闭函数。

注：
闭（下半连续）函数通过**加法、与仿射变换的复合、取上确界**不改变闭（下半连续）的性质。

## 2.4 凸集

### 基本概念

**Def 2.12** 如果过集合 $C$ 中的任意两点的直线都在 $C$ 内，则称 $C$ 为仿射集，即
$$
x_1,x_2\in C \Rightarrow \theta x_1+(1-\theta)x_2\in C,\quad \forall \theta\in\mathbb{R}.
$$

**Def 2.13** 如果连接集合 $C$ 中的任意两点的线段都在 $C$ 内，则称 $C$ 为凸集，即
$$
x_1,x_2\in C \Rightarrow \theta x_1+(1-\theta)x_2\in C,\quad \forall \theta\in[0,1].
$$

- 凸组合：$x=\sum \theta_kx_k, \sum \theta_k=1,\theta_k\ge 0.$
- 凸包：集合 $S$ 的所有点的凸组合构成的点集称为 $S$ 的凸包，记作 $\mathbf{conv}S$.

**Def 2.14**（仿射包）设 $S\subset \mathbb{R}^n$ ，称如下集合为 $S$ 的仿射包：
$$
\{x|x=\sum\theta_kx_k,x_k\in S,\sum \theta_k=1\}
$$
记为 $\mathbf{affine}\,S$.

- 锥组合：$x=\theta_1x_1+\theta_2x_2,\,\theta_1>0,\theta_2>0$.
- 凸锥：集合 $S$ 的任意点的锥组合在 $S$ 中，则称 $S$ 为凸锥。

### 2.4.2 重要的凸集

**1. 超平面与半空间**

对任意非零向量 $a$ ，集合 $\{x|a^Tx=b\}$ 称为超平面，集合 $\{x|a^Tx\leq b\}$ 称为半空间，此时 $a$ 是法向量。

**2. 球、椭球、锥**

- 欧几里得球：
$$
B(x_c,r)=\{x|\lVert x-x_c \rVert_2 \leq r\}=\{x_c+ru|\lVert u \rVert_2\leq 1\}
$$
- 椭球：
$$
\{x|(x-x_c)^TP^{-1}(x-x_c)\leq 1\}
$$
其中 $P\in \mathcal{S}_{++}^n$，另一种表示为：
$$
\{x_c+Au|\lVert u \rVert\leq 1\}
$$
其中 $A$ 为非奇异方阵。

- 范数球：
$$
\{x|\lVert x-x_c\rVert\leq r\}
$$

- 范数锥：
$$
\{(x,t)|\lVert x\rVert\leq t\}
$$
其中当范数为欧几里得范数时为二次锥。

**3. 多面体**

集合 $\{x|Ax\leq b,Cx=d\}$ 称为多面体，其中 $A\in\mathbb{R}^{m\times n}$ ，$C\in\mathbb{R}^{p\times n}$。

- 例：一个三角形
$$
A=\begin{bmatrix}
-1 & 0 \\
0 & -1 \\
1 & 1 \\
\end{bmatrix}
,
b=\begin{bmatrix}
0\\
0\\
1\\
\end{bmatrix}
,C=0,d=0.
$$

4. （半）正定锥

$\mathcal{S}^n$ 为 $n\times n$ 对称矩阵集合， $\mathcal{S}_{+}^n=\{X\in\mathcal{S}|X\succeq 0\}$ 为半正定矩阵集合，$\mathcal{S}_+^n$ 为凸锥，又称为半正定锥。

### 2.4.3 保凸运算

证明某集合为凸集又两种方法：

1. 通过定义
$$
x_1,x_2\in C,0\leq\theta\leq 1\Rightarrow\theta x_1+(1-\theta)x_2\in C.
$$
2. 通过证明该集合由凸集经过保凸运算得来。
常见的保凸运算有：取交集，仿射变换。

**Th 2.3** 任意多个凸集的交集是凸集。

**Th 2.4** 设 $f:\mathbb{R}^n\to\mathbb{R}^m$ 为仿射变换 $i.e$ $f(x)=Ax+b, A\in\mathbb{R}^{m\times n}, b\in \mathbb{R}^m$ ，则

1. 凸集在 $f$ 下的像为凸集；
2. 凸集在 $f$ 下的原像为凸集。

### 2.4.4 分离超平面定理

**Th 2.5**（分离超平面定理）若 $C$ 和 $D$ 是两个不交的凸集，则存在非零向量 $a$ 和常数 $b$ 使得
$$
a^Tx\leq b,\forall x\in C , \quad \text{且}\quad a^Tx\geq b, \forall x\in D,
$$
即超平面 $\{x|a^Tx=b\}$ 将 $C,\,D$ 分离了。

> 注：这是 Hahn-Banach 定理的几何形式。

**Th 2.6**（严格分离定理）设 $C$ 是闭凸集，点 $x_0\notin C$ ，则存在非零向量 $a$ 和常数 $b$ ，使得
$$
a^Tx<b,\,\forall x\in C\quad \text{且}\quad a^Tx_0>b.
$$

**Def 2.15**（支撑超平面）给定集合 $C$ 及其边界上一点 $x_0$ ，若 $a\neq 0$ 满足 $a^Tx\leq a^Tx_0$ ， $\forall x\in C$ ，则称集合
$$
\{x|a^Tx=a^Tx_0\}
$$
为 $C$ 在边界点 $x_0$ 处的支撑超平面。

**Th 2.7（支撑超平面定理）** 若 $C$ 是凸集，则在 $C$ 的任意边界点处都存在支撑超平面。

## 2.5 凸函数

### 2.5.1 定义

**Def 2.16 (凸函数)** 设函数 $f$ 为适当函数，若 $\mathbf{dom}\,f$ 为凸集，且
$$
f(\theta x+(1-\theta)y)\leq \theta f(x)+(1-\theta)f(y)
$$
对所有 $x,y\in\mathbf{dom}\,f,\,0\leq\theta\leq1$ 都成立，则称 $f$ 为凸函数。

> 注： 
> 1. 亦可等价为上方图为凸集。   
> 2. 不等号严格成立时称为严格凸函数。

**Def 2.17（强凸函数）** 若存在常数 $m>0$ ，使得
$$
g(x)=f(x)-\frac{m}{2}\lVert x\rVert^2
$$
为凸函数，则称 $f(x)$ 为强凸函数，其中 $m$ 为强凸参数，也称$f(x)$ 为 $m-$ 强凸函数。

这与下面的定义等价：

存在常数 $m>0$ 使得对任意的 $x,y\in \mathbf{dom}f$ 以及 $\theta\in (0,1)$ 有
$$
f(\theta x+(1-\theta)y)\leq \theta f(x)+(1-\theta)f(y)-\frac{m}{2}\theta (1-\theta)\lVert x-y\rVert^2,
$$

**Prop 2.3** 设 $f$ 为强凸函数且**存在最小值**，则 $f$ 的最小值点唯一。

### 2.5.2 凸函数判定定理

**Th 2.8** $f(x)$ 是凸函数 $\iff$ 对任意 $x\in\mathbf{dom}f,\,v\in\mathbb{R}^n,\,g:\mathbb{R}\to\mathbb{R}$ 满足
$$
g(t)=f(x+t\nu),\quad \mathbf{dom}g=\{t|x+t\nu\in\mathbf{dom}f\}
$$
是凸函数。

对于可微函数，还可以通过导数信息判断凸性。

**Th 2.9 （一阶条件）** 定义在凸集上的可微函数 $f$ ，其是凸函数当且仅当
$$
f(y)\geq f(x)+\nabla f(x)^T(y-x),\quad x,y\in\mathbf{dom} f. 
$$
> 注：这实际上描述了曲线（面）在切线（面）上方。
> 利用这一点可以得到函数的全局下界。

**Th 2.10（梯度单调性）** $f$ 可微，则 $f$ 为凸函数当且仅当 $\mathbf{dom}f$ 为凸集且 $\nabla f$ 为单调映射，即
$$
(\nabla f(x)-\nabla f(y))^T(x-y)\geq 0,\quad \forall x,y\in\mathbf{dom}f. 
$$

**Co 2.2** $f$ 可微， $\mathbf{dom}f$ 为凸集，则

(1) $f$ 严格凸当且仅当
$$
(\nabla f(x)-\nabla f(y))^T(x-y)>0,\quad x,y\in\mathbf{dom}f;
$$

(2) $f$ 为 m-强凸函数当且仅当
$$
(\nabla f(x)-\nabla f(y))^T(x-y)\geq m\lVert x-y \rVert^2,\quad x,y\in\mathbf{dom}f. 
$$

若函数二阶可微，则有

**Th 2.11**（二阶条件）设 $f$ 为定义在凸集上的二阶连续可微函数，则 $f$ 为凸函数当且仅当
$$
\nabla^2 f(x)\succeq 0 ,\quad \forall x\in\mathbf{dom}f. 
$$
若 $\nabla^2 f(x)\succ 0$ ，则称函数严格凸。

还可以使用上方图来判定凸性

**Th 2.12** 函数 $f(x)$ 为凸函数当且仅当其上方图 $\mathbf{epi} f$ 为凸集。

**Proof:**

必要性：

$\quad$ $\forall x,y\in\mathbb{R}^n,\,a,b\in\mathbb{R}\quad s.t.\quad f(x)\leq a,\,f(y)\leq b$ ，考虑 $\forall t\in(0,1)$ ，有
$$
f(tx+(1-t)y)\leq tf(x)+(1-t)f(y)\leq ta+(1-t)b
$$
故 $(tx+(1-t)y, ta+(1-t)b)\in\mathbf{epi}f$ ，从而 $f$ 的上方图为凸集。

充分性：

$\quad$ $\forall x,y\in\mathbf{dom}f$ ，$(x,f(x)),\,(y,f(y))\in\mathbf{epi}f$ ，则 $\forall t\in [0,1]$ ，有
$$
f(tx+(1-t)y)\leq tf(x)+(1-t)f(y)
$$
于是 $f$ 为凸函数。

### 2.5.3 保凸运算

保凸运算也能帮助判断一个函数是否是凸函数

**Th 2.13** 

1. $f$ 为凸函数，则 $\alpha f$ 为凸函数， 其中 $\alpha\geq 0$.
2. 若 $f_1,\,f_2$ 为凸函数，则 $f_1+f_2$ 为凸函数.
3. 若 $f$ 为 凸函数，则 $f(Ax+b)$ 为凸函数.
4. 若 $f_1,f_2,\dots , f_m$ 为凸函数，则 $f(x)=\text{max}\{f_1(x),f_2(x),\dots,f_m(x)\}$ 为凸函数.
5. 若对每个 $y\in\mathcal{A}$ ， $f(x,y)$ 为关于 $x$ 的凸函数，则
$$
g(x)=\sup_{y\in \mathcal{A}}f(x,y)
$$
$\qquad$ 是凸函数.

6. 给定函数 $g:\mathbb{R}^n\to\mathbb{R},\,h:\mathbb{R}\to\mathbb{R}$ ，令 $f(x)=h(g(x))$ 。若 $g$ 为凸函数，且 $h$ 为凸函数且单调不减，则 $f$ 是凸函数； $g$ 为凹函数， $h$ 为凸函数且单调不增，则 $f$ 是凸函数。
7. 给定 $g:\mathbb{R}^n\to\mathbb{R}^k,\,h:\mathbb{R}^k\to\mathbb{R}$，
$$
f(x)=h(g(x))=h(g_1(x),g_2(x),\cdots,g_k(x))
$$
$\qquad$ 若 $g_i$ 为凸函数，且 $h$ 凸且关于每个变量单调不减，则 $f$ 为凸函数；若 $g_i$ 为凹函数，且 $h$ 凸且关于每个变量单调不增，则 $f$ 为凸函数。

8. 若 $f(x,y)$ 关于 $(x,y)$ 整体式凸函数， $C$ 是凸集，则 
$$
g(x)=\inf_{y\in C}f(x,y)
$$
$\qquad$ 为凸函数。

9. 函数 $f:\mathbb{R}^n\to\mathbb{R}$ 的**透视函数** $g:\mathbb{R}^n\times\mathbb{R}\to\mathbb{R}$
$$ 
g(x,t)=tf(\frac{x}{t}),\,\mathbf{dom}g=\{(x,t)|\frac{x}{t}\in\mathbf{dom}f,t>0\}
$$
$\qquad$ 当 $f$ 是凸函数时 $g$ 为凸函数。

### 2.5.4 凸函数的性质

1. 连续性

**Th 2.14** 设 $f:\mathbb{R}^n\to (-\infty,+\infty]$ 为凸函数，对任意 $x_0\in\mathbf{int\,dom}f$ ，$f$ 在 $x_0$ 处连续，其中 $\mathbf{int\,dom}f$ 表示 $f$ 的定义域的内点。

**Co. 2.3** $f$ 为凸函数且 $\mathbf{dom}f$ 为开集，则 $f$ 在 $\mathbf{dom}f$ 上是连续函数。

2. 凸下水平集

**Prop. 2.4** 设 $f(x)$ 为凸函数，则 $f(x)$ 所有 $\alpha$ - 下水平集 $C_{\alpha}$ 为凸集。

3. 二次下界

**LM. 2.2**（二次下界）设 $f$ 为可微 m -强凸函数，则其满足：
$$
f(y)\geq f(x)+\nabla f(x)^T(y-x)+\frac{m}{2}\lVert y-x\rVert^2, \, \forall x,y\in \mathbf{dom} f.
$$

> 注：可微强凸函数的下水平集有界，这是由于对 $\forall x\in S_\alpha$，以及全局最小值点 $x^*$：
>$$\alpha\geq f(x)\geq f(x^*)+\frac{m}{2}\lVert x-x^*\rVert^2 \\ \lVert x-x^*\rVert \leq \sqrt{\frac{2(\alpha-f(x^*))}{m}}




