# Day 11 学习笔记：刚体动力学与 Newton–Euler 方程

> 学习资料：*Modern Robotics: Mechanics, Planning, and Control*（Lynch & Park）  
> 主题：转动惯量、空间惯性矩阵、Body Twist、六维 Newton–Euler 方程及 NumPy 验证  
> 约定：Twist 采用 **angular-first** 顺序 $V=(\omega,v)$，Wrench 采用 $F=(m,f)$。

## 1. 从运动学进入动力学

机械臂动力学通常写作

$$
M(q)\ddot q+C(q,\dot q)\dot q+g(q)=\tau.
$$

本日首先研究**单个刚体**的运动与受力关系，随后才能建立多连杆机器人的逆动力学。

## 2. 转动惯量与角动量

相对于质心 $C$ 的惯性张量定义为

$$
I_C=\int\big(\|r\|^2I_3-rr^T\big)\,dm,
$$

其中 $r$ 是质点相对于质心的位置向量。角动量（以质心为参考点）为

$$
h_b=I_C\omega_b.
$$

与线动量 $p=mv$ 不同，$I_C$ 是矩阵，因此 $h_b$ **不一定**与 $\omega_b$ 平行。

### 2.1 简单实例

- 两个质量均为 $m$ 的质点位于 $(\pm a,0,0)$：

$$
I_C=\operatorname{diag}(0,2ma^2,2ma^2).
$$

- 均匀细杆沿 $x$ 轴、质量 $m$、长度 $L$（忽略截面尺寸）：

$$
I_C=\operatorname{diag}\left(0,\frac{mL^2}{12},\frac{mL^2}{12}\right).
$$

### 2.2 主惯性轴与重复特征值

主惯性方向是 $I_C$ 的特征向量。如果

$$
I_Cu=\lambda u,\qquad \omega=\Omega u,
$$

则 $h=\lambda\omega$，从而 $\omega\times h=0$。

例如

$$
I_C=\operatorname{diag}(2,2,5),\quad \omega=(1,1,0)^T,
$$

有 $h=(2,2,0)^T=2\omega$。由于前两个特征值相同，整个 $xy$ 平面内的**非零**向量都是对应特征值 $2$ 的特征向量，而不只有坐标轴方向。

## 3. Body Twist 与动能

定义 Body Twist：

$$
V_b=\begin{bmatrix}\omega_b\\v_b\end{bmatrix}.
$$

当 Body 坐标系原点位于质心时，$v_b=R^T\dot p_C$，它确实是质心的平动速度在 Body 坐标系中的表达。

注意区分 Space Twist：其线速度分量一般不是末端坐标系原点的实际速度；对于 $T=(R,p)$，有

$$
\dot p=v_s+\omega_s\times p.
$$

在质心 Body 坐标系中，空间惯性矩阵为

$$
G=\begin{bmatrix}I_C&0\\0&mI_3\end{bmatrix},
$$

动能为

$$
K=\frac12V_b^TGV_b
=\frac12\omega_b^TI_C\omega_b+\frac12m\|v_b\|^2.
$$

**计算例子**：$m=2$、$I_C=\operatorname{diag}(1,2,3)$、$V_b=(0,0,2,1,0,0)^T$，则旋转动能为 $6\,\mathrm J$、平动动能为 $1\,\mathrm J$，总动能为 **$7\,\mathrm J$**。注意系数 $1/2$。

### 3.1 原点不在质心时

令 $c$ 表示 Body 原点 $O$ 指向质心 $C$ 的向量，以 Body 坐标表达。则

$$
v_C^b=v_b+\omega_b\times c.
$$

对应的空间惯性矩阵为

$$
G_O=\begin{bmatrix}
I_C-m[c]_\times[c]_\times&m[c]_\times\\
-m[c]_\times&mI_3
\end{bmatrix}.
$$

这里 $I_C$ 仍是**关于质心**的惯性张量。因为 $[c]_\times^T=-[c]_\times$，该矩阵是对称的。其动能仍满足 $K=\tfrac12V_b^TG_OV_b$。

例如 $c=(1,0,0)^T$、$\omega_b=(0,0,2)^T$、$v_b=0$、$m=3$，则 $v_C^b=(0,2,0)^T$，质心平动动能为 $6\,\mathrm J$；**这不包括转动动能**。

## 4. 旋转坐标系中的 Newton–Euler 方程

对质心处的合外力矩和合外力，在 Body 坐标系下有

$$
\boxed{m_b=I_C\dot\omega_b+\omega_b\times(I_C\omega_b)}
$$

$$
\boxed{f_b=m\big(\dot v_b+\omega_b\times v_b\big)}.
$$

这里的点表示 **Body 坐标系中的分量对时间求导**。出现叉乘项是因为 Body 坐标系本身在旋转。

### 4.1 平动项的物理来源

由于质心在惯性系中的速度为 $v_C^s=Rv_b$，且 $\dot R=R[\omega_b]_\times$，

$$
\begin{aligned}
a_C^s
&=\frac{d}{dt}(Rv_b)\\
&=R\big(\dot v_b+\omega_b\times v_b\big).
\end{aligned}
$$

因此 $f_b=R^T(ma_C^s)=m(\dot v_b+\omega_b\times v_b)$。

**例子**：$m=3$、$\omega_b=(0,0,2)^T$、$v_b=(1,0,0)^T$、$\dot v_b=0$，则

$$
f_b=m\omega_b\times v_b=(0,6,0)^T\ \mathrm N.
$$

虽然 Body 速度分量不变，但其惯性系方向不断变化，故仍有加速度。

### 4.2 转动项与陀螺耦合

角动量满足 $h_s=Rh_b$，所以

$$
\dot h_s=R\big(\dot h_b+\omega_b\times h_b\big).
$$

由 $m_s=\dot h_s$、$h_b=I_C\omega_b$，得到转动方程。

当 $I_C=\operatorname{diag}(I_x,I_y,I_z)$ 时，陀螺项为

$$
\omega\times(I_C\omega)=
\begin{bmatrix}
(I_z-I_y)\omega_y\omega_z\\
(I_x-I_z)\omega_x\omega_z\\
(I_y-I_x)\omega_x\omega_y
\end{bmatrix}.
$$

例如 $I_C=\operatorname{diag}(1,2,3)$、$\omega_b=(1,1,0)^T$、$\dot\omega_b=0$：

$$
h_b=(1,2,0)^T,\qquad m_b=(0,0,1)^T\ \mathrm{N\,m}.
$$

原因：$h_b$ 与 $\omega_b$ 不平行；即使 $h_b$ 在 Body 系中的分量恒定，其在惯性系中的方向仍发生变化，需要外力矩。若绕主惯性轴匀速旋转，则 $h_b\parallel\omega_b$，陀螺项为零。注意：无外力矩的自由刚体也可以存在 Body 角速度随时间变化的运动；中间主惯性轴附近的自由旋转具有不稳定性（中间轴定理）。

## 5. 六维 Newton–Euler 方程

定义

$$
V_b=\begin{bmatrix}\omega_b\\v_b\end{bmatrix},\quad
F_b=\begin{bmatrix}m_b\\f_b\end{bmatrix},\quad
\operatorname{ad}_{V_b}=\begin{bmatrix}
[\omega_b]_\times&0\\
[v_b]_\times&[\omega_b]_\times
\end{bmatrix}.
$$

则

$$
\boxed{F_b=G\dot V_b-\operatorname{ad}_{V_b}^TGV_b}.
$$

### 5.1 从六维方程展开到三维方程

对质心坐标系，

$$
G\dot V_b=\begin{bmatrix}I_C\dot\omega_b\\m\dot v_b\end{bmatrix},\qquad
GV_b=\begin{bmatrix}I_C\omega_b\\mv_b\end{bmatrix}.
$$

由于叉乘矩阵是反对称的，

$$
\operatorname{ad}_{V_b}^T=
\begin{bmatrix}
-[\omega_b]_\times&-[v_b]_\times\\
0&-[\omega_b]_\times
\end{bmatrix}.
$$

所以

$$
\begin{aligned}
\operatorname{ad}_{V_b}^TGV_b
&=\begin{bmatrix}
-\omega_b\times(I_C\omega_b)-v_b\times(mv_b)\\
-\omega_b\times(mv_b)
\end{bmatrix}\\
&=\begin{bmatrix}
-\omega_b\times(I_C\omega_b)\\
-m\omega_b\times v_b
\end{bmatrix},
\end{aligned}
$$

其中 $v_b\times(mv_b)=0$。代回六维方程即可得到第 4 节的两个三维方程。

### 5.2 为什么不能删去 $[v]_\times$？

在**质心坐标系 + 标准刚体惯性矩阵**的这一特殊计算中，$[v]_\times(mv)=0$，因此它对最终结果没有贡献；但这并不意味着 $[v]_\times$ 可以从 $\operatorname{ad}_V$ 的定义中删去。

1. **非质心坐标系**：此时 $GV$ 的线动量部分是 $m(v+\omega\times c)$，对应的叉乘一般为 $v\times m(\omega\times c)\neq0$。
2. **李代数结构**：对于 $V_1=(\omega_1,v_1)$、$V_2=(\omega_2,v_2)$，

$$
[V_1,V_2]=\operatorname{ad}_{V_1}V_2
=\begin{bmatrix}
\omega_1\times\omega_2\\
\omega_1\times v_2+v_1\times\omega_2
\end{bmatrix}.
$$

删去左下角的 $[v_1]_\times$ 将破坏正确的 $se(3)$ 李括号。

另外，**$\operatorname{Ad}_T$ 与 $\operatorname{ad}_V$ 不同**：前者负责坐标系之间的 Twist 变换，后者是李代数中的无穷小伴随算子。它们满足（$[V_s]=\dot TT^{-1}$）：

$$
\frac{d}{dt}\operatorname{Ad}_{T(t)}
=\operatorname{ad}_{V_s}\operatorname{Ad}_{T(t)}.
$$

## 6. Day 11 NumPy 编程练习

实现了四个核心函数：

| 函数                          | 功能                                   |
| --------------------------- | ------------------------------------ |
| `vec_to_so3(w)`             | 计算 $[w]_\times$                      |
| `spatial_inertia(m, I_C)`   | 构造质心处的 $6\times6$ 空间惯性矩阵             |
| `ad_twist(V)`               | 构造 $\operatorname{ad}_V$             |
| `newton_euler(G, V, V_dot)` | 计算 $G\dot V-\operatorname{ad}_V^TGV$ |

关键实现：

```python
import numpy as np


def spatial_inertia(m, I_C):
    I_C = np.asarray(I_C, dtype=float).reshape(3, 3)
    G = np.zeros((6, 6))
    G[:3, :3] = I_C
    G[3:, 3:] = m * np.eye(3)
    return G


def ad_twist(V):
    V = np.asarray(V, dtype=float).reshape(6)
    w, v = V[:3], V[3:]
    Ad = np.zeros((6, 6))
    Ad[:3, :3] = vec_to_so3(w)
    Ad[3:, :3] = vec_to_so3(v)
    Ad[3:, 3:] = vec_to_so3(w)
    return Ad


def newton_euler(G, V, V_dot):
    G = np.asarray(G, dtype=float).reshape(6, 6)
    V = np.asarray(V, dtype=float).reshape(6)
    V_dot = np.asarray(V_dot, dtype=float).reshape(6)
    return G @ V_dot - ad_twist(V).T @ G @ V
```

测试参数：

$$
m=3,\quad I_C=\operatorname{diag}(1,2,3),\quad
V=(1,1,0,1,0,0)^T,\quad \dot V=0.
$$

用户实际运行输出：

```text
6D Newton-Euler result: [ 0.  0.  1.  0.  0. -3.]
3D equations result:    [ 0.  0.  1.  0.  0. -3.]
All tests passed!
```

**结论**：所写的六维动力学函数与独立计算的质心坐标系三维方程一致。该测试验证了本组数据下的实现；更完整的验证还可以加入随机状态、非零加速度及非质心坐标系案例。

## 7. 本日核心总结

1. 刚体的角动量 $h=I_C\omega$ 不一定与角速度平行。
2. Body 坐标系旋转时，向量的惯性系导数包含 $\omega\times(\cdot)$ 运输项。
3. 六维 Newton–Euler 方程统一了平动与转动动力学。
4. 空间惯性矩阵的分块对角形式依赖于**原点位于质心**。
5. $[v]_\times$ 在特定质心计算中消失，但对一般空间动力学与 $se(3)$ 结构不可缺少。
6. 已完成 NumPy 实现和六维／三维交叉验证。

**后续衔接**：从单刚体 Newton–Euler 方程推广到多连杆机械臂的逆动力学与递推 Newton–Euler 算法。
