# ROBO_Learn · Day 4 学习笔记

**主题：SE(3)、刚体速度场、Spatial/Body Twist、Adjoint 与指数映射**  
**日期：2026-10-12—13**  
**参考：** Lynch & Park, *Modern Robotics: Mechanics, Planning, and Control*, 第 3 章（Rigid-Body Motions and Twists）  
**状态：** 理论学习完成；`se3_hat`、`se3_vee`、`adjoint`、`se3_exp` 及综合一致性测试全部通过。

> **本日主线：** 刚体**位姿** $T$ → 刚体**速度场** → 六维 **Twist** → 运动学方程 $\dot T=[V_s]T=T[V_b]$ → 恒定 Twist 下的指数映射。

---

## 1. SE(3)：描述刚体位姿

定义坐标系 $\{b\}$ 相对于 $\{a\}$ 的齐次变换

$$
T_{ab}=\begin{bmatrix}R_{ab}&p_{ab}\\0&1\end{bmatrix}\in SE(3),
\qquad R_{ab}\in SO(3),\ p_{ab}\in\mathbb R^3.
$$

- $R_{ab}$：$\{b\}$ 各坐标轴在 $\{a\}$ 中的表达。
- $p_{ab}$：$\{b\}$ 原点在 $\{a\}$ 中的位置。
- 对同一个空间点：$x_a=R_{ab}x_b+p_{ab}$，或 $\bar x_a=T_{ab}\bar x_b$，其中 $\bar x=[x^T,1]^T$。

**组合：**

$$
T_{ac}=T_{ab}T_{bc},\quad
R_{ac}=R_{ab}R_{bc},\quad
p_{ac}=R_{ab}p_{bc}+p_{ab}.
$$

**求逆：**

$$
\boxed{T_{ab}^{-1}=T_{ba}=\begin{bmatrix}R_{ab}^T&-R_{ab}^Tp_{ab}\\0&1\end{bmatrix}.}
$$

易错点：逆变换的平移部分通常**不是** $-p_{ab}$；需要先变换坐标表达。

## 2. Twist 的物理来源：经典刚体速度场

不必先从 $T$ 的导数定义 Twist。经典刚体运动学给出，刚体上任意两点 $A,B$ 的位置和速度满足

$$
\boxed{v_B=v_A+\omega\times(p_B-p_A).}
$$

其中 $\omega$ 是刚体角速度（以上量均在同一固定坐标系中表达）。这个公式来自刚体上任意两点间距离不变的约束。

令 $x=p_B$，选一个参考点 $O$，将公式展开：

$$
\dot x=v_O+\omega\times(x-p_O)
=\omega\times x+\underbrace{(v_O-\omega\times p_O)}_{v_s}.
$$

于是定义 **Spatial Twist**（采用本书的角速度在前的约定）：

$$
\boxed{V_s=\begin{bmatrix}\omega_s\\v_s\end{bmatrix}\in\mathbb R^6,\qquad
\dot x=\omega_s\times x+v_s.}
$$

**最重要的区别：** 对于刚体坐标系原点的位置 $p$，

$$
\boxed{\dot p=\omega_s\times p+v_s,\qquad
v_s=\dot p-\omega_s\times p.}
$$

- $\dot p$：**刚体坐标系原点的真实线速度**（Spatial 坐标表达）。
- $v_s$：**整个刚体速度场的常量项**，一般不等于 $\dot p$。
- 在**固定的 Spatial 坐标系**中，$v_s$ 不依赖选取刚体上的哪个点计算；但若改变 Spatial 坐标系，Twist 分量也会变化。
- 当 $\omega_s=0$ 时，所有点速度相同，$v_s=\dot p$，为瞬时纯平移。

### 2.1 用旋转轴求 Twist

若刚体绕经过 $q$ 的**固定轴纯旋转**，轴上的点静止，故

$$
0=\omega_s\times q+v_s
\quad\Rightarrow\quad
\boxed{v_s=-\omega_s\times q.}
$$

**例 1：绕经过空间原点的 z 轴纯旋转。**

$$
\omega_s=(0,0,2)^T,\quad p=(1,0,0)^T.
$$

因为旋转轴经过空间原点，$v_s=0$；但刚体坐标系原点不在轴上，故

$$
\boxed{\dot p=\omega_s\times p=(0,2,0)^T\ne0.}
$$

**例 2：绕经过 $q=(1,0,0)^T$ 的 z 向固定轴纯旋转。**

$$
\omega_s=(0,0,2)^T\quad\Rightarrow\quad
\boxed{v_s=(0,-2,0)^T.}
$$

此时轴上点 $q$ 的速度为零。**纯旋转不意味着 $v_s=0$**：只有旋转轴经过所选 Spatial 原点时才有该简化。

### 2.2 距离为什么不变？

令 $r=p_B-p_A$，则 $\dot r=\omega\times r$，从而

$$
\frac{d}{dt}\|r\|^2=2r^T\dot r=2r^T(\omega\times r)=0.
$$

两点速度可以不同，但相对速度始终垂直于连线，因而距离不变。

## 3. 从六维 Twist 到 $\mathfrak{se}(3)$ 矩阵

对 $\omega=(\omega_x,\omega_y,\omega_z)^T$，定义

$$
[\omega]_\times=\begin{bmatrix}
0&-\omega_z&\omega_y\\
\omega_z&0&-\omega_x\\
-\omega_y&\omega_x&0
\end{bmatrix},\qquad [\omega]_\times x=\omega\times x.
$$

六维 Twist 的 **hat** 形式：

$$
\boxed{[V]=\begin{bmatrix}[\omega]_\times&v\\0&0\end{bmatrix}\in\mathfrak{se}(3).}
$$

**注意区分：** $T\in SE(3)$ 的右下角是 $1$；$[V]\in\mathfrak{se}(3)$ 的右下角是 $0$。hat/vee 只是在六维向量与结构化 $4\times4$ 矩阵之间转换。

## 4. Spatial Twist 与 Body Twist

令 $T=T_{sb}=[R,p;0,1]$，$R$ 将 Body 坐标表达变换到 Spatial 坐标表达。

从 Day 2 的角速度公式出发：

$$
[\omega_s]_\times=\dot RR^T,\qquad
[\omega_b]_\times=R^T\dot R,\qquad
\boxed{\omega_s=R\omega_b.}
$$

对完整位姿矩阵，定义

$$
\boxed{[V_s]=\dot TT^{-1},\qquad [V_b]=T^{-1}\dot T.}
$$

分块相乘：

$$
\dot TT^{-1}=
\begin{bmatrix}\dot RR^T&\dot p-\dot RR^Tp\\0&0\end{bmatrix}
=\begin{bmatrix}[\omega_s]_\times&v_s\\0&0\end{bmatrix},
$$

$$
T^{-1}\dot T=
\begin{bmatrix}R^T\dot R&R^T\dot p\\0&0\end{bmatrix}
=\begin{bmatrix}[\omega_b]_\times&v_b\\0&0\end{bmatrix}.
$$

所以

$$
\boxed{v_s=\dot p-\omega_s\times p,\qquad v_b=R^T\dot p.}
$$

**理解：** $v_b$ 是**Body 原点的真实线速度在 Body 坐标系中的表达**；$v_s$ 是 Spatial 速度场的常量项。二者不是简单的同一向量在不同坐标系下的分量，因此一般 $v_s\ne Rv_b$。

### 4.1 Adjoint 变换

由 $\dot p=Rv_b$ 和 $\omega_s=R\omega_b$，

$$
\begin{aligned}
v_s&=Rv_b-\omega_s\times p\\
&=Rv_b+[p]_\times R\omega_b.
\end{aligned}
$$

合并得

$$
\operatorname{Ad}_T=\begin{bmatrix}R&0\\[p]_\times R&R\end{bmatrix}}
$$

当 $\omega=0$（纯平移）时，额外的叉乘项消失，$v_s=Rv_b$，退化为普通向量坐标变换。

**综合例子：** $R=I$，$p=(1,0,0)^T$，绕 Spatial 原点 z 轴以 $\omega_s=(0,0,2)^T$ 旋转：

$$
\dot p=(0,2,0)^T,\quad v_s=0,\quad
\omega_b=(0,0,2)^T,\quad v_b=(0,2,0)^T.
$$

$$
V_s=(0,0,2,0,0,0)^T,\qquad
V_b=(0,0,2,0,2,0)^T.
$$

## 5. 螺旋运动（Screw Motion）与螺距

对 $\omega\ne0$，刚体的瞬时运动可用螺旋轴上的点 $q$ 和螺距 $h$ 表达：

$$
\boxed{v=-\omega\times q+h\omega,\qquad
h=\frac{\omega^Tv}{\|\omega\|^2}.}
$$

- $h=0$：瞬时纯旋转，轴上点速度为零。
- $h\ne0$：旋转同时沿螺旋轴移动；每转 $1\,\mathrm{rad}$ 沿轴移动 $h$（有符号长度）。
- $\omega=0$：纯平移，以上螺距公式不适用，需单独处理。

例如 $\omega=(0,0,2)^T$，$v=(0,0,1)^T$，则 $h=1/2$，刚体绕 z 轴以 $2\,\mathrm{rad/s}$ 旋转并沿 z 轴以 $1\,\mathrm{m/s}$ 平移。

## 6. 刚体运动学方程与 SE(3) 指数映射

两种等价运动学表达为

$$
\boxed{\dot T=[V_s]T=T[V_b].}
$$

若 **Spatial Twist $V_s$ 恒定**，则

$$
\boxed{T(t)=e^{[V_s]t}T(0).}
$$

注意：$e^{[V_s]t}$ 是**从单位位姿开始累积的相对运动**；一般初始位姿不是 $I$ 时还要右乘 $T(0)$。

当 $T(0)=I$ 时，设 $V_s=(\omega_s,v_s)$，得到

$$
e^{[V_s]t}=\begin{bmatrix}R(t)&p(t)\\0&1\end{bmatrix},
\quad R(t)=e^{[\omega_s]_\times t},
\quad p(t)=\int_0^t e^{[\omega_s]_\times\tau}v_s\,d\tau.
$$

**为什么这个积分成立？** 令 $A=[\omega_s]_\times$、$b=v_s$，则 $\dot p=Ap+b$、$p(0)=0$。对候选解 $p(t)=\int_0^t e^{A\tau}b\,d\tau$ 求导得 $\dot p=e^{At}b$，同时

$$
Ap+b=\int_0^t Ae^{A\tau}b\,d\tau+b
=(e^{At}-I)b+b=e^{At}b.
$$

### 6.1 闭式表达

若 $\omega_s\ne0$，令 $\alpha=\|\omega_s\|$，$K=[\omega_s/\alpha]_\times$，由 Rodrigues 公式：

$$
R(t)=I+\sin(\alpha t)K+(1-\cos(\alpha t))K^2,
$$

$$
\boxed{p(t)=\left[tI+\frac{1-\cos(\alpha t)}{\alpha}K+
\left(t-\frac{\sin(\alpha t)}{\alpha}\right)K^2\right]v_s.}
$$

若 $\omega_s=0$，则 $R(t)=I$、$p(t)=tv_s$。

**掌握重点：** 理解旋转和平移为何会耦合、能按公式实现即可；不必现在背诵整个平移闭式表达。小角速度时直接使用上述分式会有消减误差，数值计算可考虑级数展开。

## 7. NumPy 实现：本日通过的代码

```python
import numpy as np


def skew(w):
    x, y, z = w
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=float)


def vee(W):
    return np.array([W[2, 1], W[0, 2], W[1, 0]], dtype=float)


def so3_exp(phi):
    theta = np.linalg.norm(phi)
    if theta < 1e-12:
        return np.eye(3)
    K = skew(phi / theta)
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)


def se3_hat(V):
    M = np.zeros((4, 4))
    M[:3, :3] = skew(V[:3])
    M[:3, 3] = V[3:]
    return M


def se3_vee(M):
    return np.concatenate([vee(M[:3, :3]), M[:3, 3]])


def adjoint(T):
    R, p = T[:3, :3], T[:3, 3]
    Ad = np.zeros((6, 6))
    Ad[:3, :3] = R
    Ad[3:, 3:] = R
    Ad[3:, :3] = skew(p) @ R
    return Ad


def se3_exp(V, t):
    omega, v = V[:3], V[3:]
    alpha = np.linalg.norm(omega)
    T = np.eye(4)
    if alpha < 1e-12:
        T[:3, 3] = t * v
        return T
    K = skew(omega / alpha)
    T[:3, :3] = so3_exp(omega * t)
    G = (t * np.eye(3)
         + (1 - np.cos(alpha * t)) / alpha * K
         + (t - np.sin(alpha * t) / alpha) * (K @ K))
    T[:3, 3] = G @ v
    return T
```

> 上面是根据本日提交代码整理的等价清洁版：保留原有数学实现，将 `alpha == 0` 改为容差判断，并用 `G` 表示平移系数矩阵。对于极小但非零的角速度，容差分支是近似；若追求高精度，需要小角度级数实现。

### 7.1 已通过的测试场景

1. `se3_hat` / `se3_vee` 往返一致性。

2. `adjoint(I)=I_6`；以及旋转轴不经过 Body 原点的例子。

3. `se3_exp`：纯平移、绕 Spatial 原点纯旋转、绕偏移轴纯旋转、沿轴螺旋运动。

4. **综合一致性：** 对任意合适的 $T\in SE(3)$、$V_b$ 和 $t$，
   
   $$
   \boxed{T e^{[V_b]t}T^{-1}=e^{[\operatorname{Ad}_T V_b]t}.}
   $$
   
   实际测试中 `np.allclose(lhs, rhs, atol=1e-10)` 通过。

## 8. 易混淆点：考前速查

| 容易混淆的说法                            | 正确理解                                            |
| ---------------------------------- | ----------------------------------------------- |
| “$v_s$ 就是刚体原点速度”                   | 一般不对；$\dot p=\omega_s\times p+v_s$。             |
| “纯旋转一定有 $v_s=0$”                   | 不对；旋转轴经过 Spatial 原点时才有 $v_s=0$。                 |
| “$v_b$ 和 $v_s$ 只差一个 $R$”           | 一般不对；$v_s=Rv_b+[p]_\times R\omega_b$。           |
| “角速度恒定，原点速度一定恒定”                   | 不对；$\ddot p=\omega_s\times\dot p$，只有满足相应条件时才恒定。 |
| “刚体各点速度不同，距离会改变”                   | 不对；相对速度与连线垂直，距离不变。                              |
| “$\exp([V]t)$ 的平移部分总是 $tv$”        | 只有纯平移或旋转耦合项消失等特殊情况成立。                           |
| “已知 $\omega,p$ 就能同时求 $\dot p,v_s$” | 不够；还需已知某点速度、固定旋转轴等额外运动条件。                       |

## 9. 自测题（复习时先不看答案）

1. 已知 $R,p$，写出 $T^{-1}$，解释平移项为什么是 $-R^Tp$。
2. 已知 $\omega=(0,0,2)^T$，刚体绕经过 $q=(1,0,0)^T$ 的固定轴纯旋转，求 $v_s$；再求点 $x=(2,0,0)^T$ 的速度。
3. 已知 $R=I$，$p=(1,0,0)^T$，$V_b=(0,0,2,0,2,0)^T$，求 $V_s$。
4. 对 $V=(0,0,0,1,2,3)^T$，写出 $[V]$ 和 $e^{[V]t}$。
5. 为什么 $v_s$ 不随所选刚体参考点改变，却可能随 Spatial 坐标系改变？
6. 解释恒定 Twist 时 $T(t)=e^{[V_s]t}T(0)$ 的物理意义。

<details>
<summary>参考答案（做完再展开）</summary>

1. $T^{-1}=[R^T,-R^Tp;0,1]$；反向平移须先用 $R^T$ 表达到新坐标系。
2. $v_s=(0,-2,0)^T$；$\dot x=\omega\times x+v_s=(0,2,0)^T$。
3. $V_s=(0,0,2,0,0,0)^T$。
4. $[V]=[0_{3\times3},(1,2,3)^T;0,0]$；$e^{[V]t}=[I,t(1,2,3)^T;0,1]$。
5. 对任意刚体点 $x$，$v_s=\dot x-\omega\times x$ 恒定；但坐标系平移会改变位置表达及速度场的常量项。
6. $V_s$ 在固定 Spatial 坐标系下保持恒定，指数映射给出从初始位姿演化的刚体运动。

</details>

## 10. Day 4 结论与下一步

**本日已达到的核心要求：**

- 能使用 SE(3) 描述、组合和求逆刚体位姿。
- 能从经典刚体速度关系理解 Twist 的六维结构，而非机械记忆矩阵定义。
- 能区分 $v_s$、$\dot p$ 和 $v_b$，并正确使用 Adjoint。
- 能理解恒定 Twist 下的指数映射，独立实现并通过数值测试。

**下一步：** 将 SE(3)、Twist 和指数映射用于机器人关节运动，学习关节螺旋轴、正运动学的 Product of Exponentials（PoE）表示，为 Jacobian 做准备。
