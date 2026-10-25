# ROBO_Learn — Day 6 学习笔记

## Velocity Kinematics：Space/Body Jacobian、奇异性与逆速度运动学

**参考教材**：Kevin M. Lynch & Frank C. Park, *Modern Robotics: Mechanics, Planning, and Control*，第 5 章（Jacobian），并联系第 4 章 PoE；伪逆与阻尼最小二乘作为后续数值逆运动学的预备知识。

**本日成果**：推导 Space/Body Jacobian；理解每列的几何意义及与位置 Jacobian 的区别；分析任务空间奇异性、伪逆和阻尼最小二乘；实现两个 Jacobian 函数，通过 2R、RP、随机机器人及有限差分测试。

---

## 1. 约定与速度定义

采用角速度在前、线速度在后的六维 Twist：

$$
V=\begin{bmatrix}\omega\\v\end{bmatrix},\qquad
[V]=\begin{bmatrix}[\omega]_\times&v\\0&0\end{bmatrix}.
$$

末端位姿 $T(q)=T_{sb}(q)=[R(q),p(q);0,1]$，有

$$
[V_s]=\dot T T^{-1},\qquad [V_b]=T^{-1}\dot T,
\qquad V_s=\operatorname{Ad}_{T}V_b.
$$

Adjoint 为

$$
\operatorname{Ad}_{T}=\begin{bmatrix}R&0\\[p]_\times R&R\end{bmatrix}.
$$

**重要区别**：Space Twist 的线速度分量 $v_s$ 一般不是末端原点的实际线速度 $\dot p$；二者满足

$$
\boxed{\dot p=\omega_s\times p+v_s=-[p]_\times\omega_s+v_s.}
$$

## 2. Space Jacobian：从 Space PoE 推导

正运动学：

$$
T(q)=e^{[\mathcal S_1]q_1}\cdots e^{[\mathcal S_n]q_n}M.
$$

定义 $V_s=J_s(q)\dot q$。以两关节为例，令 $A_i=e^{[\mathcal S_i]q_i}$，则 $T=A_1A_2M$，从而

$$
\begin{aligned}
\dot T T^{-1}
&=\dot A_1A_1^{-1}+A_1\dot A_2A_2^{-1}A_1^{-1}\\
&=[\mathcal S_1]\dot q_1+A_1[\mathcal S_2]A_1^{-1}\dot q_2.
\end{aligned}
$$

因此

$$
\boxed{J_{s,i}(q)=\operatorname{Ad}_{e^{[\mathcal S_1]q_1}\cdots e^{[\mathcal S_{i-1}]q_{i-1}}}\mathcal S_i,\quad J_{s,1}=\mathcal S_1.}
$$

**几何意义**：$J_{s,i}$ 是当前构型下，仅第 $i$ 个关节以单位速度运动时末端产生的 **Space Twist**。对于旋转关节，它就是当前关节轴在固定空间系中的 Screw Axis：

$$
J_{s,i}=\begin{bmatrix}\omega_{s,i}\\-\omega_{s,i}\times r_{s,i}\end{bmatrix},
$$

其中 $r_{s,i}$ 为当前关节轴上一点。它只依赖上游关节 $q_1,\ldots,q_{i-1}$。

## 3. Body Jacobian：两种推导

定义 $V_b=J_b(q)\dot q$。由 $V_b=\operatorname{Ad}_{T^{-1}}V_s$ 得到

$$
\boxed{J_b(q)=\operatorname{Ad}_{T(q)^{-1}}J_s(q).}
$$

注意与 Day 5 的零位 Screw Axis 转换区分：

$$
\mathcal B_i=\operatorname{Ad}_{M^{-1}}\mathcal S_i
\quad\text{（固定零位）},
\qquad
J_b(q)=\operatorname{Ad}_{T(q)^{-1}}J_s(q)
\quad\text{（当前构型）}.
$$

从 Body PoE 直接推导：

$$
T(q)=M e^{[\mathcal B_1]q_1}\cdots e^{[\mathcal B_n]q_n},
\qquad [V_b]=T^{-1}\dot T.
$$

两关节时 $T=MA_1A_2$（此处 $A_i=e^{[\mathcal B_i]q_i}$），则

$$
T^{-1}\dot T=A_2^{-1}[\mathcal B_1]A_2\dot q_1+[\mathcal B_2]\dot q_2.
$$

推广可得

$$
\boxed{J_{b,i}(q)=\operatorname{Ad}_{e^{-[\mathcal B_n]q_n}\cdots e^{-[\mathcal B_{i+1}]q_{i+1}}}\mathcal B_i,\quad J_{b,n}=\mathcal B_n.}
$$

**几何意义**：$J_{b,i}$ 是当前构型下，仅第 $i$ 个关节以单位速度运动时末端产生的 **Body Twist**，即当前关节轴的 Screw Axis 用**当前末端坐标系**表示。它只依赖下游关节 $q_{i+1},\ldots,q_n$。

对于旋转关节：

$$
J_{b,i}=\begin{bmatrix}\omega_{b,i}\\-\omega_{b,i}\times r_{b,i}\end{bmatrix},
$$

其中轴方向和轴上一点均在当前末端系中表示。

## 4. 平面 2R：Space 与 Body Jacobian

设 $L_1=2,L_2=1$，零位末端 $M=[I,(3,0,0)^T;0,1]$，有

$$
\begin{aligned}
\mathcal S_1&=(0,0,1,0,0,0)^T, &\mathcal S_2&=(0,0,1,0,-2,0)^T,\\
\mathcal B_1&=(0,0,1,0,3,0)^T, &\mathcal B_2&=(0,0,1,0,1,0)^T.
\end{aligned}
$$

第 2 关节当前轴经过 $r_2=(2\cos q_1,2\sin q_1,0)^T$，所以

$$
\boxed{J_s(q)=\begin{bmatrix}
0&0\\0&0\\1&1\\0&2\sin q_1\\0&-2\cos q_1\\0&0
\end{bmatrix}.}
$$

当 $q_1=\pi/2,q_2=0$：

$$
J_s=\begin{bmatrix}0&0\\0&0\\1&1\\0&2\\0&0\\0&0\end{bmatrix},
\qquad
J_b=\begin{bmatrix}0&0\\0&0\\1&1\\0&0\\3&1\\0&0\end{bmatrix}.
$$

保持 $q_2=0$、改变 $q_1$ 时，$J_b$ 不变。如果 $\dot q$ 也不变，则 **Body Twist 的数值不变**；但由于 $T(q)$ 改变，$V_s=\operatorname{Ad}_{T(q)}V_b$ 通常改变。不能说“$J_b$ 不变但 $V_b$ 必然改变”。

例如 $\dot q=(0,1)^T$ 时，$V_b=\mathcal B_2$ 始终相同，但末端原点的空间速度为

$$
\dot p=(-\sin q_1,\cos q_1,0)^T.
$$

## 5. 位置 Jacobian 与 Twist Jacobian 不同

末端位置

$$
p(q)=\begin{bmatrix}2\cos q_1+\cos(q_1+q_2)\\2\sin q_1+\sin(q_1+q_2)\end{bmatrix}.
$$

对位置求偏导可定义平面位置 Jacobian：

$$
\boxed{J_p(q)=\frac{\partial p}{\partial q}=\begin{bmatrix}
-2\sin q_1-\sin(q_1+q_2)&-\sin(q_1+q_2)\\
2\cos q_1+\cos(q_1+q_2)&\cos(q_1+q_2)
\end{bmatrix}.}
$$

它满足 $\dot p=J_p\dot q$，**不是** $J_s$ 的下三行。一般在三维中，将 $J_s=[J_\omega^T,J_v^T]^T$ 按行分块，则

$$
\boxed{J_p^{(3D)}=-[p]_\times J_\omega+J_v
=\begin{bmatrix}-[p]_\times&I_3\end{bmatrix}J_s.}
$$

常见的几何速度 Jacobian（角速度在前）定义为

$$
\begin{bmatrix}\omega_s\\\dot p\end{bmatrix}
=J_g\dot q,
\qquad
J_g=\begin{bmatrix}I_3&0\\-[p]_\times&I_3\end{bmatrix}J_s.
$$

因此 $J_g$ 与 $J_s$ 秩相同（左乘可逆矩阵），但只取平移速度的 $J_p$ 可以有不同的秩。

## 6. 奇异性：先明确任务空间

由于 $\operatorname{Ad}_{T^{-1}}$ 可逆，始终有

$$
\boxed{\operatorname{rank}J_b(q)=\operatorname{rank}J_s(q).}
$$

对平面 2R 的**位置任务**，

$$
\det J_p(q)=2\sin q_2.
$$

故 $q_2=0$ 或 $\pi$ 时，$J_p$ 秩由 2 降为 1，不能产生任意平面末端原点瞬时速度。但完整的 $6\times2$ $J_s$ 在此例中仍为秩 2：**位置任务的奇异不等于完整 Twist Jacobian 的秩下降**。

特别地，$q_1=q_2=0$ 时

$$
J_p=\begin{bmatrix}0&0\\3&1\end{bmatrix},
\qquad \dot x=0,\quad \dot y=3\dot q_1+\dot q_2.
$$

期望 $\dot p_d=(1,0)^T$ 无法由任何关节速度瞬时实现。此结论是一阶瞬时速度限制，不意味着所有有限位移都不可能；不过完全伸直时 $x=3$ 已达到最大可达 $x$ 坐标，因此向正 $x$ 的有限位移也不可达。

## 7. 伪逆与阻尼最小二乘

逆速度问题：给定 $\dot p_d$，求 $\dot q$ 使 $J\dot q\approx\dot p_d$。

Moore–Penrose 伪逆解：

$$
\boxed{\dot q^*=J^\dagger\dot p_d.}
$$

它先使 $\|J\dot q-\dot p_d\|_2$ 最小，再在所有最小残差解中选择 $\|\dot q\|_2$ 最小的一个；$JJ^\dagger\dot p_d$ 是 $\dot p_d$ 到 $\operatorname{Im}J$ 的正交投影。

对于 $J=\begin{bmatrix}0&0\\3&1\end{bmatrix}$、$\dot p_d=(1,0)^T$，有

$$
J^\dagger=\begin{bmatrix}0&3/10\\0&1/10\end{bmatrix},\qquad
\dot q^*=(0,0)^T.
$$

伪逆无法创造 Jacobian 列空间之外的速度方向。

**接近奇异**时，即使 $J$ 可逆，最小奇异值很小也会使关节速度过大。例如

$$
J=\operatorname{diag}(1,0.001),\quad \dot p_d=(0,1)^T
\ \Longrightarrow\ \dot q=(0,1000)^T,
\quad \kappa_2(J)=1000.
$$

阻尼最小二乘（Damped Least Squares, DLS）通过正则化平衡跟踪误差和关节速度：

$$
\boxed{\dot q=\arg\min_z\big(\|Jz-\dot p_d\|_2^2+\lambda^2\|z\|_2^2\big)
=(J^TJ+\lambda^2I)^{-1}J^T\dot p_d,\quad \lambda>0.}
$$

在上述例子中，$\lambda=0.01$ 时 $\dot q_2=0.001/(0.001^2+0.01^2)\approx9.90$，以跟踪误差换取较小关节速度。实际编程宜用 `np.linalg.solve` 求解线性系统，避免显式计算逆矩阵。

## 8. NumPy 实践：高效递推

输入 $S_{list},B_{list}\in\mathbb R^{6\times n}$，**每列是一个零位 Screw Axis**，$q\in\mathbb R^n$。下述实现复用 Day 5 的 `adjoint()` 和 `se3_exp()`。

```python
import numpy as np

def jacobian_space(S_list, q):
    S_list = np.asarray(S_list, dtype=float)
    q = np.asarray(q, dtype=float).reshape(-1)
    assert S_list.shape == (6, len(q))
    n = len(q)
    J = np.zeros((6, n))
    T = np.eye(4)
    for i in range(n):
        J[:, i] = adjoint(T) @ S_list[:, i]
        if i < n - 1:
            T = T @ se3_exp(S_list[:, i], q[i])
    return J


def jacobian_body(B_list, q):
    B_list = np.asarray(B_list, dtype=float)
    q = np.asarray(q, dtype=float).reshape(-1)
    assert B_list.shape == (6, len(q))
    n = len(q)
    J = np.zeros((6, n))
    T = np.eye(4)
    for i in reversed(range(n)):
        J[:, i] = adjoint(T) @ B_list[:, i]
        if i > 0:
            T = T @ se3_exp(B_list[:, i], -q[i])
    return J
```

**用户原始实现为何也正确？** 原实现计算 $J_{s,i}$ 时多包含了当前关节的指数；计算 $J_{b,i}$ 时也多包含了当前关节的指数后再取逆。这不改变结果，因为

$$
\boxed{\operatorname{Ad}_{e^{[\mathcal S_i]q_i}}\mathcal S_i=\mathcal S_i,\qquad
\operatorname{Ad}_{e^{-[\mathcal B_i]q_i}}\mathcal B_i=\mathcal B_i.}
$$

原因是 $[\mathcal S_i]$ 与其矩阵指数可交换，Body 情况同理。

原双重循环每个 Jacobian 需计算约 $n(n+1)/2$ 次 `se3_exp`，上述递推只需 $n-1$ 次，矩阵指数次数从 $O(n^2)$ 降为 $O(n)$。反向计算 Body Jacobian 时，**必须右乘负指数**以得到 $e^{-[\mathcal B_n]q_n}\cdots e^{-[\mathcal B_{i+1}]q_{i+1}}$ 的正确顺序。

如果同时需要 FK、$J_s$ 和 $J_b$，还可复用同一正向递推并使用 $J_b=\operatorname{Ad}_{T^{-1}}J_s$。

## 9. 测试与验收

在 `day6_jacobian_practice.py` 中实现两种 Jacobian 后，已通过全部测试：

- `test_2r_home`：零位 $J_s(0)=S$、$J_b(0)=B$。
- `test_2r_specific`：$q_1=\pi/2,q_2=0$ 的数值与 $J_b$ 对 $q_1$ 的独立性。
- `test_rp`：RP 机械臂中移动轴随上游旋转关节转动，以及 Body/Space 转换。
- `test_random_consistency`：40 组随机机器人 FK、一致性 $J_b=\operatorname{Ad}_{T^{-1}}J_s$ 和秩检验。
- `test_finite_difference`：Space Twist、Body Twist 与位置 Jacobian 的有限差分验证。

运行：

```bash
python day6_jacobian_practice.py
```

**本日易错点总结**：

1. $J_s$ 的下三行不是一般意义上的末端原点位置 Jacobian。
2. Body Jacobian 使用当前末端系；$B_i$ 是零位定义的常量，不等于任意构型的 $J_{b,i}$。
3. $J_s$ 从上游递推，$J_b$ 从下游递推；矩阵指数的乘法顺序不可随意交换。
4. 讨论“奇异”必须说明任务空间：完整 Twist 与仅末端位置的秩表现可能不同。
5. 伪逆求最小二乘最小范数解，不保证实现不可达速度；DLS 用正则化缓解病态放大。

## 10. 下一步：Day 7

进入第 6 章数值逆运动学：从期望末端位姿 $T_d$ 出发，使用 $SE(3)$ 的位姿误差 Twist、Jacobian 和 Newton–Raphson 迭代求关节变量；重点区分 Space/Body 误差表达与对应 Jacobian 的一致性。
