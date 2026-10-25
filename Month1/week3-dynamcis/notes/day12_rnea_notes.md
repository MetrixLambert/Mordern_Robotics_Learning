# Day 12｜递归牛顿–欧拉算法（RNEA）与 NumPy 验证

> ROBO_Learn · Week 3 · *Modern Robotics*, Chapter 8  
> 学习状态：**已完成**。4 组典型案例 + 100 组随机案例通过；最大误差约为 $10^{-15}$。  
> 练习代码：`day12_rnea_exercises.py`（仅使用 NumPy）。

## 1. 本日目标与总体思路

**逆动力学问题**：已知关节位置 $q$、速度 $\dot q$、加速度 $\ddot q$ 以及机器人惯性参数，求实现该运动所需的关节驱动力矩 $\tau$。

RNEA 分为两步：

1. **前向递推**（基座 $\rightarrow$ 末端）：求每根连杆的 Body Twist $V_i$ 和其时间导数 $\dot V_i$。
2. **反向递推**（末端 $\rightarrow$ 基座）：求各连杆关节所需传递的 Wrench $F_i$，并投影为关节力矩 $\tau_i$。

RNEA 可以直接计算 $\tau$，不必预先显式构造整个惯性矩阵 $M(q)$ 和速度耦合矩阵 $C(q,\dot q)$。

## 2. 约定：坐标系、Screw Axis 与矩阵

采用 **angular-first** 排列：

$$
V=\begin{bmatrix}\omega\\v\end{bmatrix},\qquad
F=\begin{bmatrix}m\\f\end{bmatrix}\in\mathbb R^6.
$$

- $T_{i-1,i}={}^{i-1}T_i$：子坐标系 $\{i\}$ 相对于父坐标系 $\{i-1\}$ 的位姿，将子系坐标变换到父系。
- $V_i$：连杆 $i$ 的 **Body Twist**，在自身坐标系 $\{i\}$ 中表达。
- $A_i$：**关节 $i$ 的相对运动 Screw Axis，在连杆 $i$ 自身的 Body 坐标系 $\{i\}$ 中表达**；对于本练习固定轴的单自由度关节，$A_i$ 是常量。
- $G_i$：连杆 $i$ 在坐标系 $\{i\}$ 中的 $6\times6$ 空间惯性矩阵。
- $F_i$：父连杆通过关节 $i$ 施加给连杆 $i$ 的 Wrench，表达于 $\{i\}$。

关节相对运动满足

$$
[T_{i-1,i}^{-1}\dot T_{i-1,i}]=[A_i]\dot q_i.
$$

这里的 $A_i$ **不等同于**整条机械臂 PoE 表达中统一在基座零位形表达的 $S_i$，也通常不等同于统一在末端零位形表达的 $B_i$。

对 $T=\begin{bmatrix}R&p\\0&1\end{bmatrix}$，

$$
\operatorname{Ad}_T=\begin{bmatrix}R&0\\[p]_\times R&R\end{bmatrix},\qquad
\operatorname{ad}_V=\begin{bmatrix}[\omega]_\times&0\\[v]_\times&[\omega]_\times\end{bmatrix}.
$$

特别注意 $\operatorname{Ad}_{T^{-1}}=(\operatorname{Ad}_T)^{-1}$。

## 3. 前向速度递推

从位姿连乘 $T_{s,i}=T_{s,i-1}T_{i-1,i}$ 出发，对时间求导，并使用 Body Twist 定义 $[V_i]=T_{s,i}^{-1}\dot T_{s,i}$，可得

$$
\boxed{V_i=\operatorname{Ad}_{T_{i-1,i}^{-1}}V_{i-1}+A_i\dot q_i.}
$$

**物理意义**：第一项是父连杆的运动传递到子连杆并改用子系表达；第二项是当前关节产生的相对运动。两项均在 $\{i\}$ 中表达，因此可以相加。

对于固定基座，$V_0=0$。

## 4. 前向加速度递推

对速度递推求导。由于相邻连杆坐标变换随关节运动而变化，不能把 Adjoint 当成常矩阵：

$$
\frac{d}{dt}\left(\operatorname{Ad}_{T_{i-1,i}^{-1}}V_{i-1}\right)
=\operatorname{Ad}_{T_{i-1,i}^{-1}}\dot V_{i-1}
+\operatorname{ad}_{V_i}(A_i\dot q_i).
$$

于是

$$
\boxed{\dot V_i=\operatorname{Ad}_{T_{i-1,i}^{-1}}\dot V_{i-1}
+A_i\ddot q_i+\operatorname{ad}_{V_i}(A_i\dot q_i).}
$$

第三项是由坐标变换随时间变化产生的**速度耦合项**，即使 $\ddot q_i=0$，也可能非零。

### 为什么第 1 根连杆在匀速时 $\dot V_1=0$？

固定基座时，$V_1=A_1\dot q_1$。因为 Lie bracket 反对称，

$$
\operatorname{ad}_V V=[V,V]=0,
$$

所以

$$
\operatorname{ad}_{V_1}(A_1\dot q_1)
=\operatorname{ad}_{V_1}V_1=0,
\qquad\boxed{\dot V_1=A_1\ddot q_1.}
$$

若 $\ddot q_1=0$，则 $\dot V_1=0$。但**这不表示连杆上各点的物理加速度都为零**：匀速转动的点仍有向心加速度。Body Twist 的导数与某个点在惯性系中的线加速度不是同一个量。

对于第二根连杆，$V_2$ 还包含从连杆 1 传递来的 Twist，因此 $\operatorname{ad}_{V_2}(A_2\dot q_2)$ 一般不为零。

## 5. 单刚体 Newton–Euler 方程

在连杆自身的 Body 坐标系中，净外 Wrench 满足

$$
\boxed{D_i=G_i\dot V_i-\operatorname{ad}_{V_i}^{T}G_iV_i.}
$$

其中 $D_i$ 是连杆 $i$ 所需的**净外 Wrench**，不等于仅由父关节提供的 $F_i$：它还包含子连杆的反作用。

注意：**这里是 $\operatorname{ad}_{V_i}^{T}$，不是 $\operatorname{ad}_{V_i}$**。二者通常不相同。

### 空间惯性矩阵

设连杆坐标系原点到质心的向量为 $c$，质量为 $m$，质心转动惯量为 $I_C$，令 $C=[c]_\times$，则

$$
G=\begin{bmatrix}
I_C-mCC&mC\\
-mC&mI_3
\end{bmatrix}.
$$

在本练习中，连杆坐标系设在连杆末端，**不在质心**，因此必须使用包含 $c$ 的完整空间惯性矩阵。

## 6. 反向 Wrench 递推：重点理解方向与转置

对连杆 $i$，父关节传递的 Wrench 为 $F_i$，而子连杆通过关节 $i+1$ 施加给它的是 $-F_{i+1}$（先按牛顿第三定律取反，再统一表达坐标系）。

子连杆的 Wrench 表达在 $\{i+1\}$，必须先变换到 $\{i\}$。由于功率对偶性，若

$$
V_{i+1}=\operatorname{Ad}_{T_{i,i+1}^{-1}}V_i,
$$

则对应的反向 Wrench 传递为

$$
F_{i+1}^{\{i\}}=\operatorname{Ad}_{T_{i,i+1}^{-1}}^TF_{i+1}^{\{i+1\}}.
$$

将作用力平衡写成 $D_i=F_i-F_{i+1}^{\{i\}}$，得到

$$
\boxed{
F_i=G_i\dot V_i-\operatorname{ad}_{V_i}^TG_iV_i
+\operatorname{Ad}_{T_{i,i+1}^{-1}}^TF_{i+1}.
}
$$

对于无末端外载荷的情况，递推初值为 $F_{n+1}=0$。

**物理意义**：当前关节必须承担当前连杆的惯性需求，外加所有下游连杆向上传递的载荷。

### 关节力矩投影

关节 $i$ 的相对 Twist 为 $A_i\dot q_i$。由功率关系

$$
F_i^T(A_i\dot q_i)=\tau_i\dot q_i
$$

得到

$$
\boxed{\tau_i=A_i^TF_i.}
$$

它与静力学中的 $\tau=J^TF$ 同源，都是 Twist–Wrench 功率对偶关系。

在 NumPy 中，若 `A`、`F` 均为形状 `(6,)` 的一维数组，`A @ F` 已经是内积，等价于数学中的 $A^TF$；一维数组的 `.T` 不改变形状。

## 7. NumPy 实现核心

### 前向递推

```python
T = robot.relative_transforms(q)
V_prev = np.zeros(6)
Vdot_prev = np.zeros(6)
Vs, Vdots = [], []

for i in range(2):
    H = adjoint(inv_transform(T[i]))
    A = robot.A[i]
    V_i = H @ V_prev + A * qd[i]
    Vdot_i = H @ Vdot_prev + A * qdd[i] + ad(V_i) @ (A * qd[i])
    Vs.append(V_i)
    Vdots.append(Vdot_i)
    V_prev, Vdot_prev = V_i, Vdot_i
```

### 反向递推

```python
T = robot.relative_transforms(q)
tau = np.zeros(2)
F_child = np.zeros(6)

for i in range(1, -1, -1):
    F_i = robot.G[i] @ Vdots[i] - ad(Vs[i]).T @ robot.G[i] @ Vs[i]
    if i < 1:
        F_i += adjoint(inv_transform(T[i + 1])).T @ F_child
    tau[i] = robot.A[i] @ F_i
    F_child = F_i
```

`rnea(robot, q, qd, qdd)` 只需先调用前向递推，再将结果传给反向递推。

### 本次实际排查的错误

**主要错误是相邻变换的索引**：在反向递推 `i=0`（连杆 1）时，要把连杆 2 的 Wrench 从 $\{2\}$ 变换到 $\{1\}$，必须使用 `T[1] = T_{1,2}`，而不是 `T[0] = T_{0,1}`。

正确写法：

```python
adjoint(inv_transform(T[i + 1])).T @ F_child
```

而不是：

```python
adjoint(inv_transform(T[i])).T @ F_child
```

另一个检查点是惯性项必须写 `ad(V).T`。本次用户说明原代码中曾有转置，后来被自动补全删除；最终代码已恢复。

**排错原则**：每次写 Adjoint 前，先在纸上标清该 Wrench 当前属于哪个坐标系、目标是哪个坐标系，再决定变换矩阵及其索引。

## 8. 二连杆解析动力学与数值验收

本练习参数：

| 参数                      | 数值            |
| ----------------------- | -------------:|
| 连杆长度 $L_1,L_2$          | $2,1$         |
| 质量 $m_1,m_2$            | $2,1$         |
| 质心绕 $z$ 轴转动惯量 $I_1,I_2$ | $0.25,0.10$   |
| 质心到近端关节距离 $r_1,r_2$     | $L_1/2,L_2/2$ |

忽略重力时，解析动力学为

$$
\tau=M(q)\ddot q+h(q,\dot q),
$$

其中令 $c_2=\cos q_2,\ s_2=\sin q_2$，

$$
\begin{aligned}
M_{11}&=I_1+I_2+m_1r_1^2+m_2(L_1^2+r_2^2+2L_1r_2c_2),\\
M_{12}=M_{21}&=I_2+m_2(r_2^2+L_1r_2c_2),\\
M_{22}&=I_2+m_2r_2^2,
\end{aligned}
$$

$$
h(q,\dot q)=\begin{bmatrix}
-m_2L_1r_2s_2(2\dot q_1\dot q_2+\dot q_2^2)\\
\phantom{-}m_2L_1r_2s_2\dot q_1^2
\end{bmatrix}.
$$

实际运行结果（与独立解析公式对照）：

| 工况          | RNEA $[\tau_1,\tau_2]$      | 最大绝对误差                |
| ----------- | --------------------------- | ---------------------:|
| 静止，仅关节 2 加速 | $[1.11484219,\ 0.35000000]$ | $2.220\times10^{-16}$ |
| 两关节匀速       | $[0.82459864,\ 0.92767347]$ | $6.661\times10^{-16}$ |
| 混合运动        | $[3.13836114,\ 0.76848967]$ | $1.332\times10^{-15}$ |
| 完全静止        | $[0,\ 0]$                   | $0$                   |

**验收结果：4 组典型案例 + 100 组随机案例全部通过。**

## 9. 两个重要的物理结论

**（1）只有关节 2 加速，关节 1 的力矩也可能非零。**

当 $\dot q=0$、$\ddot q=(0,1)^T$ 时，$\dot V_1=0$，但连杆 2 的惯性载荷仍通过反向递推传到关节 1。对应整体动力学中的 $\tau_1=M_{12}(q)$，反映**惯性耦合**。

**（2）所有关节匀速，也可能需要非零驱动力矩。**

当 $\ddot q=0$ 但 $\dot q\ne0$ 时，前向递推的速度耦合项和单刚体 Newton–Euler 方程的 $-\operatorname{ad}_V^TGV$ 项仍可能非零。对应整体动力学中的 $C(q,\dot q)\dot q$，反映**离心与科里奥利效应**。

## 10. 本日易错点与自检清单

- [x] $A_i$ 在**第 $i$ 根连杆自身的 Body 坐标系**中定义，不是统一的基座系。
- [x] 速度传递：$\operatorname{Ad}_{T_{i-1,i}^{-1}}$。
- [x] 加速度传递：不要漏掉 $\operatorname{ad}_{V_i}(A_i\dot q_i)$。
- [x] 单刚体惯性 Wrench：$G_i\dot V_i-\operatorname{ad}_{V_i}^TG_iV_i$，`ad` **必须转置**。
- [x] 子连杆 Wrench 向父连杆传递：$\operatorname{Ad}_{T_{i,i+1}^{-1}}^T$，注意**索引 $i+1$**。
- [x] 关节投影：$\tau_i=A_i^TF_i$；NumPy 一维向量用 `A @ F`。
- [x] $\dot V_i=0$ 不意味着刚体所有点的物理加速度为零。
- [x] 解析公式与随机数值测试均通过。

## 11. 后续衔接

当前实现不含重力、末端外部 Wrench 或关节摩擦。下一步可通过调整基座加速度初始化引入重力，并研究如何调用 RNEA 计算 $M(q)$ 的各列，从而衔接机器人正动力学与控制算法。
