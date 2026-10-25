# ROBO_Learn · Day 13 学习笔记

## 主题：从 RNEA 到关节空间动力学与正动力学仿真

**教材对应**：*Modern Robotics*, Chapter 8（开链机器人动力学）  
**本日完成状态**：理论理解、NumPy 实现与自动测试全部通过。  
**范围说明**：本日以已经实现的 Newton–Euler 递推（RNEA）为基础，数值提取惯性矩阵、速度相关项，构建正动力学并完成积分。**Euler–Lagrange 方程的系统推导留到 Day 14**；本文只列出作为独立数值验证的二连杆解析公式，不要求本日掌握其拉格朗日推导。

---

## 1. 关节空间动力学模型

对具有 $n$ 个关节的串联机械臂：

$$
\boxed{M(q)\ddot q+c(q,\dot q)+g(q)=\tau}
$$

- $q,\dot q,\ddot q\in\mathbb R^n$：关节位置、速度和加速度。
- $\tau\in\mathbb R^n$：关节驱动力矩（对移动关节则是驱动力）。
- $M(q)\in\mathbb R^{n\times n}$：**关节空间惯性矩阵**，仅依赖构型 $q$。
- $c(q,\dot q)\in\mathbb R^n$：科里奥利、离心效应构成的**速度相关项向量**。
- $g(q)\in\mathbb R^n$：重力项。

有时将速度项写作 $C(q,\dot q)\dot q$。注意 $C$ 的矩阵表示可能不唯一，而 $c=C\dot q$ 作为给定动力学模型中的速度项向量是明确的。

**本日程序的假设**：无重力、无摩擦、无外部负载。因此

$$
\boxed{M(q)\ddot q+c(q,\dot q)=\tau.}
$$

### 逆动力学与正动力学

- **逆动力学**：已知 $(q,\dot q,\ddot q)$，求实现该运动的 $\tau$。
- **正动力学**：已知 $(q,\dot q,\tau)$，求系统实际产生的 $\ddot q$。

$$
\underbrace{\tau=\operatorname{RNEA}(q,\dot q,\ddot q)}_{\text{逆动力学}},\qquad
\underbrace{\ddot q=M(q)^{-1}[\tau-c(q,\dot q)-g(q)]}_{\text{正动力学}}.
$$

正动力学在代码中应优先用 `np.linalg.solve(M, tau - c - g)`，而非显式计算 `np.linalg.inv(M)`。

---

## 2. 为什么 RNEA 可以提取惯性矩阵？

对于固定 $q$ 和 $\dot q$，关节力矩对加速度是仿射的：

$$
\tau=M(q)\ddot q+h(q,\dot q),\qquad h=c+g.
$$

令 $\dot q=0$ 且本模型无重力，则 $c(q,0)=0$，有

$$
\operatorname{RNEA}(q,0,\ddot q)=M(q)\ddot q.
$$

记 $e_k$ 为第 $k$ 个标准基向量。取 $\ddot q=e_k$，得到

$$
\boxed{M(q)[:,k]=\operatorname{RNEA}(q,0,e_k).}
$$

因此，逐列调用 $n$ 次 RNEA 即可构造完整惯性矩阵：

```python
def mass_matrix_from_rnea(robot, q):
    n = robot.n
    M = np.zeros((n, n))
    qd = np.zeros(n)
    for k in range(n):
        qdd = np.zeros(n)
        qdd[k] = 1.0
        M[:, k] = rnea(robot, q, qd, qdd)
    return M
```

**物理理解**：第 $k$ 列描述“只令第 $k$ 个关节产生单位加速度、其他关节加速度为零”时，所有关节需要提供的惯性力矩。非对角元素 $M_{ik}$ 表示关节之间的**惯性耦合**。

例如

$$
M=\begin{bmatrix}3&1\\1&2\end{bmatrix},\qquad
\ddot q=\begin{bmatrix}0\\1\end{bmatrix}
\quad\Longrightarrow\quad
\tau=\begin{bmatrix}1\\2\end{bmatrix}.
$$

虽然 $\ddot q_1=0$，关节 1 仍需施加力矩 $\tau_1=1$，以抵消关节 2 加速带来的惯性耦合。

---

## 3. 从 RNEA 提取速度相关项

在无重力模型中，令 $\ddot q=0$：

$$
\boxed{c(q,\dot q)=\operatorname{RNEA}(q,\dot q,0).}
$$

```python
def velocity_bias_from_rnea(robot, q, qd):
    return rnea(robot, q, qd, np.zeros(robot.n))
```

这解释了 Day 12 中的观察：**关节加速度为零，并不意味着驱动力矩为零**；只要关节速度非零，可能仍存在速度耦合力矩。

在当前理想刚体模型中，速度项关于 $\dot q$ 为二次齐次：

$$
\boxed{c(q,\alpha\dot q)=\alpha^2c(q,\dot q).}
$$

这是一个很方便的数值测试：将关节速度全部加倍，速度项应变为原来的四倍。

**注意重力扩展**：如果将来在 RNEA 中加入重力，则 `rnea(q, qd, 0)` 返回 $c+g$，而不是单独的 $c$；此时 $g(q)=\operatorname{RNEA}(q,0,0)$，并可计算 $c=\operatorname{RNEA}(q,\dot q,0)-g(q)$。

---

## 4. 正动力学：从力矩计算加速度

已知 $M(q)$、$c(q,\dot q)$ 和 $\tau$，解线性方程组：

$$
\boxed{M(q)\ddot q=\tau-c(q,\dot q).}
$$

```python
def forward_dynamics(robot, q, qd, tau):
    M = mass_matrix_from_rnea(robot, q)
    c = velocity_bias_from_rnea(robot, q, qd)
    return np.linalg.solve(M, tau - c)
```

对于上节示例，若改为给定 $\tau=(0,2)^T$，则

$$
\ddot q=M^{-1}\tau
=\frac15\begin{bmatrix}2&-1\\-1&3\end{bmatrix}
\begin{bmatrix}0\\2\end{bmatrix}
=\begin{bmatrix}-0.4\\1.2\end{bmatrix}.
$$

**关节 1 的驱动力矩为零，不代表其加速度为零。** 多关节系统必须整体求解动力学方程，不能把每个关节当作独立的 $I\ddot q=\tau$。

---

## 5. 惯性矩阵的基本性质

机械臂总动能可写为

$$
\boxed{K(q,\dot q)=\tfrac12\dot q^TM(q)\dot q.}
$$

若第 $i$ 根连杆的 Body Twist 为 $V_i=J_{b,i}(q)\dot q$，其空间惯性为 $G_i$，则

$$
K=\frac12\sum_i V_i^TG_iV_i
=\frac12\dot q^T\left(\sum_i J_{b,i}^TG_iJ_{b,i}\right)\dot q.
$$

因此

$$
\boxed{M(q)=\sum_i J_{b,i}(q)^TG_iJ_{b,i}(q).}
$$

由于 $G_i$ 对称，可得 $M=M^T$。对于质量、惯量非退化且关节运动能够产生实际刚体运动的机械臂，有

$$
\boxed{x^TM(q)x>0\quad(\forall x\ne0),\qquad M(q)\succ0.}
$$

其物理含义是：非零关节速度对应正的动能。正定性保证惯性矩阵可逆，也保证在给定状态与力矩时正动力学解的唯一性。

**数值检查**：

```python
np.testing.assert_allclose(M, M.T, atol=1e-10)
assert np.min(np.linalg.eigvalsh(M)) > 0
```

---

## 6. 本日二连杆模型与独立解析校验

本日代码 `Planar2R` 使用如下参数：

| 参数 | 连杆 1 | 连杆 2 |
|---|---:|---:|
| 长度 $L_i$ | 2 | 1 |
| 质量 $m_i$ | 2 | 1 |
| 质心距离近端关节 $r_i=L_i/2$ | 1 | 0.5 |
| 绕质心 $z$ 轴惯量 $I_i$ | 0.25 | 0.10 |

在该平面模型中，解析惯性矩阵为

$$
\begin{aligned}
M_{11}&=I_1+I_2+m_1r_1^2+m_2(L_1^2+r_2^2+2L_1r_2\cos q_2),\\
M_{12}=M_{21}&=I_2+m_2(r_2^2+L_1r_2\cos q_2),\\
M_{22}&=I_2+m_2r_2^2.
\end{aligned}
$$

代入数值：

$$
\boxed{M(q)=\begin{bmatrix}
6.60+2\cos q_2&0.35+\cos q_2\\
0.35+\cos q_2&0.35
\end{bmatrix}.}
$$

定义 $h=m_2L_1r_2\sin q_2=\sin q_2$，则速度项为

$$
\boxed{c(q,\dot q)=h
\begin{bmatrix}
-(2\dot q_1\dot q_2+\dot q_2^2)\\
\dot q_1^2
\end{bmatrix}.}
$$

这两组解析表达式由刚体动能与 Euler–Lagrange 方法得到，本日仅用作**独立参考答案**。推导中的相对角、绝对角速度、质心速度以及拉格朗日求导，将在 **Day 14** 正式学习。

注意 $M$ 只依赖 $q_2$：因为整个机械臂在水平面内同时旋转 $q_1$，不改变各连杆间的相对构型和动能的系数。这是该模型的旋转对称性体现。

### 本次实际测试结果

在

$$q=(0.4,0.7)^T,\qquad \dot q=(1.2,-0.8)^T$$

处，程序输出：

```text
M(q) =
[[8.12968437 1.11484219]
 [1.11484219 0.35      ]]
c(q,qd) = [0.82459864 0.92767347]
eigenvalues(M) = [0.1933937  8.28629067]
```

两个特征值均为正，与惯性矩阵正定性一致。

---

## 7. 从正动力学到状态空间模型

定义状态

$$
x=\begin{bmatrix}q\\\dot q\end{bmatrix}.
$$

则无重力机械臂可写为

$$
\boxed{\dot x=\begin{bmatrix}
\dot q\\
M(q)^{-1}[\tau-c(q,\dot q)]
\end{bmatrix}.}
$$

这一步把机器人动力学转换成控制理论中熟悉的非线性状态空间系统，为后续 PD、Computed Torque、MPC 及强化学习仿真提供基础。

### 半隐式 Euler 积分

对时间步长 $\Delta t$：

$$
\begin{aligned}
\ddot q_k&=\operatorname{ForwardDynamics}(q_k,\dot q_k,\tau_k),\\
\dot q_{k+1}&=\dot q_k+\Delta t\,\ddot q_k,\\
q_{k+1}&=q_k+\Delta t\,\dot q_{k+1}.
\end{aligned}
$$

注意位置更新使用的是**新速度** $\dot q_{k+1}$，这是半隐式 Euler 与普通显式 Euler 的区别。

```python
def simulate(robot, q0, qd0, torque_fn, dt=0.002, steps=1000):
    times = np.linspace(0, steps * dt, steps + 1)
    qs = np.zeros((steps + 1, robot.n))
    qds = np.zeros((steps + 1, robot.n))
    qs[0] = q0
    qds[0] = qd0

    for k in range(steps):
        tau = torque_fn(times[k], qs[k], qds[k])
        qdd = forward_dynamics(robot, qs[k], qds[k], tau)
        qds[k + 1] = qds[k] + dt * qdd
        qs[k + 1] = qs[k] + dt * qds[k + 1]
    return times, qs, qds
```

### 自由运动中的能量检验

在无重力、无摩擦、零输入 $\tau=0$ 的理想系统中，连续时间的动能应守恒：

$$
K(t)=\tfrac12\dot q(t)^TM(q(t))\dot q(t)=\text{常数}.
$$

但数值积分会带来离散化误差，因此数值动能可能有小幅漂移。本日结果：

```text
Free-motion kinetic energy:
initial=3.551171, final=3.551086
relative drift=2.393e-05
```

这是较小的数值漂移，不代表连续时间动力学不守恒。减小积分步长或使用更合适的积分方法可以进一步研究数值误差。

---

## 8. 工程中是否需要每次用 RNEA 构造 M？

机器人运行时 $q,\dot q$ 不断变化，动力学计算确实需要不断更新；但不同任务对应不同算法：

| 算法 | 计算目标 | 典型串联链复杂度 | 使用场景 |
|---|---|---|---|
| **RNEA** | $(q,\dot q,\ddot q)\to\tau$ | $O(n)$ | 逆动力学、力矩前馈 |
| **CRBA** | $q\to M(q)$ | $O(n^2)$ | 直接计算惯性矩阵 |
| **ABA** | $(q,\dot q,\tau)\to\ddot q$ | $O(n)$ | 高效正动力学 |

- **逆动力学控制**：通常一次 RNEA 就能求出力矩，**不需要先显式构造 $M,c,g$**。
- **本日正动力学练习**：通过 $n$ 次 RNEA 取得 $M$，再用一次 RNEA 取得偏置项，最后线性求解；适合教学与小型机械臂。
- **大型机器人仿真**：可采用 CRBA、ABA 等专用算法；含接触的仿真还需处理接触力、摩擦和约束，不能只依靠无约束正动力学。

---

## 9. 验收记录与易错点

### 已完成的程序任务

1. `mass_matrix_from_rnea`：逐列提取 $M(q)$。
2. `velocity_bias_from_rnea`：零加速度输入提取 $c(q,\dot q)$。
3. `forward_dynamics`：使用 `np.linalg.solve` 求解 $\ddot q$。
4. `simulate`：半隐式 Euler 积分，并保存全部时间步状态。

### 自动测试结果

- **100 组随机构型**：惯性矩阵与解析公式一致，且对称正定。
- **100 组随机状态**：速度项与解析公式一致，且满足速度二次齐次性。
- **100 组随机状态**：正动力学与逆动力学往返验证通过。
- **静止零输入**：机械臂保持静止。
- **自由运动**：动能仅存在很小的数值漂移。

### 常见误区

1. **RNEA 只能求力矩，不能用于构造 $M$？** 可以利用单位加速度逐列提取。
2. **$\ddot q=0$ 就意味着 $\tau=0$？** 不对，速度项与重力项仍可能非零。
3. **某关节 $\tau_i=0$ 就意味着 $\ddot q_i=0$？** 不对，惯性矩阵的非对角项体现耦合。
4. **必须显式求 $M^{-1}$？** 不必，使用线性方程求解更合适。
5. **仿真动能略有漂移就意味着模型错误？** 不一定，需要区分动力学模型误差和积分误差。
6. **`rnea(q, qd, 0)` 总等于 $c$？** 仅在无重力等相应假设下成立；一般返回偏置项 $c+g$。

---

## 10. Day 13 → Day 14 衔接

**Day 13 已完成**：从 RNEA 提取 $M,c$；正动力学求解；状态空间表达；半隐式 Euler 仿真；惯性矩阵对称正定与能量检验。

**Day 14 待系统学习**：

1. 广义坐标、广义力与 Euler–Lagrange 方程的意义。
2. 平面 2R 的相对关节角与连杆绝对角速度：$\omega_1=\dot q_1$，$\omega_2=\dot q_1+\dot q_2$。
3. 两根连杆质心位置与速度，以及平动、转动动能。
4. 从 $\mathcal L=K-U$ 独立推导 $M(q),c(q,\dot q),g(q)$。
5. 与 RNEA 数值结果交叉验证，理解两种建模方法的等价性。

**一句话总结**：Day 12 解决“给定运动需要多大力矩”，Day 13 解决“给定力矩将产生怎样的运动”，Day 14 再解释“为什么动力学方程具有这样的结构”。
