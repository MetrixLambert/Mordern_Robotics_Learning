# ROBO_Learn · Day 15 学习总结

## 动力学综合验收：物理一致性、能量守恒与数值积分

> **学习阶段**：第三周（动力学）总结；下一阶段 Day 16 将进入关节空间控制。  
> **配套代码**：`day15_dynamics_validation_exercises.py`  
> **完成情况**：TODO 1–5 全部完成，五组自动测试全部通过。

---

## 1. 从两条路线得到同一动力学模型

对无摩擦、无其他外部 Wrench 的串联机械臂，关节空间动力学为

$$
\boxed{M(q)\ddot q+C(q,\dot q)\dot q+g(q)=\tau.}
$$

- **Newton–Euler / RNEA**：从每个刚体的空间速度、加速度和 Wrench 出发，前向递推运动学，反向递推关节力矩。
- **Euler–Lagrange**：从动能 \(K\)、势能 \(U\) 和 \(\mathcal L=K-U\) 出发，推导同一个关节空间动力学方程。

Day 14 已通过随机测试验证二者一致。若 RNEA **包含重力**，可从逆动力学中分别提取各项：

$$
\begin{aligned}
g(q)&=\operatorname{RNEA}(q,0,0),\\
M(q)e_i&=\operatorname{RNEA}(q,0,e_i)-g(q),\\
c(q,\dot q)&=\operatorname{RNEA}(q,\dot q,0)-g(q).
\end{aligned}
$$

其中 \(c(q,\dot q)=C(q,\dot q)\dot q\)。**注意**：有重力时不能直接把 \(\operatorname{RNEA}(q,0,e_i)\) 当作 \(M\) 的第 \(i\) 列，必须扣除重力项。

## 2. 二连杆模型回顾

使用之前的竖直平面二连杆：\(L_1=2,L_2=1\)，\(m_1=2,m_2=1\)，\(r_1=1,r_2=0.5\)，\(I_1=0.25,I_2=0.10\)，重力加速度大小 \(g_0=9.81\)。关节角从水平正方向测量，\(q_2\) 是相对角。

$$
M(q)=\begin{bmatrix}
6.60+2\cos q_2&0.35+\cos q_2\\
0.35+\cos q_2&0.35
\end{bmatrix},
$$

$$
C(q,\dot q)=\sin q_2\begin{bmatrix}
-\dot q_2&-(\dot q_1+\dot q_2)\\
\dot q_1&0
\end{bmatrix},
$$

$$
U(q)=g_0\big[4\sin q_1+0.5\sin(q_1+q_2)\big],
$$

$$
g(q)=g_0\begin{bmatrix}
4\cos q_1+0.5\cos(q_1+q_2)\\
0.5\cos(q_1+q_2)
\end{bmatrix}.
$$

正动力学计算为

$$
\boxed{\ddot q=\operatorname{solve}\big(M(q),\tau-C(q,\dot q)\dot q-g(q)\big).}
$$

NumPy 中用 `np.linalg.solve(M, rhs)`，避免显式求逆。

## 3. 静态重力平衡与扰动

静止平衡要求 \(\dot q=\ddot q=0\)，因此

$$
\boxed{\tau=g(q^*).}
$$

在水平伸直姿态 \(q^*=(0,0)^T\)：

$$
g(q^*)=\begin{bmatrix}44.145\\4.905\end{bmatrix}\ \mathrm{N\cdot m}.
$$

施加这一恒定力矩且初始状态恰好为 \((q^*,0)\) 时，理想模型保持静止。

若位置扰动为 \(q=(0.1,0)^T\)，而 \(\tau=g(q^*)\) 不变，则

$$
\tau-g(q)=[1-\cos(0.1)]g(q^*)\ne0.
$$

**重要细节**：净关节力矩的两个分量均为正，**不代表**两个关节的角加速度都为正，因为

$$
\ddot q=M(q)^{-1}[\tau-g(q)].
$$

惯性耦合会影响各关节的实际加速度方向。

此外，**单纯重力补偿不等于位置反馈**：

- `τ = g(q*)`：仅在目标姿态精确抵消重力；位置改变后一般不再平衡。
- `τ = g(q)`：在任意姿态都抵消重力，但没有指定目标位置的恢复机制。
- `τ = g(q) - K_d q̇`：增加速度耗散，仍然没有位置误差反馈，不能保证收敛到指定 \(q_d\)。

这为 Day 16 的 PD 位置控制和后续重力补偿控制提供了动机。

## 4. 机械能守恒与功率平衡

总机械能：

$$
\boxed{E(q,\dot q)=K+U=\tfrac12\dot q^TM(q)\dot q+U(q).}
$$

由 Day 14 的 \(\dot M-2C\) 斜对称性，可得

$$
\boxed{\dot E=\tau^T\dot q.}
$$

积分后：

$$
\boxed{E(T)-E(0)=\int_0^T\tau(t)^T\dot q(t)\,dt.}
$$

- **零输入 \(\tau=0\)**：连续时间系统总机械能守恒。
- **非零输入**：机械能的变化等于关节驱动力矩所做的功。
- **有重力时**：守恒的是 \(K+U\)，不是单独的动能 \(K\)。

**课堂思考题**：从 \(q(0)=(0,0)^T\)、\(\dot q(0)=0\) 自由释放，初始 \(K(0)=U(0)=0\)。下落后 \(K>0\) 不违反能量守恒，因为质心下降使 \(U<0\)，且 \(K=-U\)。势能零点可以任意平移常数，不改变 \(g(q)=\nabla U(q)\)。

## 5. 半隐式 Euler 数值积分

已实现的离散更新：

$$
\begin{aligned}
\ddot q_k&=f(q_k,\dot q_k,\tau_k),\\
\dot q_{k+1}&=\dot q_k+\Delta t\,\ddot q_k,\\
q_{k+1}&=q_k+\Delta t\,\dot q_{k+1}.
\end{aligned}
$$

与显式 Euler 不同，位置更新使用的是**新速度** \(\dot q_{k+1}\)。这种方法是一阶数值方法；在当前非线性机器人模型中，不保证逐步严格守恒能量。

数值仿真还应注意：

1. 不修改传入的 `q0`、`qd0` 数组；
2. 每一步先在 \((t_k,q_k,\dot q_k)\) 计算输入力矩，再更新状态；
3. 时间网格必须与积分步长一致，推荐：

```python
n_steps = round(duration / dt)
assert np.isclose(n_steps * dt, duration)
t = np.arange(n_steps + 1) * dt
```

## 6. 实验结果与分析

### 6.1 五组测试结果

| 测试 | 内容 | 结果 |
|---|---|---|
| 1 | 100 组随机状态的能量表达式 | PASS |
| 2 | 半隐式 Euler 单步更新恒等式 | PASS |
| 3 | 静态重力补偿与扰动初值 | PASS |
| 4 | 自由运动能量漂移与步长比较 | PASS |
| 5 | 有输入力矩时的功—能量平衡 | PASS |

### 6.2 步长与能量漂移

| \(\Delta t\) (s) | \(\max_k|E_k-E_0|\) (J) |
|---:|---:|
| 0.004 | 0.2198 |
| 0.002 | 0.1100 |
| 0.001 | 0.05506 |

当步长减半时，最大能量漂移也近似减半：

$$
\frac{0.2198}{0.1100}\approx2.00,\qquad
\frac{0.1100}{0.05506}\approx2.00.
$$

因此，在**本次仿真时长和步长范围内**，观察到近似一阶收敛：

$$
\max_k|E_k-E_0|\approx O(\Delta t).
$$

这与半隐式 Euler 的一阶精度相符；不能据此推断任意长时间仿真中的能量误差都会保持有界。

### 6.3 输入功与能量残差

程序采用左端点功率积分：

$$
W_{k+1}=W_k+\Delta t\,\tau_k^T\dot q_k,
$$

并定义残差

$$
\boxed{R_k=E_k-E_0-W_k.}
$$

本次终端残差为：

$$
\boxed{R(T)=-5.5456\times10^{-2}\ \mathrm J.}
$$

该残差同时包含**状态积分误差**与**输入功数值求积误差**。有限步长下 \(R(T)\ne0\) 不等于物理模型违反能量平衡。进一步可通过减小 \(\Delta t\) 检查残差是否趋于零。

## 7. 代码实现回顾

本日完成以下函数：

```python
def total_energy(robot, q, qd): ...
def semi_implicit_euler_step(robot, q, qd, tau, dt): ...
def simulate(robot, q0, qd0, torque_fn, dt, duration): ...
def energy_balance_residual(robot, t, q, qd, tau): ...
def compare_time_steps(robot, q0, qd0, dt_values, duration): ...
```

主要实现思路：

- 能量：`K = 0.5 * qd.T @ M @ qd`，再加势能；
- 状态更新：先速度后位置；
- 仿真：预分配 `q`、`qd`、`tau` 数组，记录初值与每一步结果；
- 功—能量残差：累积离散功，并与机械能增量比较；
- 步长实验：对每个 `dt` 独立运行零输入仿真，比较最大能量偏差。

## 8. Day 15 的核心结论

1. **模型一致性**：RNEA 与 Euler–Lagrange 是同一物理系统的两种计算路径；可以互相交叉验证。
2. **静态平衡不等于稳定控制**：重力补偿可抵消重力，但不能替代位置反馈。
3. **能量是重要的物理验收指标**：零输入看总能量守恒，非零输入看功—能量平衡。
4. **连续模型正确不代表离散轨迹精确**：仿真误差受积分方法、步长和求积方式影响。
5. **数值收敛测试不可省略**：减小时间步长并观察误差变化，有助于区分建模错误与离散化误差。

### 下一步：Day 16

开始学习**关节空间轨迹与 PD 控制**：建立期望轨迹 \(q_d(t)\)，引入位置误差和速度误差反馈，观察闭环响应；再逐步过渡到重力补偿 PD 与 Computed Torque Control。
