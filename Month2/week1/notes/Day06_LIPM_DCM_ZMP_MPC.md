# Day 6｜LIPM、DCM 与 ZMP-MPC

> 学习主线：从线性倒立摆的发散动力学出发，理解 ZMP 如何控制质心运动；使用精确离散化和 CVXPY 实现有限支撑区域下的 MPC，并认识 DCM 可捕获性约束。

## 1. 目标与基本假设

**目标**：掌握 LIPM 的动力学、精确离散化、预测矩阵构造、ZMP-MPC 二次规划，以及终端 DCM 约束与可行性的关系。

一维 LIPM 假设：

- 质心高度恒定为 `h`，质心角动量变化等效项忽略；
- 地面水平，接触力可用 ZMP 位置 `u` 表征；
- 质心水平位置为 `x`，重力加速度为 `g`；
- 当前支撑脚的 ZMP 可行区间固定。

定义自然频率：

$$
\omega=\sqrt{\frac{g}{h}}.
$$

LIPM 方程：

$$
\boxed{\ddot x=\omega^2(x-u)}
$$

**物理解释**：若 `x > u`，质心加速度指向前方，偏离支撑位置的运动可能进一步放大。LIPM 是一个具有不稳定模态的系统，ZMP `u` 是模型中的控制输入；实际机器人需通过接触力/关节力矩实现期望 ZMP。

## 2. 状态空间形式及稳定性

取状态

$$
X=\begin{bmatrix}x\\\dot x\end{bmatrix},\qquad
\dot X=AX+Bu,
$$

其中

$$
A=\begin{bmatrix}0&1\\\omega^2&0\end{bmatrix},\qquad
B=\begin{bmatrix}0\\-\omega^2\end{bmatrix}.
$$

连续系统的特征值为 `+ω` 和 `−ω`，分别对应发散模态与收敛模态。

### DCM（Divergent Component of Motion）

定义

$$
\boxed{\xi=x+\frac{\dot x}{\omega}}.
$$

对时间求导并代入 LIPM：

$$
\boxed{\dot\xi=\omega(\xi-u)}.
$$

当 ZMP 恒定为 `u` 时：

$$
\xi(t)=u+e^{\omega t}[\xi(0)-u].
$$

因此：

- `ξ(0)=u`：发散模态被消除，DCM 保持不变；
- `ξ(0)≠u`：DCM 与 ZMP 的偏差按指数增长；
- 若 `u=ξ`，则 `ẋ=−ω(x−ξ)`，质心渐近靠近 `ξ`，速度渐近趋于零，并非瞬间停止。

### 有限支撑区域与可捕获性

若 ZMP 必须满足

$$
u\in[u_{\min},u_{\max}],
$$

且初始 DCM 位于区间外，则在**不改变支撑区域**的理想 LIPM 模型中，无法仅靠区间内 ZMP 将 DCM 稳定到某个固定位置。

例如：`h=1 m`、`g=9.81 m/s²`、`x₀=0`、`ẋ₀=0.4 m/s`：

$$
\omega\approx3.1321,\qquad
\xi_0=\frac{0.4}{3.1321}\approx0.12771\ \mathrm m.
$$

当 ZMP 允许范围为 `[-0.1, 0.1] m` 时，`ξ₀ > 0.1 m`，因此固定支撑条件下不可捕获。若 `x₀=0`，临界初速度约为

$$
\dot x_{0,\max}=\omega\cdot0.1\approx0.3132\ \mathrm{m/s}.
$$

> 这里的“不可捕获”针对固定支撑、恒高 LIPM 和给定控制范围，不是对真实机器人所有平衡动作的绝对结论；迈步会改变未来的可行支撑区域。

## 3. LIPM 的精确离散化

采用零阶保持（ZOH），采样周期 `Δt`，在每个采样区间内 ZMP 保持常数：

$$
X_{k+1}=A_dX_k+B_du_k.
$$

精确表达式：

$$
\boxed{
A_d=\begin{bmatrix}
\cosh(\omega\Delta t)&\frac{\sinh(\omega\Delta t)}{\omega}\\
\omega\sinh(\omega\Delta t)&\cosh(\omega\Delta t)
\end{bmatrix}
}
$$

$$
\boxed{
B_d=\begin{bmatrix}
1-\cosh(\omega\Delta t)\\
-\omega\sinh(\omega\Delta t)
\end{bmatrix}.
}
$$

取 `h=1 m`、`Δt=0.05 s`，数值近似：

$$
A_d\approx\begin{bmatrix}1.01228758&0.05020463\\0.49250738&1.01228758\end{bmatrix},\qquad
B_d\approx\begin{bmatrix}-0.01228758\\-0.49250738\end{bmatrix}.
$$

`A_d` 特征值约为 `1.16953309` 和 `0.85504208`，对应离散发散模态和收敛模态。

**编程提醒**：矩阵幂用 `np.linalg.matrix_power(Ad, i)`，不能用 `Ad ** i`（后者是 NumPy 的逐元素乘方）。

## 4. MPC 预测矩阵

给定预测步长 `N`，定义

$$
\mathbf X=\begin{bmatrix}X_{k+1}\\X_{k+2}\\\vdots\\X_{k+N}\end{bmatrix},\qquad
\mathbf U=\begin{bmatrix}u_k\\u_{k+1}\\\vdots\\u_{k+N-1}\end{bmatrix}.
$$

则

$$
\boxed{\mathbf X=\bar A X_k+\bar B\mathbf U},
$$

其中

$$
\bar A=\begin{bmatrix}A_d\\A_d^2\\\vdots\\A_d^N\end{bmatrix},
$$

$$
\bar B=\begin{bmatrix}
B_d&0&\cdots&0\\
A_dB_d&B_d&\cdots&0\\
\vdots&\vdots&\ddots&\vdots\\
A_d^{N-1}B_d&A_d^{N-2}B_d&\cdots&B_d
\end{bmatrix}.
$$

代码中 `A_bar` 尺寸为 `(2N, 2)`，`B_bar` 尺寸为 `(2N, N)`。

**实际验证**：逐步递推与批量矩阵预测的最大数值误差约为 `6.66×10⁻¹⁶`，表明预测矩阵构造正确。

## 5. ZMP-MPC 二次规划

典型目标函数：

$$
\boxed{
\min_{\mathbf U}\ \mathbf X^\top\bar Q\mathbf X+\mathbf U^\top\bar R\mathbf U
}
$$

约束：

$$
u_{\min}\le u_j\le u_{\max},\qquad j=0,\ldots,N-1.
$$

可取 `Q=diag(10,1)`、`R=0.01 I`、`N=20`、`Δt=0.05 s`。使用 CVXPY 的 `cp.Variable(N)`、`cp.quad_form(...)` 和 OSQP 求解。

此目标倾向于减小质心位置、速度偏差和 ZMP 控制幅值，但**有限预测时域内的优化可行不代表无限时域稳定**。尤其对不稳定的 LIPM，若没有适当的终端条件，MPC 可能给出每一步都满足 ZMP 约束、但最终发散的轨迹。

## 6. 终端 DCM 约束

从终端状态提取

$$
\xi_N=x_N+\frac{\dot x_N}{\omega}.
$$

为了要求预测末端 DCM 落入当前固定支撑范围，可添加

$$
\boxed{u_{\min}\le\xi_N\le u_{\max}}.
$$

这是关于决策变量 `U` 的**线性不等式**，因为 `X_N` 是 `U` 的仿射函数。它属于简化的终端可捕获性条件。

实验结论：

- `ẋ₀=0.2 m/s`：可以求得满足约束的解；
- `ẋ₀=0.3 m/s`：基线 MPC 与带终端 DCM 约束的结果非常接近，说明该条件在此参数下不一定显著改变最优解；
- `ẋ₀=0.4 m/s`：终端 DCM 约束下不可行，而不带终端条件的有限时域 QP 仍可能返回 `optimal`，但轨迹不具备固定支撑下的可捕获性。

**重要区分**：优化器报告 `optimal` 仅表示给定 QP 被成功求解，不等于物理系统可持续平衡。

## 7. Python 实现要点与踩坑

### 核心函数结构

```python
get_discrete_lipm(omega, dt)
build_prediction_matrices(Ad, Bd, N)
solve_mpc(x_current, A_bar, B_bar,
          Q_bar, R_bar, u_min, u_max,
          omega, terminal_dcm=False)
```

### CVXPY 求解状态

```python
problem.solve(solver=cp.OSQP)
if problem.status not in ('optimal', 'optimal_inaccurate'):
    return None, None, problem.status
return float(U.value[0]), U.value.copy(), problem.status
```

需要保持调用方与函数返回值的顺序一致；此前出现过因误把字符串状态值当成控制序列处理而导致 NumPy 类型错误的情况。

### 约束表达式

CVXPY 不支持 `lower <= U <= upper` 形式的链式比较，应拆分为：

```python
constraints = [U >= lower, U <= upper]
```

## 8. Day 6 核心结论与 Day 7 衔接

1. LIPM 的 DCM 满足 `ξ̇=ω(ξ−u)`，存在发散模态。
2. ZMP 是控制输入，但它只能位于当前接触所允许的区域内。
3. 固定支撑下，初始 DCM 超出支撑范围会带来可捕获性限制。
4. 终端 DCM 约束能够揭示仅依赖有限时域状态代价时被掩盖的可行性问题。
5. **迈步的本质是改变未来的 ZMP 可行区域**，因此下一阶段需要把落脚点也作为优化变量。

**下一课**：Day 7——Capture Point、DCM offset、单步及两步 Footstep MPC。
