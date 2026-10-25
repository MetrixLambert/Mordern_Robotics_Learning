# Day 17｜机器人 PD 控制、LaSalle 不变性原理与模型失配

## 学习目标与完成情况

- 利用机器人动力学结构证明 **PD + 精确重力补偿** 的位置镇定渐近稳定性。
- 理解 LaSalle 不变性原理、正不变集、最大不变集之间的区别。
- 区分位置镇定与时变轨迹跟踪，以及动力学补偿的两种方式。
- 用 NumPy 验证 **Computed Torque Control (CTC)** 在重力模型失配下的稳态偏差。
- **结果：4 组测试全部通过。**

## 1. 机器人动力学及结构性质

关节空间动力学为

\[
M(q)\ddot q+C(q,\dot q)\dot q+g(q)=\tau,
\]

其中 \(M(q)=M(q)^T\succ0\)，采用满足标准性质的 \(C\) 定义时，

\[
\dot M(q)-2C(q,\dot q)\quad\text{为斜对称矩阵}.
\]

因此对任意向量 \(z\)，

\[
\frac12 z^T\dot Mz=z^TCz.
\]

这是能量法分析机械臂控制的关键。

## 2. 固定目标：PD + 重力补偿

令 \(q^*\) 为固定目标，\(e=q-q^*\)，采用

\[
\tau=g(q)-K_pe-K_d\dot e,\qquad K_p,K_d\succ0.
\]

得到闭环方程

\[
M(q)\ddot e+C(q,\dot q)\dot e+K_d\dot e+K_pe=0.
\]

取 Lyapunov 函数

\[
V=\frac12\dot e^TM(q)\dot e+\frac12e^TK_pe.
\]

利用斜对称性质得到

\[
\dot V=-\dot e^TK_d\dot e\le0.
\]

这里 \(\dot V\) 仅半负定，不能直接凭负定性定理断言渐近稳定；需要进一步分析 \(\dot V=0\) 的集合。

## 3. LaSalle 不变性原理（拉萨尔不变性原理）

对自治系统 \(\dot x=f(x)\)，若在紧致正不变集 \(\Omega\) 上有 \(\dot V\le0\)，令

\[
E=\{x\in\Omega:\dot V(x)=0\}.
\]

LaSalle 原理指出：轨迹趋近于 \(E\) 中的**最大不变集** \(\mathcal M\)，而非必然趋近整个 \(E\) 中的任意平衡点。

- **正不变集**：从集合内出发，未来轨迹始终留在集合内。
- **\(E\)**：\(\dot V=0\) 的状态集合；它**不一定不变**。
- **最大不变集 \(\mathcal M\)**：包含在 \(E\) 中、按集合包含关系最大的（完整轨迹）不变子集。

### 如何判断不变性

若集合为规则等式约束 \(S=\{x:h(x)=0\}\)，可检查向量场是否处处与集合相切：

\[
Dh(x)f(x)=0,\quad x\in S.
\]

在适当的光滑性、正则性和解唯一性条件下，这可用于验证不变性。对于不等式集合，则检查边界向量场是否指向内部或切向；Lyapunov 子水平集 \(\Omega_c=\{V\le c\}\) 在 \(\dot V\le0\) 下是正不变的。

### 机械臂应用

\[
E=\{(e,\dot e):\dot e=0\}.
\]

在 \(E\) 内，要持续保持 \(\dot e=0\)，还需要 \(\ddot e=0\)。由闭环方程，

\[
M(q)\ddot e+K_pe=0.
\]

故 \(\ddot e=0\) 当且仅当 \(e=0\)，最大不变集为

\[
\mathcal M=\{(e,\dot e)=(0,0)\}.
\]

在理想模型下，利用适当紧致正不变子水平集和惯性矩阵的有界、正定性质，可证明全局渐近位置镇定。

### 课堂反例：\(E\) 可以本身就是最大不变集

对 \(\dot x_1=-x_1,\ \dot x_2=0\)，选 \(V=(x_1^2+x_2^2)/2\)，有 \(\dot V=-x_1^2\)。此时

\[
E=\mathcal M=\{(0,x_2):x_2\in\mathbb R\}.
\]

轨迹趋近 \((0,x_2(0))\)，并非必然趋近原点。**这是用户在课堂思考题中正确给出的结论。**

## 4. 时变轨迹跟踪：两种动力学补偿

设 \(e=q-q_d(t)\)。仅使用 PD 与重力补偿，误差方程会出现与 \(\ddot q_d,\dot q_d\) 相关的强迫项，不能直接沿用固定目标的证明。

### A. 保留机械臂非线性结构的补偿

\[
\tau_A=M(q)\ddot q_d+C(q,\dot q)\dot q_d+g(q)-K_pe-K_d\dot e.
\]

理想闭环为

\[
M(q)\ddot e+C(q,\dot q)\dot e+K_d\dot e+K_pe=0.
\]

仍可得到 \(\dot V=-\dot e^TK_d\dot e\)。由于参考轨迹时变，误差系统一般是非自治系统；渐近收敛论证需补充有界性及 Barbalat 等条件，不能直接套用自治 LaSalle 定理。

### B. Computed Torque Control：精确反馈线性化

\[
\tau_B=M(q)(\ddot q_d-K_d\dot e-K_pe)+C(q,\dot q)\dot q+g(q).
\]

理想闭环为

\[
\boxed{\ddot e+K_d\dot e+K_pe=0.}
\]

控制器 B 依赖精确模型来实现线性化；但**不能仅由此断言它在所有情况下都比 A 不鲁棒**。

## 5. 模型失配与稳态误差

实际控制器采用估计模型 \(\hat M,\hat C,\hat g\)，定义

\[
v=\ddot q_d-K_d\dot e-K_pe,
\]

\[
\tau=\hat M(q)v+\hat C(q,\dot q)\dot q+\hat g(q).
\]

真实误差动力学为

\[
\ddot e+K_d\dot e+K_pe=M(q)^{-1}\big[\Delta M(q)v+\Delta C(q,\dot q)\dot q+\Delta g(q)\big],
\]

其中 \(\Delta M=\hat M-M\)、\(\Delta C=\hat C-C\)、\(\Delta g=\hat g-g\)。残差可能依赖状态，并非总是独立外部扰动。

### 本次实验：仅重力补偿偏差 10%

采用 \(\hat M=M\)、\(\hat C=C\)，\(\hat g(q)=1.1g(q)\)，固定目标 \(q_d=q^*\)。于是

\[
\ddot e+K_d\dot e+K_pe=0.1M(q)^{-1}g(q).
\]

在静态平衡 \(q_{\mathrm{ss}}=q^*+e_{\mathrm{ss}}\) 处：

\[
\boxed{M(q_{\mathrm{ss}})K_pe_{\mathrm{ss}}=0.1g(q_{\mathrm{ss}}).}
\]

注意：右端的重力误差**随实际位置变化**，不是常数。此式一般为关于 \(e_{\mathrm{ss}}\) 的隐式非线性方程。

## 6. NumPy 练习与实测结果

文件：`day17_gravity_mismatch_exercises.py`。

实现的 TODO：

1. `computed_torque_regulation`：计算 \(M(q)v+c(q,\dot q)+\text{gravity_scale}\cdot g(q)\)。
2. `compare_gravity_models`：分别仿真精确模型和 10% 重力失配，\(dt=0.002\,\mathrm{s}\)、总时长 \(8\,\mathrm{s}\)。
3. `equilibrium_residual`：计算 \(M(q_{\mathrm{ss}})K_p(q_{\mathrm{ss}}-q^*)-(\text{gravity_scale}-1)g(q_{\mathrm{ss}})\)。
4. `summarize_results`：统计位置误差、速度和稳态残差。

用户本地运行输出：

| 指标 | 数值 |
|---|---:|
| 精确模型最终位置误差 | \(2.697072\times10^{-12}\,\mathrm{rad}\) |
| 10% 重力失配最终位置误差 | \(0.04592891\,\mathrm{rad}\) |
| 失配情况下最终速度范数 | \(6.531959\times10^{-12}\) |
| 失配情况下静态平衡方程残差范数 | \(5.82384\times10^{-10}\) |

**4 组测试全部 PASS。** 数值结果说明：精确模型可几乎无误差地到达目标；重力失配产生非零稳态位置偏差，且仿真终态满足推导出的静态平衡方程。

代码检查备注：TODO 3 中计算了未使用的 `qd_ss` 和 `c`，可以删去；TODO 4 的 `gravity_scale` 参数可以传给 `equilibrium_residual`，而不固定写死 `1.1`，使函数对其他失配比例也通用。

## 7. 本日结论与下一步

- **机械臂动力学的结构性质**让 PD + 重力补偿可以通过能量函数分析。
- **LaSalle 原理的重点是寻找 \(E=\{\dot V=0\}\) 内的最大不变集**，不是假定 \(E\) 本身不变。
- **轨迹跟踪不同于位置镇定**；需要考虑参考轨迹的速度和加速度。
- **Computed Torque** 在模型精确时实现线性误差动力学，但模型失配会引入状态相关残差和稳态误差。
- **Day 18**：进一步比较 CTC 与动力学前馈控制，重点关注非零初始跟踪误差、模型失配及控制性能。
