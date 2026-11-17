# ROBO_Learn — Month 2 · Day 15 总结

> 主题：Acrobot 欠驱动动力学、平衡点线性化、LQR 与部分反馈线性化（PFL）  
> 学习内容：Part 1–4（原计划 Day 15–18 的相关内容合并为 Day 15）  
> 主教材：Russ Tedrake, *Underactuated Robotics*, Chapter 3（Acrobot、Balancing、Partial Feedback Linearization）；LQR 参考 Chapter 8。  
> 后续：**Day 16 从 Energy Shaping / Swing-up 开始**。

## 1. 核心问题与结论

Acrobot 是固定肩关节、仅肘关节驱动的双连杆系统。它有两个广义自由度，但仅有一个输入。**欠驱动不等于不可控**：电机力矩可通过惯性耦合影响未驱动关节；直立平衡点附近的线性化系统可以完全可控。PFL 可以精确指定一个关节的加速度，但不能同时独立指定两个关节的加速度；被遗漏的内部动力学不能忽略。

## 2. Part 1：动力学与欠驱动约束

选取 $q=[q_1,q_2]^T$：$q_1$ 为第一连杆相对于竖直向下的绝对角，$q_2$ 为第二连杆相对于第一连杆的相对角。输入力矩为肘关节力矩 $\tau$。

$$M(q)\ddot q+C(q,\dot q)\dot q+G(q)=B\tau,\qquad B=\begin{bmatrix}0\\1\end{bmatrix}. $$

记 $h=C\dot q+G$，则

$$\begin{aligned}M_{11}\ddot q_1+M_{12}\ddot q_2+h_1&=0,\\M_{21}\ddot q_1+M_{22}\ddot q_2+h_2&=\tau.\end{aligned}$$

第一行是**未驱动关节的动力学约束**。虽然第一行没有 $\tau$，但 $\ddot q=M^{-1}(B\tau-h)$ 中第一关节的加速度仍可能随 $\tau$ 变化。对于二维正定对称 $M$，

$$\frac{\partial\ddot q_1}{\partial\tau}=-\frac{M_{12}}{\det M}. $$

本次示例参数：$m_1=m_2=1$，$l_1=1$，$l_{c1}=l_{c2}=0.5$，$I_1=I_2=1$，$g=9.81$。采用课堂所用的**近端关节转动惯量**约定，

$$M=\begin{bmatrix}I_1+I_2+m_2l_1^2+2m_2l_1l_{c2}\cos q_2&I_2+m_2l_1l_{c2}\cos q_2\\I_2+m_2l_1l_{c2}\cos q_2&I_2\end{bmatrix}. $$

$$C\dot q=\begin{bmatrix}-a\sin q_2(2\dot q_1\dot q_2+\dot q_2^2)\\a\sin q_2\dot q_1^2\end{bmatrix},\quad a=m_2l_1l_{c2}. $$

$$G=\begin{bmatrix}(m_1l_{c1}+m_2l_1)g\sin q_1+m_2gl_{c2}\sin(q_1+q_2)\\m_2gl_{c2}\sin(q_1+q_2)\end{bmatrix}. $$

**课堂思考题**：希望 $\ddot q_1=0$，若 $M_{12}\neq0$，则 $\ddot q_2=-h_1/M_{12}$，并且 $\tau=h_2-(M_{22}/M_{12})h_1$。若 $M_{12}=0$，仅在 $h_1=0$ 时该零加速度条件可自然成立；一般不能通过 $\tau$ 强制实现。若 $M_{12}$ 接近零，控制力矩可能很大。本次参数下 $M_{12}=1+0.5\cos q_2\ge0.5$，不存在该奇异性。

## 3. Part 2：平衡点、线性化与可控性

定义 $x=[q_1,q_2,\dot q_1,\dot q_2]^T$。零输入静态平衡点包括悬垂 $q^*=(0,0)$ 和直立 $q^*=(\pi,0)$，均有 $\dot q^*=0,\tau^*=0$。一般静态平衡条件为 $G(q^*)=B\tau^*$，所以非零平衡力矩时也可能有其他构型。

在零速度、零输入平衡点，令 $K_g=\partial G/\partial q|_{q^*}$，有

$$\delta\dot x=A\delta x+B_{\rm lin}\delta\tau,\quad A=\begin{bmatrix}0&I\\-M_*^{-1}K_g&0\end{bmatrix},\quad B_{\rm lin}=\begin{bmatrix}0\\M_*^{-1}B\end{bmatrix}. $$

对本次参数的直立平衡点：

$$M_*=\begin{bmatrix}4&1.5\\1.5&1\end{bmatrix},\qquad K_g=\begin{bmatrix}-19.62&-4.905\\-4.905&-4.905\end{bmatrix}. $$

$$A\approx\begin{bmatrix}0&0&1&0\\0&0&0&1\\7.0071&-1.4014&0&0\\-5.6057&7.0071&0&0\end{bmatrix},\quad B_{\rm lin}\approx\begin{bmatrix}0\\0\\-0.8571\\2.2857\end{bmatrix}. $$

开环直立平衡点不稳定，但 $\operatorname{rank}[B_{\rm lin},AB_{\rm lin},A^2B_{\rm lin},A^3B_{\rm lin}]=4$。**一个输入可以控制四维线性状态，不等于可以瞬时任意指定两个加速度。**

## 4. Part 3：LQR 局部平衡控制

代价函数

$$J=\int_0^\infty(\delta x^TQ\delta x+\delta\tau^TR\delta\tau)\,dt,$$

Riccati 方程

$$A^TP+PA-PB_{\rm lin}R^{-1}B_{\rm lin}^TP+Q=0,$$

反馈律

$$\tau=-K(x-x^*),\qquad K=R^{-1}B_{\rm lin}^TP.$$

实验使用 $Q=\operatorname{diag}(100,100,10,10)$、$R=1$。在非线性模型中比较 $2^\circ,5^\circ,10^\circ$ 的初始偏差，并测试输入饱和 $|\tau|\le\tau_{\max}$。要点：线性化闭环 Hurwitz 可推出**非线性系统局部**渐近稳定；不能保证从自然悬垂位置直接起摆，也不能忽略力矩限制。$V=\delta x^TP\delta x$ 是局部 Lyapunov 分析的自然候选。

## 5. Part 4：部分反馈线性化 PFL

完全驱动机械系统可以通过 $\tau=M(q)v+h(q,\dot q)$ 得到 $\ddot q=v$（Computed Torque / Feedback Linearization）。但 Acrobot 只有一个电机，无法独立指定两个关节的加速度。

### 5.1 Collocated PFL：指定驱动关节加速度

令 $\ddot q_2=v$，由第一行消去 $\ddot q_1$，得到

$$\boxed{\tau=\left(M_{22}-M_{21}M_{11}^{-1}M_{12}\right)v+h_2-M_{21}M_{11}^{-1}h_1}. $$

其中 $M_{22}-M_{21}M_{11}^{-1}M_{12}$ 为 Schur 补；若 $M$ 正定，则其严格为正。此输入变换在有效构型空间中没有 $M_{12}=0$ 的奇异问题。

### 5.2 Non-collocated PFL：指定未驱动关节加速度

令 $\ddot q_1=v$，由第一行消去 $\ddot q_2$，得到

$$\boxed{\tau=\left(M_{21}-\frac{M_{22}M_{11}}{M_{12}}\right)v+h_2-\frac{M_{22}}{M_{12}}h_1},\qquad M_{12}\neq0. $$

当 $v=0$ 时，退化为 Part 1 的思考题。多自由度推广中对应强惯性耦合/解耦矩阵的秩条件。

### 5.3 内部动力学与零动态

令输出 $y=q_1-q_{1d}$，Non-collocated PFL 可令 $\ddot y=v$。即使选择 $v=-k_py-k_d\dot y$ 使输出渐近收敛，剩余的 $q_2$ 动力学仍可能不稳定。**输出线性化不等于全状态稳定。**当 $y=\dot y=0$ 时剩余的内部动力学称为零动态，与非线性系统的 minimum-phase 问题相关。

## 6. 已完成实验与待验证项目

- 实现 `dynamics_terms(q,dq)`、`acceleration(q,dq,tau)`；验证第一行欠驱动动力学残差约为零。
- 在直立平衡点采用中心差分计算 $(A,B_{\rm lin})$；检查特征值和可控性秩。
- 用 `scipy.linalg.solve_continuous_are` 计算 LQR；用 `solve_ivp` 对原始非线性模型仿真；比较初始扰动、$R$ 权重及力矩饱和。
- 分别实现两种 PFL，检查 `ddq[1] ≈ v`（collocated）和 `ddq[0] ≈ v`（non-collocated）。

> 以上是课堂提供的实验方案与理论预期；不把尚未由学习者运行的实验记为已验证成功。

## 7. Day 15 五条核心结论

1. 欠驱动意味着输入维数不足以瞬时独立指定所有广义加速度，并不意味着线性化系统不可控。
2. 惯性耦合 $M_{12}$ 使肘关节力矩能够影响未驱动肩关节的加速度。
3. LQR 能稳定直立平衡点的局部邻域，无法自动解决大范围起摆。
4. PFL 精确线性化选定的输出，而不是整个欠驱动系统；collocated 与 non-collocated 的可逆性条件不同。
5. 必须区分输出误差动力学与内部/零动态的稳定性。

## 8. 下一步：Day 16

从 *Underactuated Robotics* §3.6（Swing-up Control）开始，重点研究能量整形。教材 §3.6.3 给出 Acrobot 的方案：先用 collocated PFL 令 $\ddot q_2=v$，再利用总能量误差构造 $v=-k_1q_2-k_2\dot q_2+k_3\dot q_1(E-E_d)$。教材明确指出该控制器的额外 PD 项妨碍了对目标同宿轨道渐近收敛的直接证明；实际效果需要仿真检查。随后研究 Swing-up 与 LQR 的切换。
