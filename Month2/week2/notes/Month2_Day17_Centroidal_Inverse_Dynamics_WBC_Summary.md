# ROBO_Learn — Month 2 · Day 17 总结

**主题**：高自由度腿足机器人的质心动力学、逆动力学 QP 与 Whole-Body Control（WBC）  
**课程调整**：将原 Day 18 的逆动力学 QP 内容并入 Day 17 Part 2；Day 17 分为三部分。  
**教材对应**：Russ Tedrake, *Underactuated Robotics*, Chapter 5（Highly-articulated Legged Robots）；QP/WBC 具体数学形式结合此前课程和工程实践补充。教材 §5.4 WBC 在所附版本中尚未展开。

## Part 1：从浮动基座动力学到质心动力学

### 完整动力学

在局部广义坐标（或采用广义速度 `v` 的适当表示）下：

\[
M(q)\dot v+h(q,v)=S^\top\tau+J_c(q)^\top\lambda.
\]

- `S` 为驱动选择矩阵；浮动基座无直接驱动。
- `J_c` 是接触 Jacobian；`λ` 为地面施加于机器人的接触力或 wrench。
- 固定无滑动接触模式下，约束为 `J_c \dot v+\dot J_c v=0`（需在一致的速度坐标中解释）。

### 质心线动量与角动量

\[
m\ddot p_c=mg_{\rm vec}+\sum_i f_i,\qquad
\dot L_c=\sum_i(p_i-p_c)\times f_i+\sum_i\mu_i.
\]

其中 `μ_i` 为接触力偶；点接触模型中通常取零。质心动量矩阵 `A_G(q)` 满足

\[
\begin{bmatrix}l_c\\L_c\end{bmatrix}=A_G(q)v.
\]

**关键认识**：质心动力学描述整体运动，但不自动保证接触力、驱动力矩及内部关节运动的可实现性。

### 课堂思考

`\ddot p_c=0`、`\dot L_c=0` 不足以保证机器人静止：接触力可能违反摩擦/单边接触/力矩约束；内部关节也可能移动。零动量、固定接触的瞬时内部速度满足

\[
v\in\ker A_G(q)\cap\ker J_c(q).
\]

这个空间可能非平凡，也可能只有零向量。

## Part 2：逆动力学 QP

### 决策变量

\[
z=\begin{bmatrix}\dot v^\top&\tau^\top&\lambda^\top\end{bmatrix}^\top.
\]

### 硬约束

1. **完整动力学**：`[M, -S^T, -J_c^T] z = -h`。
2. **固定接触加速度**：`J_c \dot v+\dot J_c v=0`。
3. **单边接触**：`f_{iz}\ge0`（地面法向取正）。
4. **摩擦约束**：`\sqrt{f_{ix}^2+f_{iy}^2}\le\mu f_{iz}`；标准线性 QP 可使用保守的多面体内近似。
5. **电机力矩界**：`\tau_{\min}\le\tau\le\tau_{\max}`。

### 软目标

质心加速度由运动学给出：

\[
\ddot p_c=J_{\rm CoM}\dot v+\dot J_{\rm CoM}v.
\]

典型优化：

\[
\min_z\;\tfrac12\|J_{\rm CoM}\dot v+\dot J_{\rm CoM}v-\ddot p_c^d\|_{W_c}^2
+\tfrac{\epsilon_\tau}{2}\|\tau\|^2
+\tfrac{\epsilon_\lambda}{2}\|\lambda\|^2
\]

并满足上述硬约束。整理为标准形式：

\[
\min_z\tfrac12 z^\top H z+g^\top z,\quad A_{eq}z=b_{eq},\quad A_{ineq}z\le b_{ineq}.
\]

**软目标不可实现 ≠ QP 不可行**：只有硬约束集合为空才必然不可行。

### 课堂数值例子

质量 `m=20 kg`，摩擦系数 `μ=0.2`，保持竖直加速度为零：

\[
|\ddot x_c|\le\mu g=1.962\;\mathrm{m/s^2}.
\]

期望 `3 m/s²` 超出摩擦极限。若仅优化水平加速度跟踪误差，且其他约束不更严格，最优值为 `1.962 m/s²`；加入其他约束/正则化后可能更小。课堂选择题选 **B**，正确。

## Part 3：Whole-Body Control

### 统一的任务空间加速度

对任务变量 `y_i=f_i(q)`，在合适坐标表示下：

\[
\ddot y_i=J_i\dot v+\dot J_i v.
\]

期望任务加速度可采用前馈 + PD：

\[
\ddot y_i^{cmd}=\ddot y_i^d+K_{d,i}(\dot y_i^d-\dot y_i)+K_{p,i}(y_i^d-y_i).
\]

加权 WBC 的典型目标：

\[
\min_z\sum_i\tfrac12w_i\|J_i\dot v+\dot J_i v-\ddot y_i^{cmd}\|^2+\tfrac12\|z\|_R^2.
\]

物理约束保持为硬约束；质心、躯干、摆动足、关节姿态等通常为任务目标。

### 加权 QP 与层级 QP

两个冲突任务 `a=2`、`a=-1`，加权优化为

\[
\min_a\tfrac12w_A(a-2)^2+\tfrac12w_B(a+1)^2
\Rightarrow a^*=\frac{2w_A-w_B}{w_A+w_B}.
\]

- **加权 QP**：有限权重代表折中；高权重不等于绝对优先。
- **层级 QP**：低优先级任务不能恶化高优先级任务的最优结果；若存在零空间，低优先级任务可在其中实现。
- **接触约束**：在假设固定、无滑动接触的模式下应作为硬约束，而非仅用低误差的软目标代替。

## 与后续课程的联系

- **逆动力学 QP/WBC**：在当前状态、当前接触模式下优化瞬时加速度、力矩和接触力。
- **轨迹优化（Day 18）**：同时优化未来多个时刻的状态和控制，使其满足跨时间的动力学约束。
- **MPC（后续）**：反复求解有限时域轨迹优化，只执行首段控制并滚动更新。

## Day 17 验收

- 能区分质心平衡条件与全身静止/接触可行性。
- 能构造包含动力学、接触、摩擦和力矩限制的逆动力学 QP。
- 能区分硬约束、软任务、加权优先级和严格层级优先级。
- 用户已学习过 WBC 数值实现，因此本日不重复实验。

**结论**：质心动力学负责描述整体目标，完整动力学与接触约束决定可行性，WBC 在这些约束下协调多个任务；下一步从瞬时优化过渡到有限时域轨迹优化。
