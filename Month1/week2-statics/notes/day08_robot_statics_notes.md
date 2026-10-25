# Day 8 — Robot Statics（机器人静力学）

参考：*Modern Robotics*，Chapter 5.2。

## 1. Wrench：力与力矩

Wrench 将作用于刚体的力系表示为合力矩和合力（**力矩在前**）：

$$
F=\begin{bmatrix}m\\f\end{bmatrix}\in\mathbb R^6.
$$

其中，$m$ 是**关于所选坐标系原点**的合力矩，$f$ 是合力。若一个力 $f$ 作用于位置 $r$，则

$$
m=r\times f.
$$

多个力及独立力偶可以等效为

$$
f=\sum_i f_i,\qquad m_O=\sum_i r_i\times f_i+\sum_j m_{c,j}.
$$

更换力矩参考点会改变 $m$，但不会改变物理力系。

## 2. Twist、Wrench 与功率

采用角速度在前的 Twist 约定：

$$
V=\begin{bmatrix}\omega\\v\end{bmatrix}.
$$

在相互匹配的坐标表示下，瞬时功率为

$$
P=F^TV=m^T\omega+f^Tv.
$$

这里的 $v$ 是 Twist 的线性分量；对 Space Twist 而言，它不一定是末端原点的实际线速度。

## 3. Jacobian 转置与关节广义力

由速度运动学和功率关系：

$$
V=J(q)\dot q,\qquad P=F^TV=\tau_{\mathrm{ext}}^T\dot q.
$$

因此，外力产生的关节广义力为

$$
\boxed{\tau_{\mathrm{ext}}=J(q)^TF.}
$$

若忽略重力等其他载荷，静力平衡所需的驱动力矩为

$$
\tau_{\mathrm{motor}}=-\tau_{\mathrm{ext}}.
$$

若仅考虑作用在末端点的三维力 $f$，使用末端**位置 Jacobian**：

$$
\dot p=J_p(q)\dot q,\qquad \tau_{\mathrm{ext}}=J_p(q)^Tf.
$$

若改用 Space Jacobian，须将末端点的纯力表示为关于 Space 原点的完整 Wrench：

$$
F_s=\begin{bmatrix}p\times f\\f\end{bmatrix},\qquad J_s^TF_s=J_p^Tf.
$$

## 4. Wrench 的坐标变换与对偶性

设 $T_{sb}=(R,p)$ 表示 Body 坐标系相对于 Space 坐标系的位姿，则

$$
V_s=\operatorname{Ad}_T V_b,\qquad
\operatorname{Ad}_T=\begin{bmatrix}R&0\\[p]_\times R&R\end{bmatrix}.
$$

由于功率不依赖坐标系选择：

$$
F_s^TV_s=F_b^TV_b,
$$

得到 Wrench 的**对偶坐标变换**：

$$
\boxed{F_s=\operatorname{Ad}_T^{-T}F_b.}
$$

展开为

$$
f_s=Rf_b,\qquad m_s=Rm_b+p\times Rf_b.
$$

因此，纯平移不会改变合力的坐标分量，但通常会改变关于原点的力矩；若 $f=0$（纯力偶），平移不影响其力矩。

Space 和 Body Jacobian 满足

$$
J_s=\operatorname{Ad}_T J_b.
$$

故两种表示计算的关节广义力一致：

$$
\boxed{J_s^TF_s=J_b^TF_b.}
$$

## 5. NumPy 编程练习

实现了三个函数：

- `wrench_transform(T, F_b)`：通过 `np.linalg.solve(AdT.T, F_b)` 计算 $F_s$。
- `joint_torque(J, F)`：计算 $J^TF$。
- `verify_static_duality(T, J_b, F_b)`：验证两种坐标表示下的关节广义力一致。

测试结果：

- Wrench 平移变换：$F_b=(0,0,5,0,10,0)^T$，$p=(2,0,0)^T$，得到 $F_s=(0,0,25,0,10,0)^T$。
- 平面 2R 机械臂：$q=(0,\pi/2)^T$，末端受力 $f=(10,0)^T$，得到 $\tau_{\mathrm{ext}}=(-10,-10)^T$。
- 静力学对偶验证误差：$2.459\times10^{-15}$，属于浮点舍入误差。
- 纯力偶在不同平移下保持不变。

**Day 8 完成。**
