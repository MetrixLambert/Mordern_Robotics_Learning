# Day 9 — Jacobian 奇异性与可操作度

## 1. 速度可操作度椭球

任务速度与关节速度满足：

$$
\dot{x}=J(q)\dot{q},\qquad \|\dot{q}\|_2\leq 1.
$$

当 $J$ 行满秩时，末端可实现的速度集合为：

$$
\dot{x}^{T}(JJ^{T})^{-1}\dot{x}\leq 1.
$$

对 Jacobian 进行奇异值分解：

$$
J=U\Sigma V^{T}.
$$

速度椭球的主轴方向是 $U$ 的列向量，半轴长度是对应奇异值 $\sigma_i$。若某个奇异值为零，椭球退化，说明末端失去某一方向的瞬时运动能力。

## 2. 力可操作度椭球

由静力学关系：

$$
\tau=J^{T}f,\qquad \|\tau\|_2\leq 1,
$$

可得力的可行集合：

$$
f^{T}JJ^{T}f\leq 1.
$$

当 $J$ 为可逆方阵时，力椭球与速度椭球主轴方向相同，但半轴长度分别为：

$$
\text{速度：}\ \sigma_i,\qquad \text{力：}\ \frac{1}{\sigma_i}.
$$

因此，某方向运动能力越强，单位关节力矩约束下该方向的静态施力能力越弱。奇异时，力的可行集合可能沿某些方向无界；这不代表实际机构能承受无限载荷。

## 3. 可操作度与条件数

Yoshikawa 可操作度指标：

$$
w(q)=\sqrt{\det(JJ^{T})}=\prod_i\sigma_i.
$$

它与速度椭球的面积或体积成正比。对于方阵 Jacobian：

$$
w(q)=|\det J(q)|.
$$

Jacobian 条件数：

$$
\kappa(J)=\frac{\sigma_{\max}}{\sigma_{\min}}.
$$

- $w$：反映速度椭球的整体大小。
- $\kappa$：反映不同方向运动能力的均衡程度；$\kappa=1$ 表示各向同性。
- 两个指标的最优构型不一定相同。

## 4. 平面 2R 机械臂

设连杆长度 $L_1=2$、$L_2=1$，末端位置 Jacobian 为：

$$
J_p=\begin{bmatrix}
-2\sin q_1-\sin(q_1+q_2)&-\sin(q_1+q_2)\\
2\cos q_1+\cos(q_1+q_2)&\cos(q_1+q_2)
\end{bmatrix}.
$$

其行列式和可操作度：

$$
\det J_p=2\sin q_2,\qquad w(q)=2|\sin q_2|.
$$

因此，最大可操作度为 $2$，在以下构型取得：

$$
q_2=\frac{\pi}{2}+k\pi.
$$

它与 $q_1$ 无关。完全伸直（$q_2=0$）时，Jacobian 奇异，速度椭球退化为线段。

## 5. NumPy 练习结果

完成了三个函数：

- `ellipse_axes(J)`：利用 `np.linalg.svd()` 求速度椭球与力椭球的主轴、半轴。
- `yoshikawa_measure(J)`：计算奇异值乘积。
- `jacobian_condition_number(J)`：计算最大与最小奇异值之比，奇异时返回 `np.inf`。

**全部测试通过**，生成 `day09_ellipses.png` 和 `day09_metrics.png`。

在 $q_2\in[-\pi,\pi]$ 的网格搜索中：

$$
\max w=2,\qquad q_2\approx-1.57080;
$$

$$
\min\kappa\approx1.61807,\qquad q_2\approx-2.30383.
$$

**结论：** 最大化可操作度与最小化条件数不是同一个优化目标。

## 6. 延伸：冗余机械臂

对于冗余机械臂，可以利用 Jacobian 零空间，在保持末端瞬时任务速度不变的同时调整构型：

$$
\dot q=J^\dagger\dot x_d+(I-J^\dagger J)z.
$$

例如，选取 $z=k\nabla_q w$，可用于局部提升可操作度（需考虑关节约束和奇异性）。
