# ROBO_Learn · Day 14 学习总结

## Euler–Lagrange 方程、关节空间动力学与重力建模

> **学习主线**：从能量描述推导关节空间动力学；理解 \(M(q)\)、\(C(q,\dot q)\)、\(g(q)\) 的来源；通过 RNEA 交叉验证解析结果。  
> **参考**：Lynch & Park, *Modern Robotics*（机器人动力学章节）；配套练习 `day14_lagrange_exercises.py`。  
> **完成情况**：TODO 1–5 均已完成；四组主要测试各 100 次随机采样全部通过；额外的正/逆动力学一致性测试也已通过。

---

## 1. 广义坐标、广义力与虚功

对 \(n\) 自由度机械臂，选取关节位置作为广义坐标：

$$
q=(q_1,\ldots,q_n)^T,\qquad \dot q,\ddot q\in\mathbb R^n.
$$




广义力 \(Q\) 由虚功定义：

$$
\delta W=Q^T\delta q.
$$




对于关节驱动力矩，\(Q=\tau\)。如果末端受到空间 Wrench \(F\)，末端虚位移（瞬时 Twist 的对应微小运动）由雅可比映射，故

$$
\delta W=F^T J(q)\delta q
 =(J(q)^TF)^T\delta q,\qquad
\boxed{\tau_{\mathrm{ext}}=J(q)^TF}.
$$

这与前几天学过的 Jacobian 转置、力与运动的对偶关系相衔接。需要保持 \(J,F\) 的参考坐标系和排列顺序一致。

## 2. Lagrangian 与 Euler–Lagrange 方程

定义动能 \(K\)、势能 \(U\)，Lagrangian 为

$$
\mathcal L(q,\dot q)=K(q,\dot q)-U(q).
$$




Euler–Lagrange 方程：

$$
\boxed{
\frac{d}{dt}\frac{\partial\mathcal L}{\partial\dot q_i}
-\frac{\partial\mathcal L}{\partial q_i}=Q_i,\quad i=1,\ldots,n.
}
$$




机器人通常可将动能写为

$$
\boxed{K(q,\dot q)=\frac12\dot q^TM(q)\dot q},
$$




其中惯性矩阵 \(M(q)\) 对称且正定（独立关节坐标、非退化物理模型）。

若 \(U=U(q)\)，代入 Euler–Lagrange 方程可得到

$$
\boxed{M(q)\ddot q+c(q,\dot q)+g(q)=\tau},
$$

其中 \(c\) 汇集速度二次项，\(g(q)=\nabla_q U(q)\)。也常写成

$$
\boxed{M(q)\ddot q+C(q,\dot q)\dot q+g(q)=\tau},
\quad c(q,\dot q)=C(q,\dot q)\dot q.
$$




**提醒**：\(c\) 是向量；\(C\) 是矩阵。仅要求 \(C\dot q=c\) 时，\(C\) 的表示一般不唯一；Christoffel 构造给出一种具有重要结构性质的标准选取。

### 2.1 一自由度热身

若

\[
K=\tfrac12(2+\cos q)\dot q^2,\qquad U=3(1-\cos q),
\]

则

\[
\frac{\partial\mathcal L}{\partial\dot q}=(2+\cos q)\dot q,
\]

\[
\frac{d}{dt}\frac{\partial\mathcal L}{\partial\dot q}
=(2+\cos q)\ddot q-\sin q\,\dot q^2,
\]

\[
\frac{\partial\mathcal L}{\partial q}
=-\tfrac12\sin q\,\dot q^2-3\sin q.
\]

因此

\[
\boxed{(2+\cos q)\ddot q-\tfrac12\sin q\,\dot q^2+3\sin q=\tau}.
\]

易错点：对 \((2+\cos q)\dot q\) 求时间导数时，必须记得 \(q=q(t)\)，应用乘积法则和链式法则。

---

## 3. 二连杆平面机械臂模型

**约定**：机械臂在 \(XY\) 竖直平面运动；\(X\) 水平向右，\(Y\) 竖直向上；关节旋转轴为 \(+Z\)；\(q_1\) 相对水平正方向，\(q_2\) 为第二连杆相对第一连杆的角度。

| 参数          | 数值            | 含义               |
| ----------- | -------------:| ---------------- |
| \(L_1,L_2\) | \(2,1\)       | 连杆长度             |
| \(m_1,m_2\) | \(2,1\)       | 连杆质量             |
| \(r_1,r_2\) | \(1,0.5\)     | 关节到对应质心距离        |
| \(I_1,I_2\) | \(0.25,0.10\) | 绕质心 \(Z\) 轴的转动惯量 |
| \(g_0\)     | \(9.81\)      | 重力加速度大小          |

为避免与重力向量 \(g(q)\) 混淆，下文用 \(g_0\) 表示重力加速度的标量大小。

两连杆质心位置：

\[
p_{C_1}=\begin{bmatrix}r_1\cos q_1\\r_1\sin q_1\end{bmatrix},
\]

\[
p_{C_2}=
\begin{bmatrix}
L_1\cos q_1+r_2\cos(q_1+q_2)\\
L_1\sin q_1+r_2\sin(q_1+q_2)
\end{bmatrix}.
\]

角速度：

\[
\omega_1=\dot q_1,\qquad \omega_2=\dot q_1+\dot q_2.
\]

总动能由两连杆的平动和转动动能组成：

\[
K=\tfrac12m_1\|\dot p_{C_1}\|^2
+\tfrac12I_1\dot q_1^2
+\tfrac12m_2\|\dot p_{C_2}\|^2
+\tfrac12I_2(\dot q_1+\dot q_2)^2.
\]

整理成 \(K=\frac12\dot q^TM(q)\dot q\) 得到：

\[
\begin{aligned}
M_{11}&=I_1+I_2+m_1r_1^2
+m_2(L_1^2+r_2^2+2L_1r_2\cos q_2),\\
M_{12}=M_{21}&=I_2+m_2(r_2^2+L_1r_2\cos q_2),\\
M_{22}&=I_2+m_2r_2^2.
\end{aligned}
\]

代入数值：

\[
\boxed{
M(q)=\begin{bmatrix}
6.60+2\cos q_2&0.35+\cos q_2\\
0.35+\cos q_2&0.35
\end{bmatrix}.
}
\]

**理解**：惯性矩阵的非对角元素反映关节之间的惯性耦合；即使 \(\ddot q_2=0\)，\(\ddot q_1\) 也可能需要在第二关节施加力矩。

---

## 4. Christoffel 符号与速度耦合项

从 \(K=\frac12\sum_{jk}M_{jk}(q)\dot q_j\dot q_k\) 出发，Euler–Lagrange 方程中的速度二次项可写为

\[
c_i(q,\dot q)=\sum_{j=1}^n\sum_{k=1}^n
\Gamma_{ijk}(q)\dot q_j\dot q_k,
\]

其中**第一类 Christoffel 符号**定义为

\[
\boxed{
\Gamma_{ijk}=
\frac12\left(
\frac{\partial M_{ij}}{\partial q_k}
+\frac{\partial M_{ik}}{\partial q_j}
-\frac{\partial M_{jk}}{\partial q_i}
\right).
}
\]

满足 \(\Gamma_{ijk}=\Gamma_{ikj}\)。定义标准速度耦合矩阵

\[
\boxed{C_{ij}(q,\dot q)=\sum_k\Gamma_{ijk}(q)\dot q_k},
\]

于是 \(c=C\dot q\)。

对当前二连杆，令

\[
h(q_2)=m_2L_1r_2\sin q_2=\sin q_2.
\]

关键系数包括

\[
\Gamma_{112}=\Gamma_{121}=-h,\quad
\Gamma_{122}=-h,\quad
\Gamma_{211}=h,
\]

其余与 \(c\) 有关的系数为零。由此

\[
\boxed{
C(q,\dot q)=
\begin{bmatrix}
-h\dot q_2&-h(\dot q_1+\dot q_2)\\
h\dot q_1&0
\end{bmatrix}
}
\]

以及

\[
\boxed{
c(q,\dot q)=h
\begin{bmatrix}
-2\dot q_1\dot q_2-\dot q_2^2\\
\dot q_1^2
\end{bmatrix}.
}
\]

这里的速度项并非单纯的“阻尼”：它来自机械臂几何结构随关节位置变化而产生的惯性耦合。

---

## 5. 结构性质：\(\dot M-2C\) 的斜对称性

由于 \(M\) 只依赖 \(q_2\)，

\[
\dot M(q)=
\begin{bmatrix}
-2h\dot q_2&-h\dot q_2\\
-h\dot q_2&0
\end{bmatrix}.
\]

因此

\[
\boxed{
\dot M-2C=
\begin{bmatrix}
0&h(2\dot q_1+\dot q_2)\\
-h(2\dot q_1+\dot q_2)&0
\end{bmatrix}.
}
\]

显然

\[
\boxed{(\dot M-2C)^T=-(\dot M-2C)}.
\]

这一性质依赖于前述 Christoffel 方式构造的 \(C\)；不能对任意满足 \(C\dot q=c\) 的矩阵表示直接作同样断言。

对任意斜对称矩阵 \(S\)，有 \(x^TSx=0\)，因此

\[
\boxed{\dot q^T(\dot M-2C)\dot q=0}.
\]

### 5.1 能量平衡与无源性

定义机械能

\[
E(q,\dot q)=\frac12\dot q^TM(q)\dot q+U(q).
\]

求导并代入 \(M\ddot q+C\dot q+g=\tau\)，且 \(g=\nabla U\)：

\[
\begin{aligned}
\dot E
&=\dot q^TM\ddot q+\tfrac12\dot q^T\dot M\dot q+\dot q^Tg\\
&=\dot q^T\tau+\tfrac12\dot q^T(\dot M-2C)\dot q\\
&=\boxed{\dot q^T\tau}.
\end{aligned}
\]

在无摩擦、无其他外部非保守力的理想模型中：

- \(\tau=0\)：机械能守恒；
- \(\dot E=\tau^T\dot q\)：输入力矩 \(\tau\) 与输出关节速度 \(\dot q\) 构成机械功率端口；
- 这也是机器人控制中能量法、无源性分析和 Lyapunov 方法的重要起点。

与 Day 13 联系：无重力、零输入仿真中动能应守恒；数值积分导致的微小漂移不代表连续时间动力学违反守恒。

---

## 6. 势能与重力项

在竖直平面、\(Y\) 向上时：

\[
U(q)=m_1g_0r_1\sin q_1+
m_2g_0\left[L_1\sin q_1+r_2\sin(q_1+q_2)\right].
\]

由 \(g(q)=\nabla_q U(q)\)：

\[
\boxed{
g(q)=
\begin{bmatrix}
(m_1r_1+m_2L_1)g_0\cos q_1
+m_2g_0r_2\cos(q_1+q_2)\\
m_2g_0r_2\cos(q_1+q_2)
\end{bmatrix}.
}
\]

代入参数：

\[
\boxed{
g(q)=g_0\begin{bmatrix}
4\cos q_1+0.5\cos(q_1+q_2)\\
0.5\cos(q_1+q_2)
\end{bmatrix}.
}
\]

物理理解：第一关节需要补偿两根连杆的重力力矩；第二关节只需要补偿第二连杆的重力力矩。这里 \(g(q)\) 是**平衡重力所需的广义驱动力矩项**，而非重力本身施加的广义力；重力施加的广义力是 \(-g(q)\)。

### 完整二连杆动力学

\[
\boxed{
\underbrace{\begin{bmatrix}
6.60+2\cos q_2&0.35+\cos q_2\\
0.35+\cos q_2&0.35
\end{bmatrix}}_{M(q)}\ddot q
+
\underbrace{\sin q_2
\begin{bmatrix}
-2\dot q_1\dot q_2-\dot q_2^2\\
\dot q_1^2
\end{bmatrix}}_{c(q,\dot q)}
+
\underbrace{g_0
\begin{bmatrix}
4\cos q_1+0.5\cos(q_1+q_2)\\
0.5\cos(q_1+q_2)
\end{bmatrix}}_{g(q)}
=\tau.
}
\]

正动力学：

\[
\boxed{\ddot q=M(q)^{-1}\bigl[\tau-C(q,\dot q)\dot q-g(q)\bigr]}.
\]

数值实现时优先使用 `np.linalg.solve(M, rhs)`，而不是显式计算逆矩阵。

---

## 7. 为什么 RNEA 加重力只需修改前向递推？

### 7.1 物理出发点

对质点：

\[
ma=f_{\mathrm{act}}+m g_{\mathrm{vec}}
\quad\Longrightarrow\quad
f_{\mathrm{act}}=m(a-g_{\mathrm{vec}}).
\]

因此计算驱动力时，可以用“真实加速度减去重力加速度”的**等效加速度**。

对固定基座机械臂，真实基座加速度为零。若世界坐标系中

\[
g_{\mathrm{vec}}=(0,-g_0,0)^T,
\]

在 RNEA 前向递推中使用算法初值

\[
\boxed{
\dot V_0=
\begin{bmatrix}0_{3\times1}\\-g_{\mathrm{vec}}\end{bmatrix}
=
(0,0,0,0,g_0,0)^T.
}
\]

注意：这不是说物理基座真的在加速，而是引入重力的计算技巧。

### 7.2 前向、反向递推的分工

前向递推：

\[
V_i=\operatorname{Ad}_{T_{i-1,i}^{-1}}V_{i-1}+A_i\dot q_i,
\]

\[
\dot V_i=
\operatorname{Ad}_{T_{i-1,i}^{-1}}\dot V_{i-1}
+A_i\ddot q_i+\operatorname{ad}_{V_i}(A_i\dot q_i).
\]

改变 \(\dot V_0\) 后，等效重力加速度通过坐标变换传递到每根连杆的 \(\dot V_i\)。

反向递推仍为

\[
F_i=G_i\dot V_i-\operatorname{ad}_{V_i}^TG_iV_i
+\operatorname{Ad}_{T_{i,i+1}^{-1}}^TF_{i+1},
\qquad \tau_i=A_i^TF_i,
\]

（此处假设无其他显式外部 Wrench）。由于 \(\dot V_i\) 已包含等效重力效应，\(G_i\dot V_i\) 便已包含相应的重力补偿 Wrench，所以**无需再次修改反向递推**。

另一种等价实现：保持真实基座加速度 \(\dot V_0=0\)，在反向递推中对每根连杆显式减去重力 Wrench \(F_{g,i}\)：

\[
F_i=G_i\dot V_i-\operatorname{ad}_{V_i}^TG_iV_i
-F_{g,i}+\operatorname{Ad}_{T_{i,i+1}^{-1}}^TF_{i+1}.
\]

两种方法在固定基座、均匀重力场且坐标系处理一致时等价。**不要同时用等效基座加速度和显式重力 Wrench 重复计算重力。**

---

## 8. NumPy 实践与测试结果

配套文件：`day14_lagrange_exercises.py`。

实现内容：

1. `lagrange_mass_matrix(robot, q)`：解析惯性矩阵；
2. `christoffel_C(robot, q, qd)`：速度耦合矩阵；
3. `potential_energy(robot, q)`：重力势能；
4. `gravity_vector(robot, q)`：势能梯度；
5. `lagrange_inverse_dynamics(robot, q, qd, qdd)`：逆动力学；
6. `lagrange_forward_dynamics(robot, q, qd, tau)`：正动力学。

已完成的验证：

| 验证            | 方法                                                  | 结果      |
| ------------- | --------------------------------------------------- | ------- |
| \(M(q)\)      | 与无重力 RNEA 的单位加速度列比较                                 | 100 次通过 |
| \(C\dot q\)   | 与无重力 RNEA 的零加速度输出比较                                 | 100 次通过 |
| \(\dot M-2C\) | 有限差分计算 \(\dot M\)，验证斜对称                             | 100 次通过 |
| \(g(q)\)      | 有限差分势能梯度 + 有重力 RNEA                                 | 100 次通过 |
| 完整逆动力学        | Lagrange 与有重力 RNEA 比较                               | 100 次通过 |
| 正/逆动力学        | 给定目标 \(\ddot q\)，逆动力学求 \(\tau\)，再正动力学还原 \(\ddot q\) | 额外测试通过  |

**实现上的小提醒**：如果函数提供了 `gravity` 参数，内部调用 `gravity_vector(robot, q, gravity=gravity)`，而不是省略参数；否则默认 \(9.81\) 的测试虽然能通过，但调用者传入其他重力加速度时结果会不一致。正动力学和逆动力学两个函数都应注意这一点。

## 9. 与前后学习内容的衔接

- **Day 11–12**：从刚体空间惯性与 Newton–Euler 方程出发，建立 RNEA。
- **Day 13**：由 RNEA 数值提取 \(M\)、速度偏置项，并进行正动力学与数值积分。
- **Day 14**：从能量与 Euler–Lagrange 方程独立推导 \(M,C,g\)，与 RNEA 相互验证；理解斜对称性、能量平衡及重力的两种引入方法。
- **后续控制学习**：重力补偿、计算力矩控制、PD+重力补偿、能量/Lyapunov 分析，都建立在这里的动力学结构之上。

### 自查问题

1. 为什么 \(M(q)\) 必须对称？为什么通常正定？
2. 为什么 \(C(q,\dot q)\) 本身不是唯一的，但 \(C\dot q\) 可以唯一确定？
3. 如何由 Christoffel 符号推导 \(\dot M-2C\) 的斜对称性？
4. 为什么 \(\tau=0\) 时理想机械臂的总机械能守恒，但动能不一定守恒？
5. 为什么 RNEA 可以通过设定 \(\dot V_0=(0,-g_{\mathrm{vec}})^T\) 加入重力，而反向递推不变？
