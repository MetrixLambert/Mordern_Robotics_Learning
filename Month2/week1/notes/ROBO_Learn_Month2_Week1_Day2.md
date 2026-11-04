# ROBO_Learn 学习笔记｜Month 2 · Week 1 · Day 2

**主题：Contact Constraints & Constrained Dynamics（接触约束与约束动力学）**  
**学习状态：理论学习及问答验收完成；MuJoCo 接触实验尚未进行（计划 Day 3）。**

## 1. 从浮动基座到接触约束

无接触浮动基座机器人动力学为

$$
M(q)\dot v+h(q,v)=S^T\tau.
$$

当足端与地面保持**固定、不滑动且不脱离**的接触时，足端世界坐标位置满足完整约束：

$$
p_c(q)=p_c^*.
$$

在局部坐标且 $v=\dot q$ 时，定义接触 Jacobian：

$$
J_c(q)=\frac{\partial p_c(q)}{\partial q},\qquad J_c(q)\dot q=0.
$$

对于一般浮动基座速度 $v$（可能使用基座 twist 而非姿态参数导数），若 $\dot q=N(q)v$，则相应速度 Jacobian 为 $J_c^{(v)}=(\partial p_c/\partial q)N(q)$，速度约束写作 $J_c^{(v)}v=0$。

**物理含义：**接触点固定不代表整个机器人静止；广义速度只能位于 $\ker J_c$ 中。若 $q\in\mathbb R^n$ 且 $J_c$ 的秩为 $r$，瞬时允许运动空间的维数为 $n-r$。

## 2. 平面单连杆接触例子

设

$$
q=[x_b,z_b,\theta_b,q_1,q_2]^T,\qquad
\alpha=\theta_b+q_1,
$$

$$
p_c(q)=\begin{bmatrix}x_b+l\cos\alpha\\z_b+l\sin\alpha\end{bmatrix}.
$$

注意：$\alpha$ 是连杆相对于**世界水平轴**的角度，因此 $\alpha=0$ 时连杆水平向右。

接触 Jacobian：

$$
J_c=\begin{bmatrix}
1&0&-l\sin\alpha&-l\sin\alpha&0\\
0&1&l\cos\alpha&l\cos\alpha&0
\end{bmatrix}.
$$

在 $l=1,\alpha=0,\dot\theta_b=0,\dot q_1=1,\dot q_2=0$ 时：

$$
J_c=\begin{bmatrix}1&0&0&0&0\\0&1&1&1&0\end{bmatrix},
\quad \dot x_b=0,\quad\dot z_b=-1.
$$

因为 $\operatorname{rank}J_c=2$，五维广义速度中有 **3 维**满足接触约束。

## 3. 加速度层面的接触约束

对速度约束求导：

$$
\boxed{J_c(q)\ddot q+\dot J_c(q,\dot q)\dot q=0.}
$$

这里 $\dot J_c\dot q$ 不可随意忽略；它描述随构型变化产生的速度相关加速度。

对于上面的连杆：

$$
\ddot p_c=\begin{bmatrix}
\ddot x_b-l\cos\alpha\,\dot\alpha^2-l\sin\alpha\,\ddot\alpha\\
\ddot z_b-l\sin\alpha\,\dot\alpha^2+l\cos\alpha\,\ddot\alpha
\end{bmatrix}.
$$

在 $l=1,\alpha=0,\dot\alpha=1,\ddot\alpha=0$ 时，要求 $\ddot p_c=0$，得

$$
\boxed{\ddot x_b=1,\qquad\ddot z_b=0.}
$$

若错误地忽略 $\dot J_c\dot q$，会求得 $\ddot x_b=\ddot z_b=0$，此时真实足端加速度却是 $[-1,0]^T$，违反固定接触条件。

## 4. 接触力为何以 $J_c^T\lambda$ 进入动力学

由虚功关系：

$$
\delta p_c=J_c\delta q,\qquad
\delta W=\lambda^T\delta p_c
=(J_c^T\lambda)^T\delta q.
$$

因此接触广义力为

$$
\boxed{Q_c=J_c^T\lambda.}
$$

接触动力学：

$$
\boxed{M(q)\ddot q+h(q,\dot q)=S^T\tau+J_c(q)^T\lambda.}
$$

这里采用 $v=\dot q$ 的局部表示；若采用一般广义速度，则对应 $M\dot v+h=S^T\tau+J_c^{(v)T}\lambda$，加速度约束也需对 $J_c^{(v)}v=0$ 求导。

$\tau$ 是直接指定的关节驱动力矩；$\lambda$ 是为满足接触模式而出现的约束反作用力，不能任意独立指定。

**简单例子：**质量 $m=2\,\mathrm{kg}$，重力加速度 $g=10\,\mathrm{m/s^2}$，竖直向上接触力 $\lambda_z=10\,\mathrm N$，则

$$
\ddot z=\frac{\lambda_z-mg}{m}=-5\,\mathrm{m/s^2}.
$$

若接触点固定，必须有 $\ddot z=0$，因此 $\lambda_z=20\,\mathrm N$。

## 5. 增广约束动力学与接触力求解

将动力学与加速度约束联立：

$$
\boxed{
\begin{bmatrix}M&-J_c^T\\J_c&0\end{bmatrix}
\begin{bmatrix}\ddot q\\\lambda\end{bmatrix}
=
\begin{bmatrix}S^T\tau-h\\-\dot J_c\dot q\end{bmatrix}.}
$$

这是一种增广/KKT 形式。消去 $\ddot q$：

$$
J_cM^{-1}J_c^T\lambda
=J_cM^{-1}(h-S^T\tau)-\dot J_c\dot q.
$$

当 $M$ 正定且 $J_c$ 满行秩时：

$$
\boxed{
\lambda=(J_cM^{-1}J_c^T)^{-1}
\left[J_cM^{-1}(h-S^T\tau)-\dot J_c\dot q\right].}
$$

$\Lambda_c=(J_cM^{-1}J_c^T)^{-1}$ 可理解为接触空间等效惯性矩阵。

**课堂例题：**

$$
M=\begin{bmatrix}2&0\\0&1\end{bmatrix},\quad
J_c=\begin{bmatrix}1&1\end{bmatrix},\quad
S^T\tau=\begin{bmatrix}0\\3\end{bmatrix},\quad h=0,\quad\dot J_c=0.
$$

联立

$$
2\ddot q_1=\lambda,\quad\ddot q_2=3+\lambda,\quad
\ddot q_1+\ddot q_2=0
$$

得

$$
\boxed{\lambda=-2,\quad\ddot q_1=-1,\quad\ddot q_2=1.}
$$

此解满足等式约束，但若 $\lambda$ 被定义为单边接触的**正法向力**，则负值不具有物理可行性。抽象的 $J_c=[1,1]$ 本身并不能在没有法向约定时直接解释为地面接触。

## 6. 单边接触、互补条件与摩擦锥

对普通非粘附的刚性地面接触，定义间隙 $\phi(q)$ 与法向力 $f_n$：

$$
\boxed{\phi(q)\ge0,\qquad f_n\ge0,\qquad\phi(q)f_n=0.}
$$

- $\phi>0$：分离，$f_n=0$。
- $\phi=0$：可能存在非负法向接触力。
- $f_n<0$：意味着地面必须拉住物体，普通单边接触无法实现；通常应重新判断接触模式。

库仑静摩擦的可行条件：

$$
\boxed{|f_t|\le\mu f_n\quad\text{（二维）}},
\qquad
\boxed{\sqrt{f_x^2+f_y^2}\le\mu f_z,\ f_z\ge0\quad\text{（三维水平地面）}}.
$$

对于 $\mu=0.5$：

| 接触力 $(f_x,f_z)$ / N | 摩擦与法向可行性 | 理由 |
|---|---|---|
| $(20,100)$ | 可行 | $20\le50$ |
| $(60,100)$ | 不可行 | $60>50$，超出摩擦锥 |
| $(0,-10)$ | 不可行 | 法向力为负 |

这些条件是**局部接触力可行性条件**，并不单独保证整机静止或整体动力学平衡。

> 进一步说明：固定接触的等式动力学可以给出负法向力；这不是说真实地面产生了负摩擦力，而是该假定接触模式不成立。实际离地后，法向接触力和摩擦力均为零。

## 7. 本日核心公式速查

1. 固定接触位置：$p_c(q)=p_c^*$。
2. 接触速度约束：$J_c\dot q=0$。
3. 接触加速度约束：$J_c\ddot q+\dot J_c\dot q=0$。
4. 接触广义力：$Q_c=J_c^T\lambda$。
5. 约束动力学：$M\ddot q+h=S^T\tau+J_c^T\lambda$。
6. 单边接触：$\phi\ge0,\ f_n\ge0,\ \phi f_n=0$。
7. 静摩擦：$\|f_t\|\le\mu f_n$。

## 8. 易错点与 Day 3 衔接

- **接触点固定 $\ne$ 机器人整体固定**：允许的瞬时速度位于 $\ker J_c$。
- **不能遗漏 $\dot J_c\dot q$**：否则可能算出违反接触运动学的加速度。
- **接触力不等于可自由指定的控制输入**：需要与动力学和约束联立求解。
- **等式约束的解不一定物理可行**：还需检查法向力、摩擦锥以及是否发生离地或滑动。
- **角度约定要一致**：单连杆例子的 $\theta_b+q_1=0$ 对应水平向右。

**Day 2 结论：理论验收通过。**  
**Day 3 预告：**使用 MuJoCo 观察接触力、摩擦系数、滑动与离地，联系理论中的约束方程、互补条件和摩擦锥。当前尚未进行代码实验。
