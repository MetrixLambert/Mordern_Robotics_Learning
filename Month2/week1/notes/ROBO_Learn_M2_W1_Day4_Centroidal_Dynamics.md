# ROBO_Learn｜Month 2 · Week 1 · Day 4

# 质心动力学（Centroidal Dynamics）学习总结

> 本节重点：质心线动量、整体角动量及其推导、Centroidal Momentum Matrix（CMM）、双足接触力分配。附录专门整理刚体力学与向量运算常用公式，便于复习。
>
> 约定：所有向量默认在**同一个世界坐标系**表达；重力加速度向量记为 $\boldsymbol g=[0,0,-g]^T$，其中标量 $g>0$。

## 1. 为什么使用质心动力学

具有 $n$ 个驱动关节的浮动基座机器人，其广义速度 $v\in\mathbb R^{n+6}$，完整动力学为

$$
M(q)\dot v+h(q,v)=S^T\tau+J_c(q)^T\lambda.
$$

质心动力学不直接追踪每个关节的全部运动，而是关注机器人**整体线动量与整体角动量**。它降低了平衡、接触力规划问题的描述维度，但**不能独立决定所有关节运动**。

## 2. 质心、线动量与外力

机器人由 $N$ 个刚体构成，第 $i$ 个刚体质量 $m_i$、质心位置 $p_i$：

$$
m=\sum_{i=1}^N m_i,\qquad
p_G=\frac{1}{m}\sum_{i=1}^N m_i p_i.
$$

总线动量：

$$
l_G=\sum_i m_i\dot p_i=m\dot p_G.
$$

对时间求导，并利用内部相互作用力成对抵消：

$$
\boxed{\dot l_G=m\ddot p_G=\sum_{k=1}^{N_c}f_k+m\boldsymbol g.}
$$

- $f_k$ 是外界在第 $k$ 个接触点施加在机器人上的力。
- **内部关节力矩不能直接改变系统整体线动量**；它们可通过改变姿态、接触力间接影响整体运动。
- 完全腾空、忽略空气阻力时：$\ddot p_G=\boldsymbol g$。摆腿不能改变整体质心的自由落体轨迹。

**课堂例题**：$m=60\,\mathrm{kg}$，$g=10\,\mathrm{m/s^2}$，双脚竖直支撑力为 $350$ N 与 $250$ N，则 $\ddot z_G=(350+250-600)/60=0$。若要 $\ddot z_G=1\,\mathrm{m/s^2}$，则总支撑力应为 $660$ N。

## 3. 整体角动量：定义及完整推导

### 3.1 质点角动量

质量为 $m$ 的质点，位置 $p$、速度 $\dot p$，关于参考点 $O$（位置 $p_O$）的角动量为

$$
L_O=(p-p_O)\times m\dot p.
$$

叉乘反映了**力臂方向与线动量方向**共同决定的转动效应。

### 3.2 一个刚体：自转与轨道角动量

对第 $i$ 个刚体，令 $p_i$ 为其质心，$\rho$ 为内部质量微元相对该质心的位置，则

$$
p=p_i+\rho,\qquad
\dot p=\dot p_i+\omega_i\times\rho.
$$

关于机器人整体质心 $G$ 的角动量：

$$
L_G^{(i)}=\int_{\mathcal B_i}(p-p_G)\times\dot p\,dm.
$$

代入并展开：

$$
\begin{aligned}
L_G^{(i)}={}&\int(p_i-p_G)\times\dot p_i\,dm\\
&+\int(p_i-p_G)\times(\omega_i\times\rho)\,dm\\
&+\int\rho\times\dot p_i\,dm\\
&+\int\rho\times(\omega_i\times\rho)\,dm.
\end{aligned}
$$

根据质心定义 $\int\rho\,dm=0$，中间两项消失；第一项为 $(p_i-p_G)\times m_i\dot p_i$。第四项使用向量三重积：

$$
\rho\times(\omega_i\times\rho)
=(\rho^T\rho)\omega_i-\rho\rho^T\omega_i.
$$

定义该刚体**关于自身质心**的转动惯量矩阵（在世界系表达）：

$$
I_i=\int_{\mathcal B_i}\left[(\rho^T\rho)I_3-\rho\rho^T\right]dm.
$$

因此

$$
\boxed{L_G^{(i)}=I_i\omega_i+(p_i-p_G)\times m_i\dot p_i.}
$$

- $I_i\omega_i$：**自转角动量**。
- $(p_i-p_G)\times m_i\dot p_i$：**轨道角动量**（该刚体质心相对整体质心运动产生）。

将所有刚体相加：

$$
\boxed{L_G=\sum_{i=1}^N\left[I_i\omega_i+(p_i-p_G)\times m_i\dot p_i\right].}
$$

**易错点**：即使 $p_i=p_G$，该刚体仍可能因 $\omega_i\ne0$ 而具有 $I_i\omega_i\ne0$ 的角动量。

### 3.3 角动量变化率

对整个机器人求角动量导数，内部力矩相互抵消；选择随整体质心移动的参考点时，额外项 $-\dot p_G\times m\dot p_G=0$，所以

$$
\boxed{\dot L_G=\sum_{k=1}^{N_c}(p_k-p_G)\times f_k+\sum_{k=1}^{N_c}\tau_k.}
$$

这里 $\tau_k$ 是**外部接触力偶矩**（不是机器人内部关节驱动力矩）。理想点接触通常取 $\tau_k=0$。均匀重力场下，重力关于整体质心的合力矩为零。

**必须区分**：

| 量 | 含义 | 公式 |
|---|---|---|
| $L_G$ | 当前整体角动量 | 各刚体自转与轨道角动量之和 |
| $\dot L_G$ | 角动量变化率 | 外部合力矩 |

无外力矩意味着 $\dot L_G=0$，**不意味着** $L_G=0$。内部反作用飞轮可改变机身与飞轮各自角动量，但无法改变整体总角动量。

## 4. 六维质心动量与 CMM

定义六维整体动量（本笔记采用**线动量在上、角动量在下**的排序）：

$$
h_G=\begin{bmatrix}l_G\\L_G\end{bmatrix}\in\mathbb R^6.
$$

广义速度 $v=[v_b^T,\dot q_j^T]^T\in\mathbb R^{n+6}$，Centroidal Momentum Matrix 满足

$$
\boxed{h_G=A_G(q)v,\qquad A_G(q)\in\mathbb R^{6\times(n+6)}.}
$$

对第 $i$ 个刚体的质心线速度和角速度：

$$
\dot p_i=J_{v,i}(q)v,\qquad \omega_i=J_{\omega,i}(q)v.
$$

利用叉乘矩阵 $[a]_\times b=a\times b$：

$$
\boxed{
A_G(q)=\begin{bmatrix}
\displaystyle\sum_i m_iJ_{v,i}\\[5pt]
\displaystyle\sum_i\left(I_iJ_{\omega,i}+m_i[p_i-p_G]_\times J_{v,i}\right)
\end{bmatrix}.}
$$

特别地，若 $\dot p_G=J_G(q)v$，则 CMM 的线动量部分满足

$$
A_G^{\mathrm{lin}}(q)=mJ_G(q).
$$

求导：

$$
\dot h_G=A_G(q)\dot v+\dot A_G(q,v)v
=\begin{bmatrix}
\sum_k f_k+m\boldsymbol g\\
\sum_k(p_k-p_G)\times f_k+\sum_k\tau_k
\end{bmatrix}.
$$

**直觉**：Jacobian 把广义速度映射到几何速度；CMM 把广义速度映射到**整体动量**。这是 WBC、动量控制与接触力优化的基础。

**一维课堂练习**：两个质点 $m_1=2$、$m_2=3$，$v=[\dot x_1,\dot x_2]^T$，则 $A=[2,3]$。当 $\dot x_1=1,\dot x_2=-1$ 时，总动量 $l_x=-1\,\mathrm{kg\,m/s}$，总质心速度 $\dot x_G=-1/5=-0.2\,\mathrm{m/s}$。

## 5. 双足接触力分配与可行性

二维水平地面、仅考虑左右脚的**竖直接触力**，忽略水平力与外部接触力偶矩。静态平衡要求

$$
\begin{cases}
f_{L,z}+f_{R,z}=mg,\\
(x_L-x_G)f_{L,z}+(x_R-x_G)f_{R,z}=0.
\end{cases}
$$

若 $x_L\ne x_R$：

$$
\boxed{
\begin{aligned}
f_{L,z}&=mg\frac{x_R-x_G}{x_R-x_L},\\
f_{R,z}&=mg\frac{x_G-x_L}{x_R-x_L}.
\end{aligned}}
$$

普通地面单边接触要求

$$
f_{L,z}\ge0,\qquad f_{R,z}\ge0.
$$

因此，在上述简化假设下，质心投影需位于 $[x_L,x_R]$ 内，才能找到满足非负法向力的静态分配。

**课堂例题**：$m=60$ kg，$g=10\,\mathrm{m/s^2}$，$x_L=-0.2$ m、$x_R=0.2$ m：

| 质心水平位置 $x_G$ | 左脚 $f_{L,z}$ | 右脚 $f_{R,z}$ | 是否可行 |
|---|---:|---:|---|
| $0$ m | $300$ N | $300$ N | 是 |
| $0.1$ m | $150$ N | $450$ N | 是 |
| $0.3$ m | $-150$ N | $750$ N | 否（左脚需拉住地面） |

这展示了：**满足动力学等式不代表满足接触物理约束**。注意，真实机器人还可能存在脚掌面接触力矩、水平摩擦力、外部支撑等因素；上述结论仅针对本题简化模型。

---

# 附录 A｜向量与刚体力学常用公式速查

> 这一部分是专门补充的复习资料，建议做题时随时翻阅。

## A1. 向量点积与叉积

设 $a=[a_x,a_y,a_z]^T$，$b=[b_x,b_y,b_z]^T$。

**点积**：

$$
a\cdot b=a^Tb=a_xb_x+a_yb_y+a_zb_z=\|a\|\|b\|\cos\theta.
$$

**叉积**：

$$
a\times b=\begin{bmatrix}
a_yb_z-a_zb_y\\a_zb_x-a_xb_z\\a_xb_y-a_yb_x
\end{bmatrix},\qquad \|a\times b\|=\|a\|\|b\|\sin\theta.
$$

叉积方向遵循右手定则；结果垂直于 $a,b$。

**常用计算律**：

$$
\begin{aligned}
a\times b&=-b\times a &&\text{（反交换律）},\\
a\times(b+c)&=a\times b+a\times c &&\text{（分配律）},\\
(\alpha a)\times b&=\alpha(a\times b),\\
a\times a&=0,\\
a\times(b\times c)&=b(a\cdot c)-c(a\cdot b),\\
(a\times b)\times c&=b(a\cdot c)-a(b\cdot c).
\end{aligned}
$$

**注意**：叉积**没有交换律，也没有结合律**，通常 $a\times(b\times c)\ne(a\times b)\times c$。

**标量三重积**：

$$
a\cdot(b\times c)=b\cdot(c\times a)=c\cdot(a\times b).
$$

## A2. 反对称矩阵（skew matrix）

$$
[a]_\times=\begin{bmatrix}
0&-a_z&a_y\\
a_z&0&-a_x\\
-a_y&a_x&0
\end{bmatrix},\qquad [a]_\times b=a\times b.
$$

常用恒等式：

$$
\begin{aligned}
[a]_\times^T&=-[a]_\times,\\
[a]_\times[b]_\times&=ba^T-(a^Tb)I_3,\\
[a]_\times^2&=aa^T-\|a\|^2I_3.
\end{aligned}
$$

由此可得：

$$
\rho\times(\omega\times\rho)
=-[\rho]_\times^2\omega
=\big(\|\rho\|^2I_3-\rho\rho^T\big)\omega.
$$

这是推导刚体转动惯量矩阵的关键一步。

## A3. 力矩、角动量与惯量

**力对参考点 $O$ 的力矩**：

$$
\tau_O=(p_f-p_O)\times f.
$$

**质点角动量**：

$$
L_O=(p-p_O)\times m\dot p.
$$

**刚体关于自身质心的角动量**：

$$
L_C=I_C\omega.
$$

**惯量张量定义**：

$$
I_C=\int\big(\|\rho\|^2I_3-\rho\rho^T\big)dm.
$$

**惯量张量坐标变换**（$R$ 把刚体坐标系向量映射到世界坐标系）：

$$
I_{\mathrm{world}}=R I_{\mathrm{body}}R^T.
$$

**平行轴定理**（参考点 $O$ 与质心 $C$ 的位移为 $d=p_C-p_O$）：

$$
I_O=I_C+m\big(\|d\|^2I_3-dd^T\big).
$$

**单刚体关于质心的旋转动力学**（世界系表达）：

$$
\tau_C=\frac{d}{dt}(I_{\mathrm{world}}\omega)
=I_{\mathrm{world}}\dot\omega+\omega\times(I_{\mathrm{world}}\omega).
$$

若在刚体坐标系表达，则是熟悉的欧拉方程

$$
\tau_{\mathrm{body}}=I_{\mathrm{body}}\dot\omega_{\mathrm{body}}+\omega_{\mathrm{body}}\times(I_{\mathrm{body}}\omega_{\mathrm{body}}).
$$

**注意**：不能一般地把 $\tau=I\dot\omega$ 当作三维刚体的完整旋转动力学；还存在陀螺耦合项 $\omega\times(I\omega)$。

## A4. 刚体上点的速度与加速度

若刚体上点 $P$ 相对于质心 $C$ 的位置为 $\rho=p_P-p_C$：

$$
\boxed{v_P=v_C+\omega\times\rho.}
$$

$$
\boxed{a_P=a_C+\dot\omega\times\rho+\omega\times(\omega\times\rho).}
$$

其中 $\dot\omega\times\rho$ 是切向加速度项，$\omega\times(\omega\times\rho)$ 是向心加速度项。

## A5. 对运动参考点求角动量导数

定义关于运动点 $G$ 的角动量

$$
L_G=\sum_a(p_a-p_G)\times m_a\dot p_a.
$$

一般有

$$
\dot L_G=\sum\tau_G^{\mathrm{ext}}-\dot p_G\times P,
\qquad P=\sum_a m_a\dot p_a,
$$

其中式中的 $\sum\tau_G^{\mathrm{ext}}$ 是全部外力关于 $G$ 的合力矩。当 $G$ 为**系统整体质心**时，$P=m\dot p_G$，于是额外项为零，得到 $\dot L_G=\sum\tau_G^{\mathrm{ext}}$。

## A6. 维度与单位检查

| 物理量 | 单位 |
|---|---|
| 质量 $m$ | kg |
| 速度 $v$ | m/s |
| 线动量 $l$ | kg·m/s |
| 力 $f$ | N = kg·m/s² |
| 力矩 $\tau$ | N·m |
| 转动惯量 $I$ | kg·m² |
| 角速度 $\omega$ | rad/s（量纲 s⁻¹） |
| 角动量 $L$ | kg·m²/s |
| 角动量变化率 $\dot L$ | N·m |

---

## 本节知识自检

1. 为什么内部关节力矩不能直接改变整体质心加速度？
2. 为什么一个刚体的质心位于机器人整体质心时，该刚体的角动量仍可能非零？
3. $L_G$ 与 $\dot L_G$ 有什么区别？
4. 为什么 $A_G^{\mathrm{lin}}=mJ_G$？
5. 为什么双足静态平衡需要同时检查等式和法向力非负约束？

**下一节**：Day 4.5 — MuJoCo 编程实战：从零定义 MJCF、加载模型、仿真循环、状态读取、控制器与可视化。
