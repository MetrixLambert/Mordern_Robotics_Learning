# ROBO_Learn · Month 2 · Week 1 · Day 1

**主题：浮动基座与欠驱动系统（Floating-Base & Underactuated Systems）**  
**学习状态：理论学习与问答验收完成；本单元未安排独立代码实验**  
**主教材：** Russ Tedrake, *Underactuated Robotics*, Chapter 1（尤其 §1.2 Definitions、§1.3 Feedback Equivalence）  
**后续衔接：** Chapter 5（腿足机器人）、Chapter 17（接触与混合系统）、Appendix B（多体动力学）  
**辅助参考：** Lynch & Park, *Modern Robotics*（Jacobian、动力学）；Murray, Li & Sastry, *A Mathematical Introduction to Robotic Manipulation*（刚体与约束建模）

## 1. 与 Month 1 的衔接

固定基座、完全驱动机械臂通常采用

$$
M(q)\ddot q+h(q,\dot q)=\tau.
$$

在模型准确、力矩不受限等理想条件下，Computed Torque Control 令 $\tau=M(q)v+h(q,\dot q)$，可实现 $\ddot q=v$。

浮动基座机器人额外具有整体平移和旋转自由度，而关节电机不能直接对这些基座自由度施加广义力。因此，无法简单地将全部广义加速度独立指定。接触力和自然动力学成为控制问题的核心。

## 2. 欠驱动的定义与判据

考虑二阶系统

$$
\ddot q=f(q,\dot q,u,t),\qquad q\in\mathbb R^n,\quad u\in\mathbb R^m.
$$

在给定状态和时间，若对任意期望加速度 $\ddot q_d\in\mathbb R^n$，都存在允许输入 $u$ 使 $f(q,\dot q,u,t)=\ddot q_d$，则系统在该状态下**完全驱动**；否则为**欠驱动**。

对于输入不受约束、质量矩阵可逆的机械系统

$$
M(q)\ddot q+h(q,\dot q)=B(q)u,
$$

有

$$
\ddot q=-M^{-1}h+M^{-1}Bu,
\qquad \operatorname{rank}(M^{-1}B)=\operatorname{rank}(B).
$$

因此，**在上述假设下**，$\operatorname{rank}(B)<n$ 表示欠驱动；$\operatorname{rank}(B)=n$ 表示完全驱动。输入饱和等约束也可能使加速度无法任意指定。

**区别：** 欠驱动描述瞬时加速度的独立指定能力；可控性描述在一段时间内是否可以通过输入实现状态转移。欠驱动不必然意味着不可控。

## 3. 练习一：惯性耦合不一定带来可控性

给定

$$
M=\begin{bmatrix}2&1\\1&2\end{bmatrix},\quad B=\begin{bmatrix}0\\1\end{bmatrix},\quad M\ddot q=Bu.
$$

输入矩阵秩为 1，小于广义坐标维度 2，因此欠驱动。求逆得到

$$
\ddot q=\begin{bmatrix}-1/3\\2/3\end{bmatrix}u.
$$

当 $u=3$，$\ddot q=[-1,2]^T$。第一个自由度虽然没有电机，仍受惯性耦合影响。

令 $x=[q_1,q_2,\dot q_1,\dot q_2]^T$，

$$
A=\begin{bmatrix}0&I_2\\0&0\end{bmatrix},\qquad
\mathcal B=\begin{bmatrix}0\\0\\-1/3\\2/3\end{bmatrix}.
$$

由于 $A^2=0$，$\operatorname{rank}[\mathcal B,A\mathcal B,A^2\mathcal B,A^3\mathcal B]=2<4$，系统**不可控**。其加速度始终满足 $\ddot q_1=-\frac12\ddot q_2$。

## 4. 练习二：欠驱动但可控的二阶系统

考虑抽象系统

$$
\ddot q_1=u,\qquad \ddot q_2=q_1.
$$

输入只有一维，不能任意指定二维加速度，故欠驱动。采用 $x=[q_1,q_2,\dot q_1,\dot q_2]^T$：

$$
A=\begin{bmatrix}
0&0&1&0\\
0&0&0&1\\
0&0&0&0\\
1&0&0&0
\end{bmatrix},\qquad
\mathcal B=\begin{bmatrix}0\\0\\1\\0\end{bmatrix}.
$$

可控性矩阵为

$$
\mathcal C=[\mathcal B,A\mathcal B,A^2\mathcal B,A^3\mathcal B]
=\begin{bmatrix}
0&1&0&0\\
0&0&0&1\\
1&0&0&0\\
0&0&1&0
\end{bmatrix},
$$

故 $\operatorname{rank}(\mathcal C)=4$，**系统可控**。状态耦合 $\ddot q_2=q_1$ 使输入通过动态演化影响第二自由度。

**本次易错点：** 初次构造 $A$ 时遗漏了 $A_{41}=1$；状态方程的耦合项不能漏掉。该例是用于区分概念的抽象二阶系统，不作为典型自由漂浮机械系统的物理模型。

## 5. 浮动基座的构型与广义速度

固定基座机械臂仅需关节坐标 $q_j$；浮动基座还需描述躯干/骨盆相对于世界系的位姿 $T_{WB}\in SE(3)$。

三维浮动基座机器人具有 6 个基座自由度。对于 $n$ 个独立关节，构型空间局部维度为 $n+6$；可以在局部使用

$$
q=\begin{bmatrix}q_b\\q_j\end{bmatrix},\quad q_b\in\mathbb R^6,\ q_j\in\mathbb R^n.
$$

这里 $q_b$ 可由世界系下基座位置与局部姿态参数组成。由于三维角速度一般不等于姿态参数的导数，后续动力学宜使用

$$
v=\begin{bmatrix}v_b\\\dot q_j\end{bmatrix}\in\mathbb R^{n+6},\qquad \dot q=N(q)v.
$$

$v_b$ 的线速度、角速度必须声明参考系和分量顺序。对平面基座 $(x_b,z_b,\theta_b)$，可方便地取 $v=\dot q$。

## 6. 驱动选择矩阵

无环境接触、无基座推进器且全部 $n$ 个关节独立驱动时，

$$
M(q)\dot v+h(q,v)=S^T\tau,
\quad S=\begin{bmatrix}0_{n\times6}&I_n\end{bmatrix},
\quad S^T\tau=\begin{bmatrix}0_{6\times1}\\\tau\end{bmatrix}.
$$

其中 $M\in\mathbb R^{(n+6)\times(n+6)}$，$h\in\mathbb R^{n+6}$，$S\in\mathbb R^{n\times(n+6)}$，$\tau\in\mathbb R^n$。

由于 $\operatorname{rank}(S^T)=n<n+6$，该模型欠驱动。**前六项为零仅表示电机没有直接施加基座广义力，不代表基座加速度为零，也不代表没有重力。**

**平面练习：** 基座 $(x_b,z_b,\theta_b)$ 加两个驱动关节，总自由度为 5，浮动基座自由度为 3，

$$
S=\begin{bmatrix}0_{2\times3}&I_2\end{bmatrix},\qquad S^T\tau=\begin{bmatrix}0_{3\times1}\\\tau\end{bmatrix}.
$$

## 7. 浮动基座动力学的分块与 Schur 补

对平面机器人取 $q=[q_b^T,q_j^T]^T$、$v=\dot q$，且无接触：

$$
\begin{bmatrix}M_{bb}&M_{bj}\\M_{jb}&M_{jj}\end{bmatrix}
\begin{bmatrix}\ddot q_b\\\ddot q_j\end{bmatrix}
+\begin{bmatrix}h_b\\h_j\end{bmatrix}
=\begin{bmatrix}0\\\tau\end{bmatrix}.
$$

其中 $M_{bb}\in\mathbb R^{3\times3}$、$M_{bj}\in\mathbb R^{3\times2}$、$M_{jb}=M_{bj}^T$、$M_{jj}\in\mathbb R^{2\times2}$。

第一行给出

$$
\ddot q_b=-M_{bb}^{-1}(M_{bj}\ddot q_j+h_b).
$$

代入关节方程，得到

$$
\underbrace{(M_{jj}-M_{jb}M_{bb}^{-1}M_{bj})}_{\bar M}\ddot q_j
+h_j-M_{jb}M_{bb}^{-1}h_b=\tau.
$$

$\bar M$ 是关于 $M_{bb}$ 的 **Schur 补**，表示消去自由基座加速度后的等效关节惯性。正定质量矩阵保证 $M_{bb}$ 可逆且 $\bar M$ 正定。

### 数值验收：固定基座与浮动基座对比

$$
M=\begin{bmatrix}2&1\\1&2\end{bmatrix},\quad S=[0\ \ 1],\quad M\ddot q=S^T\tau.
$$

第一行 $2\ddot q_1+\ddot q_2=0$，所以 $\ddot q_1=-\frac12\ddot q_2$。等效惯性

$$
\bar M=2-1\cdot\frac12\cdot1=\frac32,
\qquad \ddot q_2=\frac23\tau.
$$

如果**另行施加固定基座约束** $\ddot q_1=0$，基座约束反力负责平衡第一行，关节动力学则变为

$$
2\ddot q_2=\tau,\qquad \ddot q_2=\frac12\tau.
$$

同样的关节加速度目标 $\ddot q_2=1$：自由基座所需 $\tau=1.5$；固定基座所需 $\tau=2$。若按固定基座模型给真实浮动系统施加 $\tau=2$，实际关节加速度为 $4/3$，不是 1。

**物理结论：** 固定基座不是简单地忽略基座坐标，而是增加外部运动约束；其反力和等效惯性都会改变。

## 8. 与接触控制的衔接

机器人站在地面时，动力学需要增加接触广义力：

$$
\boxed{M(q)\dot v+h(q,v)=S^T\tau+J_c(q)^T\lambda.}
$$

$J_c$ 为接触 Jacobian，$\lambda$ 为接触力/接触 wrench（具体维度取决于接触模型）。接触力并非任意独立输入，还需满足接触运动学、法向非穿透和摩擦等条件。

这将连接后续的质心动力学、接触力分配、QP 和 Whole-Body Control。无外力矩时，空中摆腿可以改变躯干姿态，却不能任意改变系统总角动量；内部关节力也不能改变质心的自由飞行轨迹。

## 9. 本日关键结论与验收

1. **欠驱动 $\not\Rightarrow$ 不可控**：瞬时加速度映射与有限时间状态可达性是不同概念。
2. **浮动基座额外增加 6 个自由度（三维）**，即使每个关节都有电机，基座也通常没有直接驱动。
3. **零基座驱动力 $\not\Rightarrow$ 零基座加速度**：惯性耦合、重力和接触均会影响基座运动。
4. **浮动与固定基座的等效惯性不同**：可通过 Schur 补解释。
5. **接触会改变运动约束和广义力平衡**：不能仅将接触力当成自由控制输入。

**实际完成：** 欠驱动秩判据；两个可控性分析例题；平面浮动基座自由度与驱动矩阵；分块动力学及 Schur 补；固定/浮动基座数值对比。理论问答已验收，未进行本日独立代码实验。

**建议复习：** 闭卷写出 $S$、浮动基座分块动力学和 Schur 补；重新计算欠驱动但可控例子的 $A$ 与可控性矩阵。

**下一单元：Month 2 · Week 1 · Day 2** —— 接触点位置约束、接触 Jacobian、速度/加速度约束、Lagrange 乘子与接触力。
