# ROBO_Learn · Day 5 学习总结

**主题：双足机器人平衡、CoP/ZMP、LIPM、DCM 与 Capture Point**  
**进度：Day 5 完成；下一步 Day 6——离散时间 LIPM 与 MPC**

## 1. 平衡与支撑区域

- **支撑多边形（Support Polygon）**：平面接触点（或接触区域）的凸包。静态情况下，质心的竖直投影落在支撑区域内是维持平衡的重要必要条件；实际可行性还受摩擦、关节力矩等约束影响。
- **动态平衡**不能只看质心投影：还要考虑质心速度、加速度和角动量。
- **CoP（Center of Pressure）**：地面法向压力分布的加权中心；在离散竖直压力模型下，
  $$x_{\mathrm{CoP}}=\frac{\sum_i x_iF_{z,i}}{\sum_iF_{z,i}},\qquad y_{\mathrm{CoP}}=\frac{\sum_i y_iF_{z,i}}{\sum_iF_{z,i}}.$$

## 2. ZMP 的定义与两条计算路径

**ZMP（Zero Moment Point）**：在选定支撑平面上，使接触合力矩沿平面切向的两个分量为零的点。标准水平地面上，要求 $\tau_{z,x}=\tau_{z,y}=0$；绕竖直轴的力矩可以不为零。

### 路径 A：由地面接触力系求 ZMP

对于水平支撑面 $z=0$，以平面上的固定原点 $O$ 为参考：

$$f_r=\sum_i f_i,\qquad M_O=\sum_i p_i\times f_i+\sum_i\tau_i.$$

当 $F_z=(f_r)_z\ne0$ 时，

$$\boxed{x_z=-\frac{M_{O,y}}{F_z},\qquad y_z=\frac{M_{O,x}}{F_z}}.$$

推导依据是力矩换点：

$$M_{p_z}=M_O-p_z\times f_r.$$

**存在性与可行性**：$F_z\ne0$ 时，在给定水平平面上有唯一的数学 ZMP；$F_z=0$ 时可能不存在或不唯一。数学上算出的 ZMP 不一定在支撑区域内；但在纯压缩、共面接触且无额外切向接触力偶的标准条件下，实际 CoP/ZMP 位于接触区域凸包内。所需 ZMP 超出支撑区域意味着当前接触模式不能按该模型实现。

### 路径 B：由质心运动与角动量求 ZMP

设质心位置 $p_c=(x_c,y_c,z_c)$，质心角动量为 $L_c$，重力沿 $-z$ 方向。质心线动量与角动量定理：

$$m\ddot p_c=f_r-mg e_z,$$

$$\boxed{\dot L_c=(p_z-p_c)\times f_r+\tau_z}.$$

**角动量方程的来源**：对所有接触点 $p_i$，

$$\dot L_c=\sum_i(p_i-p_c)\times f_i+\sum_i\tau_i.$$

插入 $p_z$：$p_i-p_c=(p_i-p_z)+(p_z-p_c)$，并定义 $\tau_z=\sum_i(p_i-p_z)\times f_i+\sum_i\tau_i$，即得上式。质心处的均匀重力不产生力矩。该换点恒等式对任意参考点成立，选为 ZMP 时 $\tau_{z,x}=\tau_{z,y}=0$。

联立后得到：

$$\boxed{x_z=x_c-\frac{z_c\ddot x_c+\dot L_{c,y}/m}{g+\ddot z_c}}$$

$$\boxed{y_z=y_c-\frac{z_c\ddot y_c-\dot L_{c,x}/m}{g+\ddot z_c}}.$$

**理解**：路径 A 使用接触力和接触力矩；路径 B 使用质心加速度及角动量变化率。它们是同一 ZMP 定义在不同可用信息下的等价表达，并非两种 ZMP。

## 3. LIPM：线性倒立摆模型

LIPM（Linear Inverted Pendulum Model）的核心假设：

1. 支撑平面水平且接触条件可维持；
2. 质心高度固定 $z_c=h>0$，故 $\ddot z_c=0$；
3. 忽略质心角动量的水平变化率，$\dot L_{c,x}=\dot L_{c,y}=0$；
4. 将 ZMP 视为可在支撑区域内调节的**虚拟控制输入**，不显式建模关节、接触力分配与驱动约束。

于是二维水平方向解耦：

$$\boxed{\ddot x=\omega_0^2(x-u),\qquad \omega_0=\sqrt{g/h},\quad u=x_{\mathrm{ZMP}}}$$

$$\ddot y=\omega_0^2(y-y_{\mathrm{ZMP}}).$$

状态空间模型：

$$X=\begin{bmatrix}x\\\dot x\end{bmatrix},\qquad
\dot X=\underbrace{\begin{bmatrix}0&1\\\omega_0^2&0\end{bmatrix}}_AX+
\underbrace{\begin{bmatrix}0\\-\omega_0^2\end{bmatrix}}_Bu.$$

- 特征值 $\lambda=\pm\omega_0$：**开环不稳定**。
- 可控性矩阵 $[B,AB]$ 满秩：理想无约束模型可控。
- 若 $u=0$，则 $x(t)=c_1e^{\omega_0t}+c_2e^{-\omega_0t}$，通常含发散项。

### 简单 ZMP 状态反馈

目标为 $x_d$、目标速度为零。令 $e=x-x_d$，指定

$$\ddot e+k_d\dot e+k_pe=0,\qquad k_p,k_d>0.$$

解得

$$\boxed{u=x+\frac{k_p}{\omega_0^2}(x-x_d)+\frac{k_d}{\omega_0^2}\dot x}.$$

注意：若 $x>x_d$ 且 $\dot x=0$，要让质心向左加速，ZMP 需要移到质心右侧。

**输入饱和的关键结论**：真实 ZMP 必须在支撑区域内。将上述反馈截断为 $u\in[u_{\min},u_{\max}]$ 后，不能直接沿用无约束系统的全局渐近稳定性结论；某些初始状态甚至根本无法在不迈步时恢复平衡。

## 4. ZMP 作为输入，实际如何实现？

ZMP **不是**电机直接施加的控制量。真实浮动基座机器人的动力学为

$$M(q)\ddot q+h(q,\dot q)=S^T\tau+J_c(q)^Tf_c.$$

物理链条是：**关节力矩 $\tau$ → 机器人动力学及接触约束 → 接触力/力矩 → 实际 ZMP**。

- 踝关节力矩能够改变足底前后压力分布，从而改变 CoP/ZMP；在简化准静态条件下，力矩变化量级满足 $|\Delta\tau_{\mathrm{ankle}}|\approx F_z|\Delta x_z|$，符号取决于坐标约定。
- 常见分层架构：**LIPM/MPC（规划期望 CoM、ZMP）→ WBC（动力学与接触约束下分配加速度/接触力）→ 关节力矩或阻抗控制**。
- 强化学习策略可能直接输出关节位置目标或力矩，不显式计算 ZMP，但仍受接触动力学与可行支撑区域约束。

## 5. DCM：不稳定模态分解

DCM 全称 **Divergent Component of Motion（发散运动分量）**，不是 Divergent Control Mode。

定义

$$\boxed{\xi=x+\frac{\dot x}{\omega_0}},\qquad
\eta=x-\frac{\dot x}{\omega_0}.$$

得到

$$\boxed{\dot\xi=\omega_0(\xi-u),\qquad\dot\eta=-\omega_0(\eta-u)}.$$

因此：

- $\xi$ 对应 $+\omega_0$ 不稳定模态；
- $\eta$ 对应 $-\omega_0$ 稳定模态；
- 固定 ZMP 时 $\xi(t)-u=(\xi(0)-u)e^{\omega_0t}$。

这就是标准的线性系统模态分解在机器人平衡控制中的应用。

## 6. Capture Point：捕获点与可恢复性

在理想 LIPM 中，若希望选择一个**固定** ZMP $u$，使质心最终渐近静止，必须消除不稳定模态：

$$\xi(0)=u.$$

于是得到 Capture Point（CP）：

$$\boxed{x_{\mathrm{CP}}=x+\frac{\dot x}{\omega_0}}.$$

在经典固定高度 LIPM 中，CP 与 DCM 的表达式相同；DCM 强调不稳定状态变量，CP 强调使系统最终静止的支撑位置。

对于固定支撑区间 $[u_{\min},u_{\max}]$，理想模型下不迈步渐近捕获的条件为

$$\boxed{u_{\min}\le x+\frac{\dot x}{\omega_0}\le u_{\max}}.$$

这是状态平面 $(x,\dot x)$ 上两条平行直线之间的可捕获区域（在上述理想模型及控制假设下）。若捕获点超出支撑区域，往往需要**迈步改变支撑区域**，或利用未纳入 LIPM 的其他机制。

### 本日练习与答案

给定 $h=1\,\mathrm m$、$g\approx10\,\mathrm{m/s^2}$、$x(0)=0$、支撑区间 $[-0.1,0.1]\,\mathrm m$，求不迈步时允许的初始水平速度范围。

$$\omega_0=\sqrt{10},\qquad -0.1\le\frac{\dot x(0)}{\sqrt{10}}\le0.1.$$

故

$$\boxed{-\frac{\sqrt{10}}{10}\le\dot x(0)\le\frac{\sqrt{10}}{10}\ \mathrm{m/s}}$$

即 $|\dot x(0)|\le0.316\,\mathrm{m/s}$。边界情形在理想模型中仍可渐近捕获。

## 7. 模型适用范围与扩展

**ZMP 与 LIPM 的假设不能混淆**：ZMP 是接触力系的力矩定义；LIPM 才额外假设固定质心高度、忽略角动量变化率。

| 模型 / 框架 | 主要区别 | 适用方向 |
|---|---|---|
| LIPM | 固定质心高度，忽略质心角动量变化 | 简化平衡、ZMP 规划、MPC |
| VHIP | 允许质心高度变化，系数可随时间变化 | 起伏运动、非恒高步态 |
| LIPM + Flywheel | 保留部分角动量调节效应 | 躯干/摆臂辅助平衡 |
| ALIP | 使用位置和关于支撑点的角动量等状态 | 欠驱动行走、落脚点控制 |
| SLIP | 腿部弹性与弹跳运动 | 跑步、跳跃 |
| HZD | 完整非线性动力学 + 虚拟约束 + 混合系统分析 | 周期步态与轨道稳定性 |

**HZD**：Hybrid Zero Dynamics（混合零动态）；**PZD**：Partial Zero Dynamics（部分零动态）。HZD 是非线性控制与步态稳定性分析框架，并非与 LIPM 完全同层次的简化模型。其基本思路是用虚拟约束 $y(q)=h_0(q)-h_d(\theta(q))$ 定义零动态流形，并结合触地状态重置与 Poincaré 映射分析周期轨道稳定性。

LIPM 在平坦地面、中低速、质心高度变化小、角动量变化较小时较实用；在深蹲、跳跃、剧烈摆臂、复杂地形等情况下需更丰富的模型或完整动力学。

## 8. 今日关键认识与下一步

1. **ZMP 是接触力系的结果，LIPM 中把它抽象为虚拟控制输入。**
2. **LIPM 是开环不稳定但无约束可控的二阶系统；支撑区域约束决定实际可恢复性。**
3. **DCM 是不稳定模态坐标；Capture Point 将质心位置和速度合成为一个直观的落脚/支撑位置指标。**
4. **真实机器人通过关节力矩、足底压力分布和全身控制实现期望 ZMP。**
5. **HZD 等方法提供另一条基于完整非线性动力学的步态控制路线。**

**Day 6 预告**：LIPM 精确离散化、离散状态空间模型、ZMP 约束、有限时域预测，以及 MPC 如何在约束下规划质心/ZMP 轨迹。
