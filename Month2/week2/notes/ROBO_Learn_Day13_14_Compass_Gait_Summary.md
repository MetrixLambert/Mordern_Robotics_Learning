# ROBO_Learn · Month 2 · Day 13–14 总结

> 主题：Compass Gait 的单支撑动力学、碰撞与换腿、混合动力学、Poincaré Map 和周期步态搜索。模型及实验以 `run_13_2.py`、`run_14_1.py`–`run_14_4.py` 为准。

## 1. 模型和符号约定

- 平面两腿等长 `l`；髋部集中质量 `M_h`，每条腿中点集中质量 `m`；腿杆自身转动惯量忽略。
- 世界 `x` 向右、`y` 向上；地面向右下倾斜 `γ>0`，方程 `y_ground=-x tan γ`。
- `q1` 为**支撑脚→髋部**相对世界竖直向上的绝对角度，`q2` 为**摆动脚→髋部**相对世界竖直向上的绝对角度；向右倾斜为正。
- 单支撑时支撑脚固定在原点，`q=[q1,q2]^T`。坡度不显式进入单支撑的 `M,C,G`，但进入足端触地几何。

位置：

\[
 p_h=l\begin{bmatrix}\sin q_1\\\cos q_1\end{bmatrix},\quad
 c_1=\frac l2\begin{bmatrix}\sin q_1\\\cos q_1\end{bmatrix},\quad
 c_2=l\begin{bmatrix}\sin q_1\\\cos q_1\end{bmatrix}-\frac l2\begin{bmatrix}\sin q_2\\\cos q_2\end{bmatrix}.
\]

摆动脚位置和 Jacobian：

\[
 p_{sw}=l\begin{bmatrix}\sin q_1-\sin q_2\\\cos q_1-\cos q_2\end{bmatrix},\qquad
 J_{sw}=l\begin{bmatrix}\cos q_1&-\cos q_2\\-\sin q_1&\sin q_2\end{bmatrix}.
\]

## 2. Day 13：单支撑 Euler–Lagrange 动力学

动能、势能：

\[
 T=\frac{l^2}{2}\left[\left(M_h+\frac54m\right)\dot q_1^2+\frac14m\dot q_2^2-m\cos(q_1-q_2)\dot q_1\dot q_2\right],
\]
\[
 V=gl\left[\left(M_h+\frac32m\right)\cos q_1-\frac12m\cos q_2\right].
\]

记 `δ=q1-q2`，单支撑动力学为

\[
 M(q)\ddot q+C(q,\dot q)\dot q+G(q)=0,
\]

\[
 M(q)=l^2\begin{bmatrix}M_h+5m/4&-(m/2)\cos\delta\\-(m/2)\cos\delta&m/4\end{bmatrix},
\]
\[
 C(q,\dot q)\dot q=\frac{ml^2}{2}\begin{bmatrix}-\sin\delta\,\dot q_2^2\\\sin\delta\,\dot q_1^2\end{bmatrix},\quad
 G(q)=gl\begin{bmatrix}-(M_h+3m/2)\sin q_1\\(m/2)\sin q_2\end{bmatrix}.
\]

状态 `x=[q1,q2,qdot1,qdot2]^T`，数值积分用 `np.linalg.solve(M,-(Cqdot+G))`，不显式求逆。摆动脚竖直地面间隙：

\[
 h_{sw}=y_{sw}+x_{sw}\tan\gamma.
\]

触地事件为 `h_sw=0` 且从正向负穿越。理想单支撑阶段机械能 `E=T+V` 守恒。

**实验 13.2**：`l=1 m, M_h=10 kg, m=1 kg, γ=3°`，初态 `[15°,-15°,0.8,-0.5]`。第一次触地约 `t=0.105535 s`，机械能数值漂移约 `2.84e-14 J`。说明耦合动力学和保守性数值实现一致。

## 3. Day 14.1：碰撞冲量和支撑腿切换

理想瞬时碰撞的通用方程：

\[
 M(q)(\dot q^+-\dot q^-)=J_c^T\Lambda.
\]

若碰撞后接触点无滑动、完全非弹性：`J_c qdot^+=0`。由此得到

\[
 \Lambda=-(J_cM^{-1}J_c^T)^{-1}J_c\dot q^-,
\]
\[
 \dot q^+=\left[I-M^{-1}J_c^T(J_cM^{-1}J_c^T)^{-1}J_c\right]\dot q^-.
\]

**关键**：Day 13 的二维质量矩阵已隐含**旧支撑脚固定**，不能在旧脚释放、新脚接触的冲击瞬间直接使用。临时引入浮动基座 `qbar=[x_h,y_h,q1,q2]^T`、四维质量矩阵 `Mbar` 和新脚 Jacobian：

\[
 J_2=\begin{bmatrix}1&0&0&-l\cos q_2\\0&1&0&l\sin q_2\end{bmatrix}.
\]

在**旧脚无冲击冲量、立即释放**的模型假设下，求解

\[
 \begin{bmatrix}\bar M&-J_2^T\\J_2&0\end{bmatrix}
 \begin{bmatrix}\dot{\bar q}^+\\\Lambda\end{bmatrix}
 =\begin{bmatrix}\bar M\dot{\bar q}^-\\0\end{bmatrix}.
\]

然后进行腿角色交换 `q^+=R q^-`，`qdot^+=R [qdot1_old^+,qdot2_old^+]`，`R=[[0,1],[1,0]]`。**角色交换不是额外一次物理碰撞。**

物理可行性需另外检查：法向冲量非负、切向冲量满足摩擦锥、旧脚能够离地、碰撞不增能、触地前脚不穿透地面等。

**实验 14.1**：第一次碰撞水平/竖直冲量约 `(-2.131,+7.613) N·s`，碰撞后新脚速度近零，动能损失约 `2.704 J`。

## 4. Day 14.2：多步混合动力学

将连续流和离散重置交替：

\[
 \dot x=f(x)\quad(h_{sw}>0),\qquad x^+=\Delta(x^-)\quad(\text{有效触地}).
\]

**实验 14.2**：第 1 次碰撞有效；第 2 次候选碰撞在 `t≈0.302748 s` 时计算出法向冲量约 `-1.840 N·s`、旧脚离地速度约 `-0.1668 m/s`，因此按当前单向、无滑动、立即换腿模型判定无效并停止。这**不证明被动 Compass Gait 不可能持续下坡行走**；只是该初值和接触模式无法继续。

## 5. Day 14.3–14.4：Poincaré Map 与周期轨道

以**碰撞后、完成换腿**的状态作为截面。等长腿触地几何给出 `q1+q2=2γ`（非退化构型），因此选取三维截面坐标：

\[
 z=[q_1,\dot q_1,\dot q_2]^T,\quad q_2=2\gamma-q_1.
\]

一步映射 `z_{k+1}=P(z_k)` 包含连续积分、触地、冲量重置及换腿。周期步态是固定点 `P(z*)=z*`；局部稳定性由 `DP(z*)` 的特征值判断（谱半径小于 1 为局部渐近稳定的充分条件，通常针对光滑且有效的截面映射）。

- **14.3 事件驱动搜索**：使用 `scipy.optimize.root`，102 组初始猜测未找到候选固定点。数值搜索失败不等于不存在周期轨道。
- **14.4 Direct Shooting**：把 `y=[q1(0),qdot1(0),qdot2(0),T]` 作为优化变量；使用 `least_squares` 求解 3 个周期性残差 + 1 个终端触地间隙残差；再检查整个步内间隙、冲量方向、摩擦、离地及碰撞耗能。
- **14.4 数值结果**：`γ=7°` 时找到当前模型假设下可行的周期轨道：`z*≈[-0.24543210,1.42611425,0.17759131]^T`，`T≈0.556 s`，固定点残差约 `5.66e-16`。冲量法向约 `+13.8145 N·s`，碰撞损失约 `10.3117 J`。
- 数值线性化得到 Poincaré 特征值约 `[9.3467,0.00953,0.33178]`，谱半径 `≈9.35>1`：**找到周期轨道，但它局部不稳定**。

## 6. 本阶段结论与待办

1. Compass Gait 与 Rimless Wheel 都通过重力势能下降补偿非弹性碰撞耗能，可能存在被动周期步态。
2. 单次碰撞有效 ≠ 多步有效；周期轨道存在 ≠ 周期轨道稳定；局部稳定 ≠ 大范围鲁棒。
3. 碰撞方程和 Day 11 的广义冲量—动量方程相同；新增困难是**接触约束切换**与**腿角色交换**。
4. `run_14_4.py` 的可行性检查是针对理想点接触、给定摩擦系数和采样间隙的数值检验，并非完整的接触互补性证明。
5. **选做、暂缓**：搜索稳定被动周期轨道，扩大参数范围或使用延拓法，并通过扰动多步仿真验证吸引域。
6. **下一步 Day 15**：给 Compass Gait 引入髋关节力矩，研究欠驱动系统的驱动矩阵、可控与欠驱动的区别，以及单支撑阶段的主动控制。

## 7. 对应脚本

- `run_13_2.py`：单支撑动力学与能量守恒。
- `run_14_1.py`：浮动基座碰撞冲量与腿角色交换。
- `run_14_2.py`：多步混合仿真与无效接触检测。
- `run_14_3.py`：事件驱动 Poincaré Map 和固定点搜索。
- `run_14_4.py`：Direct Shooting、可行性筛选、周期轨道及局部稳定性。

> 备注：Day 13/14 使用的绝对角度定义应始终与模型脚本保持一致；髋部力矩的输入矩阵符号也必须依据相同坐标定义推导。
