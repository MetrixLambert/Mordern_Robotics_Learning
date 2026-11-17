# ROBO_Learn · Day 12 学习总结

**主题：Rimless Wheel（无轮缘轮）——被动步行、混合动力学与周期步态稳定性**  
**状态：已完成（实验 12.1、12.2）**

## 1. 本日学习目标

通过最简化的被动步行模型，掌握从**连续动力学 → 碰撞重置 → 一步返回映射（Poincaré Map）→ 周期轨道与局部稳定性**的分析流程。区别站立控制中的平衡点稳定性与步行中的周期轨道稳定性，并学习 `scipy.integrate.solve_ivp` 的事件驱动混合系统仿真。

## 2. 模型与坐标约定

- 轮心为集中质量 $m$，辐条无质量，长度 $l$；相邻辐条夹角 $2\alpha$。
- 斜坡倾角 $\gamma>0$，沿下坡方向前进。
- $\theta$ 是当前支撑辐条相对于**斜面法线**的有符号角度，正方向为前进转动；一步内 $\theta\in[-\alpha,+\alpha]$。
- 支撑点不滑动；旧支撑点在新辐条触地时释放；新支撑点处发生理想完全非弹性碰撞。
- 碰撞前后上标 $-$、$+$ 表示同一碰撞事件的前后状态。$\omega=\dot\theta$。

**高度投影的区别：** 质心相对于当前支撑点的世界竖直高度为

$$h=l\cos(\theta+\gamma),$$

而相对于斜面的法向距离为 $h_\perp=l\cos\theta$。**重力势能必须使用世界竖直高度**，不能直接用斜面法向距离。

## 3. 连续动力学（Flow）

动能、势能和拉格朗日函数：

$$T=\tfrac12ml^2\dot\theta^2,\qquad V=mgl\cos(\theta+\gamma),\qquad L=T-V.$$

在无驱动、无阻尼、理想固定支撑点假设下，重力已包含于势能，支撑约束力不做虚功，故 Euler–Lagrange 方程的非保守广义力 $Q_\theta=0$：

$$\frac{d}{dt}\frac{\partial L}{\partial\dot\theta}-\frac{\partial L}{\partial\theta}=0
\quad\Longrightarrow\quad
\boxed{\ddot\theta=\frac gl\sin(\theta+\gamma)}.$$

定义状态 $x=[\theta,\omega]^\top$，则

$$\boxed{\dot x=f(x)=\begin{bmatrix}\omega\\(g/l)\sin(\theta+\gamma)\end{bmatrix}}.$$

## 4. 碰撞 Guard 与 Reset Map

**Guard：** 下一辐条触地时 $\theta=+\alpha$，且 $\omega>0$；可定义事件函数

$$g(x)=\theta-\alpha,$$

仅检测从负到正的零点穿越。

**角度重置：** 碰撞瞬间轮子物理姿态连续，但新支撑辐条取代旧支撑辐条成为坐标参考，故

$$\boxed{\theta^+=-\alpha=\theta^--2\alpha}.$$

这**不是**轮子在瞬间发生 $2\alpha$ 的真实转动。

**角速度重置：** 关于新支撑点 $P$，在理想瞬时碰撞中，新支撑点的冲量对 $P$ 的角冲量为零，重力等有限力的角冲量可忽略，并假定旧支撑点不施加额外冲量，故关于 $P$ 的角动量守恒：

$$H_P^-=H_P^+.$$

碰撞前轮心绕旧支撑点运动，关于新支撑点的角动量为 $ml^2\cos(2\alpha)\omega^-$；碰撞后绕新支撑点运动，为 $ml^2\omega^+$。因此

$$\boxed{\omega^+=\cos(2\alpha)\omega^-}.$$

令 $c=\cos(2\alpha)$，则完整 Reset Map：

$$\boxed{\Delta(x^-)=\begin{bmatrix}-\alpha\\c\omega^-\end{bmatrix}}.$$

碰撞的动能保留比例为 $T^+/T^-=c^2$（本模型 $\alpha=15^\circ$ 时为 $0.75$），即每次碰撞损失约 25% 的碰撞前动能。**不是所有动量守恒**：地面冲量会改变线动量；守恒的是特定假设下关于新支撑点的角动量。

## 5. Poincaré Map：从一步到下一步

设 $\omega_k$ 为第 $k$ 次碰撞**后**的角速度。两次碰撞之间机械能守恒：

$$\frac12ml^2\omega_k^2+mgl\cos(\gamma-\alpha)
=\frac12ml^2(\omega^-)^2+mgl\cos(\gamma+\alpha).$$

利用 $\cos(\gamma-\alpha)-\cos(\gamma+\alpha)=2\sin\gamma\sin\alpha$，定义

$$K=\frac{4g}{l}\sin\gamma\sin\alpha,$$

得到 $\omega^- =\sqrt{\omega_k^2+K}$（选择向前运动的正速度分支）。结合碰撞重置：

$$\boxed{\omega_{k+1}=P(\omega_k)=c\sqrt{\omega_k^2+K}}.$$

该映射只在轮子确实能越过势能障碍并到达下一触地事件时适用，不能不加检查地用于任意初始速度。

### 5.1 周期步态不动点

周期步态要求 $P(\omega^*)=\omega^*$，于是

$$\boxed{\omega^*=\sqrt{\frac{c^2K}{1-c^2}}}.$$

此时每一步的重力势能下降恰好补偿碰撞能量损失。**周期步态不是静止平衡点**：完整状态在时间上周期变化，只是每一步碰撞后的截面状态相同。

### 5.2 局部稳定性

$$P'(\omega)=\frac{c\omega}{\sqrt{\omega^2+K}},\qquad
\boxed{P'(\omega^*)=c^2}.$$

在不动点附近，误差 $e_k=\omega_k-\omega^*$ 满足

$$e_{k+1}\approx c^2e_k.$$

若 $|c^2|<1$ 且周期步态可实现，则周期轨道在 Poincaré 截面上局部渐近稳定。多维双足模型中相应条件推广为一步映射 Jacobian 的特征值模均小于 1（或谱半径小于 1）。

### 5.3 势能障碍

对于 $0<\gamma<\alpha$，从 $\theta=-\alpha$ 到 $\theta=-\gamma$ 需先爬升至势能峰值，最低起步角速度阈值为

$$\boxed{\omega_{\min}=\sqrt{\frac{2g}{l}\,[1-\cos(\alpha-\gamma)]}}.$$

严格越过障碍需 $\omega_0>\omega_{\min}$；等号对应理想情况下到达势能顶点时速度为零。

## 6. 实验 12.1：直接迭代一步映射

参数：$g=9.81\ \mathrm{m/s^2}$，$l=1\ \mathrm m$，$\alpha=15^\circ$，$\gamma=5^\circ$。选择不同初始角速度（如 $0.65,0.95,1.40,2.00\ \mathrm{rad/s}$），重复迭代 $P(\omega)$。

**实验现象：** 所有可行初始速度的迭代轨迹都收敛到同一不动点；误差的半对数曲线近似为直线。

**关键结果：**

- 不动点 $\omega^*\approx1.62956334\ \mathrm{rad/s}$。
- 局部稳定性乘子 $c^2=0.75$。
- 势能障碍阈值 $\omega_{\min}\approx0.545960\ \mathrm{rad/s}$。

## 7. 实验 12.2：事件驱动的完整混合仿真

采用 `scipy.integrate.solve_ivp`：

1. 从 $x_0=[-\alpha,\omega_0]^\top$ 出发积分 $\dot x=f(x)$。
2. 事件函数 `touchdown(t, x) = x[0] - alpha`；设置 `touchdown.terminal = True`、`touchdown.direction = 1`。
3. 通过 `sol.t_events[0]` 和 `sol.y_events[0]` 读取触地时刻与碰撞前状态。
4. 执行 $x^+=[-\alpha,c\omega^-]^\top$，以重置后的状态重新启动积分器。
5. 逐步记录 $\theta(t)$、$\omega(t)$ 和碰撞后的 $\omega_k$，与解析映射对比。

**本次实际运行输出：**

```text
Fixed-point post-impact speed: 1.62956334 rad/s
Initial post-impact speed:     2.00000000 rad/s
Final post-impact speed:       1.63087107 rad/s
Stability multiplier:          0.75000000
Max ODE-vs-map discrepancy:    1.256e-11 rad/s
Simulated 20 steps; final time: 6.1843 s
```

**图像观察：**

- $\theta(t)$ 在 $-15^\circ$ 到 $+15^\circ$ 间重复，因参考辐条切换而呈锯齿形重置。
- $\omega(t)$ 在每一步连续演化，在碰撞瞬间向下跳变；步态逐渐趋于周期重复。
- $\omega_k$ 从 $2.0$ 收敛至约 $1.63\ \mathrm{rad/s}$。
- ODE 积分与解析 Poincaré Map 的最大差异约 $1.256\times10^{-11}\ \mathrm{rad/s}$，与数值积分精度一致。

## 8. 关键认识与后续衔接

1. **平衡点 vs. 极限环：** 站立控制研究 $x(t)\to x^*$；步行稳定性研究运动是否趋向重复的周期轨道。
2. **混合系统的三要素：** 连续动力学（Flow）、切换条件（Guard）、离散重置（Reset）。
3. **Poincaré 截面降维：** 把连续时间周期运动的稳定性转化为离散一步映射的不动点稳定性。
4. **解析映射的适用边界：** 只有在下一次触地可达、理想碰撞假设成立时才适用；实际接触和复杂双足模型通常需要数值积分与事件检测。
5. **程序实现重点：** `solve_ivp` 的 `events`、`terminal`、`direction`、`t_events`、`y_events`，以及碰撞后重新启动积分。

**Day 12 已完成。** 下一步进入 **Day 13：Compass Gait（罗盘步态双足机器人）**，将上述方法推广到具有摆动腿和支撑腿、多个广义坐标、腿角色交换及更复杂碰撞重置的双足步行模型。
