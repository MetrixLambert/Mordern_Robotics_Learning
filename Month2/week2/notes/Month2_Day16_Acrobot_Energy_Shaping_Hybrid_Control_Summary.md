# ROBO_Learn · Month 2 · Day 16 总结

## 主题：Acrobot Energy Shaping、Swing-up 与 LQR 分段控制

**主教材**：Russ Tedrake, *Underactuated Robotics*, Chapter 3（Energy-based swing-up、PFL）；衔接 Day 15 的 Acrobot 动力学、LQR 与 PFL。

**学习进度**：Day 16 Part 1–3 已完成；后续从 Day 17 进入高自由度腿足机器人动力学（教材 Chapter 5）。

## 1. 问题与目标

Acrobot 有两个广义坐标但仅肘关节驱动：

\[
q=[q_1,q_2]^\top,\qquad M(q)\ddot q+h(q,\dot q)=B\tau,\quad B=[0,1]^\top.
\]

- 自然悬垂平衡：\(q=(0,0)\)。
- 直立目标：\(x^*=(\pi,0,0,0)\)。
- Day 15 的 LQR 仅提供直立附近的局部稳定化，不能保证从悬垂位置起摆。
- Day 16 研究：先由非线性能量控制实现大范围起摆，再切换到 LQR 捕获直立平衡。

## 2. 机械能与功率平衡

总能量：

\[
E=\tfrac12\dot q^\top M(q)\dot q+U(q),\qquad
U=-(m_1l_{c1}+m_2l_1)g\cos q_1-m_2gl_{c2}\cos(q_1+q_2).
\]

对能量求导并代入动力学：

\[
\dot E=\dot q^\top B\tau-\dot q^\top C\dot q+\tfrac12\dot q^\top\dot M\dot q.
\]

标准机械系统的 Coriolis 矩阵满足 \(\dot M-2C\) 反对称，因此

\[
\dot q^\top(\dot M-2C)\dot q=0
\quad\Rightarrow\quad
\boxed{\dot E=\dot q^\top B\tau=\dot q_2\tau}.
\]

**要点**：不是 \(\dot M\) 项单独为零，而是它与 Coriolis/离心力项抵消。\(\tau\dot q_2\) 是肘关节电机的瞬时功率。

本次模型参数：\(m_1=m_2=l_1=I_1=I_2=1\)，\(l_{c1}=l_{c2}=0.5\)，\(g=9.81\)。下垂静止能量 \(-19.62\,\mathrm J\)，直立静止目标能量 \(E_d=19.62\,\mathrm J\)，所需增加的机械能为 \(39.24\,\mathrm J\)。

## 3. Direct Energy Feedback

令 \(\tilde E=E-E_d\)，选取 \(V_E=\frac12\tilde E^2\)。控制律：

\[
\boxed{\tau=-k_E\dot q_2\tilde E},\qquad
\dot V_E=-k_E\dot q_2^2\tilde E^2\le0.
\]

此控制器直接调节电机功率，能够使能量误差不增，但 \(\dot V_E\le0\) 本身不能保证 \(E\to E_d\)，更不能保证 \(x\to x^*\)。当 \(E\approx E_d\) 时，\(\tau\approx0\)，系统近似无驱动、无阻尼，动能和势能可以持续交换，关节可能继续摆动或旋转。完全静止的下垂平衡点也不会自行起摆。

## 4. 教材形式：Collocated PFL + Energy Shaping

先通过 Collocated PFL 使 \(\ddot q_2=v\)：

\[
\tau=\left(M_{22}-\frac{M_{21}M_{12}}{M_{11}}\right)v+h_2-\frac{M_{21}}{M_{11}}h_1.
\]

教材风格的虚拟控制输入：

\[
\boxed{v=-k_1q_2-k_2\dot q_2+k_3\dot q_1(E-E_d)}.
\]

- \(-k_1q_2\)：限制第二关节偏离伸展构型；
- \(-k_2\dot q_2\)：在虚拟加速度动力学中引入第二关节阻尼；
- \(k_3\dot q_1\tilde E\)：根据第一关节的摆动方向及能量误差，借惯性耦合调整其摆动。

由欠驱动方程 \(M_{11}\ddot q_1+M_{12}v+h_1=0\)，能量反馈对第一关节加速度的贡献为

\[
\ddot q_{1,E}=-\frac{M_{12}}{M_{11}}k_3\dot q_1\tilde E.
\]

在 \(M_{12}>0\) 的构型，能量不足时这一贡献近似负阻尼（增强摆动），能量过剩时近似正阻尼（抑制摆动）。**这只是该反馈项的局部动力学解释，不是总能量单调收敛的证明**。真实功率仍由 \(\dot E=\tau\dot q_2\) 决定。完整 PFL+PD 控制律没有直接的全局收敛保证。

## 5. Part 1–2 实验与结论

实验代码：`month2_day16_acrobot_swingup.py`（独立运行）；测试两种控制器，并核对能量功率恒等式。

用户一次实际运行的结果：

| 指标 | PFL | Direct |
|---|---:|---:|
| 末态总能量 | 10.0158 J | 19.6200 J |
| 末态能量误差 | -9.6042 J | 约 0 J |
| 最大力矩 | 38.708 Nm | 50 Nm |
| 力矩饱和样本 | 0 | 20 |
| 末态 \(q_1\) 包角误差 | +114.56° | -45.97° |
| 末态 \(q_2\) 包角误差 | +145.46° | -166.12° |

**解释**：Direct 的能量跟踪优于该组 PFL 参数，但两者均未达到直立静止状态。Direct 能量到位后可能继续沿保守系统轨道运动；PFL 的 PD 项约束相对关节运动，但不能保证整个系统耗散或捕获直立。

## 6. Part 3：Swing-up → LQR 分段控制

定义控制模式 \(\sigma\in\{\mathrm{swing},\mathrm{LQR}\}\)：

\[
\tau=\begin{cases}
\tau_{\rm swing}(x),&\sigma=\mathrm{swing},\\
-K\delta x,&\sigma=\mathrm{LQR}.
\end{cases}
\]

切换不仅看能量，还看**两关节角度误差与角速度**。实验采用滞回阈值：

| 条件 | 进入 LQR | 退出 LQR |
|---|---:|---:|
| 两关节包角误差 | 各小于 22° | 任一达到 40° |
| 两关节角速度绝对值 | 各小于 3 rad/s | 任一达到 5 rad/s |

实现为固定控制周期 \(\Delta t=0.02\,\mathrm s\)、零阶保持力矩、RK4 积分；力矩限制 \(|\tau|\le50\,\mathrm{Nm}\)。这组阈值是经验捕获区域，不是经证明的 LQR 吸引域。

参数调优：原 PFL 虚拟控制增益 \((k_1,k_2,k_3)=(10,5,1)\)；经验搜索后采用 \((6,3,1)\)。实际程序使用 \(\operatorname{wrap}(q_2)\) 替代未包角的 \(q_2\)，因此与最初的连续角度反馈存在差别：

\[
v=-6\operatorname{wrap}(q_2)-3\dot q_2+\dot q_1(E-E_d).
\]

在当前模型、初始 \(q_1=2.865^\circ\)、18 s 仿真条件下：

- **PFL + LQR**：成功捕获直立；首次进入 LQR 约 2.88 s，后续曾退出，约 6.18 s 再次切入后保持稳定。
- **Direct + LQR**：能量接近目标值，但未进入所设 LQR 捕获区域，因此没有完成直立稳定。
- 成功是对特定参数和初始条件的数值验证，不代表全局稳定或鲁棒性保证。

Part 3 完整代码：`month2_day16_part3_hybrid.py`；输出图：`day16_part3_final.png`。可用 `--tune` 搜索参数，或 `--pfl-gains 6 3 1` 指定增益。

## 7. 本日理论结论

1. **能量正确不等于状态正确**：\(E=E_d\) 定义能量水平集，不唯一确定四维状态。
2. **PFL 阻尼不等于整体机械能耗散**：真实能量变化取决于电机功率 \(\tau\dot q_2\)。
3. **PFL 和 Direct 是不同目标的起摆策略**：前者同时调节相对关节运动，后者直接调节总能量。
4. **LQR 负责最终局部稳定**；切换时必须考虑角度和速度，不能只看能量。
5. **进入人为阈值区域不等于进入真实吸引域**；前两次切换失败可能与惯性、力矩饱和、非线性耦合及捕获区域选择有关，不能只凭曲线归因。
6. **混合控制**不只来自地面碰撞，也可以来自控制器离散模式的切换；滞回可缓解边界抖振。

## 8. 与 Day 17 的衔接

Day 15–16 完成了简化欠驱动机器人的主动控制。Day 17 返回多自由度腿足机器人，以 Tedrake Chapter 5 为主线，从“可移动的推力器/足端”思想实验切入，将浮动基座动力学、接触力、质心线动量与质心角动量统一起来。随后 Day 18–19 构造逆动力学 QP 和 WBC。

---

**文件位置**：`ROBO_Learn_notes/Month2_Day16_Acrobot_Energy_Shaping_Hybrid_Control_Summary.md`
