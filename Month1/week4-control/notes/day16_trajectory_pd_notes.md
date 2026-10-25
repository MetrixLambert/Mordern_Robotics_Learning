# Day 16｜关节空间轨迹规划与 PD 跟踪控制

## 1. 学习目标与验收

- 使用五次多项式生成静止到静止的关节轨迹及其一、二阶导数。
- 在竖直平面二连杆模型上实现纯 PD、逆动力学前馈 + PD、实时重力补偿 + PD。
- 区分「期望轨迹是闭环解」与「误差渐近收敛」。
- 对比跟踪误差与关节力矩，理解模型补偿的作用及局限。

**验收：5 组自动测试全部通过。**

## 2. 动力学模型与误差定义

沿用 Week 3 的二连杆模型：

$$
M(q)\ddot q+C(q,\dot q)\dot q+g(q)=\tau,\qquad M(q)\succ0.
$$

轨迹规划输出 \(q_d(t),\dot q_d(t),\ddot q_d(t)\)。在本笔记统一定义

$$
e=q-q_d,\qquad \dot e=\dot q-\dot q_d.
$$

因此 PD 反馈写为 \(-K_pe-K_d\dot e\)，与程序中的 \(K_p(q_d-q)+K_d(\dot q_d-\dot q)\) 完全一致。

## 3. 五次多项式轨迹

设运动时间 \(T>0\)，\(s=t/T\)，满足

$$
q_d(0)=q_0,\quad q_d(T)=q_f,\quad
\dot q_d(0)=\dot q_d(T)=0,\quad
\ddot q_d(0)=\ddot q_d(T)=0.
$$

取五次多项式可得

$$
\begin{aligned}
q_d(t)&=q_0+(q_f-q_0)(10s^3-15s^4+6s^5),\\
\dot q_d(t)&=\frac{q_f-q_0}{T}(30s^2-60s^3+30s^4),\\
\ddot q_d(t)&=\frac{q_f-q_0}{T^2}(60s-180s^2+120s^3).
\end{aligned}
$$

向量公式逐关节成立。区间外保持端点位置，速度和加速度为零。五次多项式使加速度在端点连续地归零，但不要求 jerk 在端点为零。

## 4. 三种控制器

### 4.1 纯 PD

$$
\boxed{\tau_{\rm PD}=-K_pe-K_d\dot e.}
$$

不补偿重力和运动中的惯性、速度耦合。对常值目标 \(q^*\)，静态平衡满足

$$
K_p(q_{\rm eq}-q^*)+g(q_{\rm eq})=0,
$$

故通常有静态位置偏差。重力为保守力，并非静止时持续做功的外部输入。

### 4.2 逆动力学前馈 + PD（FF+PD）

$$
\boxed{\tau_{\rm FF+PD}=M(q_d)\ddot q_d+C(q_d,\dot q_d)\dot q_d+g(q_d)-K_pe-K_d\dot e.}
$$

前馈在**期望状态**处计算。模型准确且初始状态严格匹配时，\(q=q_d\) 是闭环系统的一条解轨迹；这不等于证明任意初始误差都收敛。

其精确误差动力学为

$$
M(q)\ddot e+C(q,\dot q)\dot e+K_d\dot e+K_pe=\Delta,
$$

$$
\Delta=[M(q_d)-M(q)]\ddot q_d+[C(q_d,\dot q_d)-C(q,\dot q)]\dot q_d+g(q_d)-g(q).
$$

所以前馈没有完全消除偏离轨迹时的非线性误差项。

### 4.3 实时重力补偿 + PD（Gravity+PD）

$$
\boxed{\tau_{\rm G+PD}=g(q)-K_pe-K_d\dot e.}
$$

重力在**实际位置**处补偿。对常值目标，可通过机器人动力学的斜对称结构进行 Lyapunov 分析；对时变目标，缺少惯性和速度耦合前馈，运动中通常存在跟踪误差。轨迹结束后仍有误差，不代表最终稳态误差一定非零。

## 5. 数值实验结果

仿真使用半隐式 Euler，三种控制器采用相同轨迹、初始状态和增益，轨迹结束后继续保持目标 1 秒。

| 控制器 | 跟踪 RMSE (rad) | 最终误差范数 (rad) | 最大单关节力矩绝对值 (N·m) |
|---|---:|---:|---:|
| PD | 0.24639 | 0.34992 | 80.53 |
| FF+PD | 0.00012 | 0.00004 | 63.98 |
| Gravity+PD | 0.07767 | 0.05689 | 65.83 |

观察：纯 PD 重力偏差明显；Gravity+PD 改善静态补偿但存在动态跟踪误差；FF+PD 在准确模型、匹配初始状态下几乎零误差。残余误差包含离散积分误差。上述力矩最大值不是能耗指标，且结果只针对本组实验参数。

## 6. Lyapunov 分析中的关键区分

**零误差轨迹可行性**：将 \(q=q_d,\dot q=\dot q_d\) 代入方程，只能检验零误差是否构成闭环解。

**误差稳定性**：必须推导非零 \(e,\dot e\) 的完整动态，构造误差 Lyapunov 函数并分析其导数；不能通过先假设 \(\dot q=0\) 来证明一般的时变轨迹跟踪稳定性。

对于 FF+PD，可尝试

$$
V=\tfrac12\dot e^TM(q)\dot e+\tfrac12e^TK_pe,
$$

利用 \(\dot M-2C\) 的斜对称性质得

$$
\dot V=-\dot e^TK_d\dot e+\dot e^T\Delta.
$$

\(\Delta\) 在零误差处消失，但一般不为零，因此不能直接据此声称全局渐近稳定。

## 7. 与后续内容的联系

采用实际状态进行模型补偿的 Computed Torque Control（CTC）：

$$
\tau=M(q)(\ddot q_d-K_d\dot e-K_pe)+C(q,\dot q)\dot q+g(q).
$$

在模型准确时得到

$$
\ddot e+K_d\dot e+K_pe=0,
$$

这是反馈线性化，而非单纯的期望轨迹前馈。若使用估计模型 \(\hat M,\hat C,\hat g\)，则误差动态会出现模型残差；这是后续鲁棒、自适应及学习补偿的切入点。

## 8. 代码实现小结

完成五个 TODO：`quintic_trajectory`、`control_pd`、`control_ff_pd`、`control_gravity_pd`、`compare_controllers`。五组测试均 PASS。建议在轨迹函数开头显式检查 `T > 0`。

**Day 16 完成。下一步：Day 17，机器人 PD 与重力补偿的稳定性结构。**
