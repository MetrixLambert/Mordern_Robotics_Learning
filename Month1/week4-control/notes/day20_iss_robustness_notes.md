# Day 20｜输入到状态稳定性（ISS）与鲁棒稳定性

## 1. 学习目标

- 区分无扰动渐近稳定性、输入到状态稳定性（ISS）与稳态跟踪精度。
- 利用 Lyapunov 方法理解有界外部扰动下的状态界。
- 识别状态相关模型误差如何改变闭环特征根，以及为什么名义 ISS 不自动保证鲁棒稳定性。
- 用 NumPy 仿真验证稳定边界、恒定扰动稳态误差及周期输入响应。

## 2. 从 CTC 到等效扰动系统

理想计算力矩控制（CTC）产生的误差动力学为

$$
\ddot e+K_d\dot e+K_pe=0.
$$

引入模型失配或外部扰动后，可写为

$$
\boxed{\ddot e+K_d\dot e+K_pe=w(t).}
$$

令 $x=[e^\top,\dot e^\top]^\top$，则

$$
\dot x=Ax+Bw,\qquad
A=\begin{bmatrix}0&I\\-K_p&-K_d\end{bmatrix},\quad
B=\begin{bmatrix}0\\I\end{bmatrix}.
$$

若 $A$ 为 Hurwitz，则无扰动系统指数稳定。对于实际机器人，等效扰动可能依赖于状态，不能不加条件地视为独立外部输入。

## 3. ISS 的含义

系统 $\dot x=f(x,w)$ 满足 ISS，是指存在 $\mathcal{KL}$ 函数 $\beta$ 和 $\mathcal K$ 函数 $\gamma$，使

$$
\boxed{\|x(t)\|\le\beta(\|x(0)\|,t)+\gamma(\|w\|_{\infty,[0,t]})}.
$$

第一项刻画初始误差的衰减，第二项刻画输入幅值对状态界的影响。

**ISS 不要求非零持续扰动下误差趋于零。** 例如

$$
\ddot e+k_d\dot e+k_pe=w_0,\qquad k_p,k_d>0,
$$

对于恒定 $w_0$，有

$$
\boxed{e_{\mathrm{ss}}=\frac{w_0}{k_p},\qquad \dot e_{\mathrm{ss}}=0.}
$$

因此，增大 $k_p$ 可以减小该例中的稳态位置误差，但实际控制中还需考虑力矩饱和、噪声和未建模动态。

## 4. ISS-Lyapunov 推导

对标量误差系统，令

$$
A=\begin{bmatrix}0&1\\-k_p&-k_d\end{bmatrix},\qquad B=\begin{bmatrix}0\\1\end{bmatrix}.
$$

当 $k_p,k_d>0$ 时，$A$ Hurwitz。取 $Q=Q^\top\succ0$，求解

$$
A^\top P+PA=-Q,\qquad P=P^\top\succ0.
$$

对 $V=x^\top Px$，有

$$
\begin{aligned}
\dot V
&=-x^\top Qx+2x^\top PBw\\
&\le-\lambda_{\min}(Q)\|x\|^2+2\|PB\|\|x\|\,|w|\\
&\le-\frac{\lambda_{\min}(Q)}2\|x\|^2
+\frac{2\|PB\|^2}{\lambda_{\min}(Q)}|w|^2.
\end{aligned}
$$

这给出标准的 ISS-Lyapunov 不等式。

### 为什么简单机械能函数不够直接？

取

$$
V_0=\frac12\dot e^2+\frac12k_pe^2,
\qquad
\dot V_0=-k_d\dot e^2+\dot e\,w.
$$

Young 不等式给出

$$
\dot V_0\le-\frac{k_d}{2}\dot e^2+\frac{w^2}{2k_d}.
$$

该不等式缺少关于 $e^2$ 的负项，**不能直接**得到完整状态的标准 ISS-Lyapunov 衰减不等式。这不意味着系统不满足 ISS，而是所选函数的导数估计不够强。可改用 $V=x^\top Px$，或在机械能函数中添加足够小的交叉项 $\varepsilon e\dot e$。

## 5. 状态相关模型误差与鲁棒稳定性

考虑

$$
\ddot e+k_d\dot e+k_pe=\delta_pe+\delta_d\dot e+d(t).
$$

整理为

$$
\boxed{\ddot e+(k_d-\delta_d)\dot e+(k_p-\delta_p)e=d(t).}
$$

定义

$$
k_{p,\mathrm{eff}}=k_p-\delta_p,\qquad
k_{d,\mathrm{eff}}=k_d-\delta_d.
$$

特征多项式为

$$
s^2+k_{d,\mathrm{eff}}s+k_{p,\mathrm{eff}}.
$$

对该标量二阶线性时不变系统，相对于外部输入 $d(t)$ 满足 ISS 的充要条件是

$$
\boxed{\delta_p<k_p,\qquad \delta_d<k_d.}
$$

- $\delta_p\uparrow k_p$：等效刚度趋于零，位置恢复能力变弱。
- $\delta_d\uparrow k_d$：等效阻尼趋于零，衰减速度变慢。
- 越过任一稳定边界：对应齐次系统不再 Hurwitz，不能满足 ISS。

对于满足上述严格条件的恒定扰动 $d(t)=d_0$，

$$
\boxed{e_{\mathrm{ss}}=\frac{d_0}{k_p-\delta_p}.}
$$

**关键区别：** 名义系统对独立输入 $w$ 满足 ISS，并不保证将状态相关项 $w=\delta_pe+\delta_d\dot e+d$ 接回闭环后仍然稳定。

## 6. NumPy 实验：实现与结果

完成的四项 TODO：

1. `closed_loop_coefficients()`：返回 $(k_p-\delta_p,k_d-\delta_d)$。
2. `is_iss()`：检查两个等效系数是否严格大于零。
3. `constant_disturbance_equilibrium()`：稳定时返回 $d_0/(k_p-\delta_p)$；不满足 ISS 时抛出 `ValueError`。
4. `compare_robustness_cases()`：执行六组对比仿真。

积分方式：半隐式 Euler；默认 $\Delta t=0.001\,\mathrm{s}$，时长 $12\,\mathrm{s}$，初始误差 $e(0)=0.3$，$\dot e(0)=0$；基准增益 $k_p=25,k_d=10$。

| 算例 | $(\delta_p,\delta_d)$ | 扰动（实际运行） | $e(12)$ (rad) | $\max|e|$ (rad) |
|---|---|---|---:|---:|
| `nominal` | $(0,0)$ | 恒定 $0.2$ | $+0.008000$ | $0.300000$ |
| `near_p_limit` | $(24,0)$ | 恒定 $0.2$ | $+0.230058$ | $0.300000$ |
| `p_unstable` | $(30,0)$ | 恒定 $0.2$ | $+99.947037$ | $99.947037$ |
| `near_d_limit` | $(0,9)$ | 默认恒定 $0.2$ | $+0.007277$ | $0.300000$ |
| `d_unstable` | $(0,12)$ | 默认恒定 $0.2$ | $-35830.487678$ | $35830.487678$ |
| `periodic` | $(5,2)$ | $0.2\sin(3t)$ | $-0.002257$ | $0.300000$ |

**测试结果：** 5/5 组测试通过（包括有界正弦输入下响应有限的检查）。

### 结果解读与数值注意事项

- `nominal`：理论稳态 $e_{\mathrm{ss}}=0.2/25=0.008$ rad，与最终值一致。
- `near_p_limit`：等效刚度仅为 $1$，理论稳态应为 $0.2$ rad。12 秒最终值 $0.230058$ rad，说明此时仍存在尚未完全衰减的瞬态，**不能将 $e(12)$ 当作稳态值**。对应特征根为 $-5\pm\sqrt{24}$，慢极点约 $-0.101$，收敛明显变慢。
- `p_unstable`：等效刚度 $-5$，存在正实特征根，误差发散。
- `near_d_limit`：等效阻尼为 $1>0$，仍稳定但衰减较慢；其输入为默认恒定 $0.2$，不是零扰动。
- `d_unstable`：等效阻尼为 $-2$，产生负阻尼失稳。
- `periodic`：闭环满足 ISS，周期输入导致有界响应；末端时刻误差值不是周期响应幅值。

**代码与题目说明的细微差异：** `p_unstable`、`near_d_limit` 和 `d_unstable` 的题目注释要求零扰动，但提交的调用未指定 `kind="zero"`，因此使用了 `simulate()` 默认的恒定扰动。若希望实验严格对应注释，可将三处调用改为 `kind="zero"`。这不改变稳定边界的理论结论，但会影响具体仿真数值。

## 7. 与机器人学习的联系

在机器人 CTC 中，真实等效扰动可写为

$$
w=M(q)^{-1}\left[(\hat M-M)v+(\hat c-c)+(\hat g-g)+d\right],
\qquad v=\ddot q_d-K_d\dot e-K_pe.
$$

该 $w$ 通常依赖于状态、期望轨迹和控制指令。若引入神经网络残差补偿，不能仅因为名义误差系统满足 ISS，就断言任意学习补偿均安全。应进一步约束或分析学习残差的幅值、增长率及其与反馈闭环的耦合，例如通过鲁棒 Lyapunov 分析或小增益条件。

## 8. 核心结论

1. **ISS 不等于零稳态误差。** 有界持续输入可以产生非零但受控的误差。
2. **ISS 是关于完整状态的输入—状态界。** 单独看到速度耗散项不足以直接证明完整状态 ISS。
3. **状态相关模型误差会改变闭环动力学。** 对本次标量模型，$k_p>\delta_p$、$k_d>\delta_d$ 是严格稳定边界。
4. **靠近稳定边界不等于性能良好。** 系统虽然仍稳定，但收敛可能极慢，稳态误差也可能显著放大。
5. **仿真末端值不一定是稳态值。** 需要结合特征根、理论平衡点和仿真时长判断。

**Day 20：全部测试通过。下一阶段：将 ISS 和鲁棒性分析推广到更贴近机器人动力学的控制与学习补偿问题。**
