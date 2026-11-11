# ROBO_Learn · Day 9 学习总结：浮动基座 WBC-QP、点接触与平足接触

> 学习主题：在二维浮动基座单腿模型中，用二次规划（QP）同时处理动力学、固定接触、摩擦、足底力矩与加速度任务。  
> 本文记录本次实际推导、程序输出和结论；不把单步 QP 的结果等同于闭环站立稳定性。

## 1. 本日目标与完成情况

- 建立二维浮动基座动力学方程，区分关节驱动力矩与地面接触 wrench。
- 分别实现 **Point contact** 与 **Flat foot** 两种接触模型。
- 通过 `dynamics residual`、`contact residual` 验证等式约束。
- 分析任务加速度之间的动力学耦合及 QP 权重折中。
- 推导 CoP 公式和脚底支撑约束，并用不同初始站姿验证静态平衡条件。
- 扫描水平加速度命令，观察 CoP 饱和与最优解的分段变化。

**结论：Day 9 核心实验完成。** CoP 安全裕度与更复杂的权重调参不再继续，后续转向 MuJoCo 中的闭环 WBC。

## 2. 二维浮动基座模型

广义坐标、速度与加速度：

\[
q=[x,z,\theta,q_1,q_2]^\top,\qquad v=\dot q,\qquad a=\ddot q.
\]

- \(x,z\)：躯干质心位置；\(\theta\)：躯干角度。
- \(q_1,q_2\)：髋、膝相对关节角。
- 两段腿长 \(l_1=l_2=0.5\,\mathrm m\)；腿部质量忽略。
- 躯干质量 \(m=10\,\mathrm{kg}\)，转动惯量 \(I=0.6\,\mathrm{kg\,m^2}\)，\(g=9.81\,\mathrm{m/s^2}\)。

定义绝对杆角：

\[
\alpha_1=\theta+q_1,\qquad \alpha_2=\theta+q_1+q_2.
\]

脚底中心位置：

\[
p_f=\begin{bmatrix}
x+l_1\sin\alpha_1+l_2\sin\alpha_2\\
z-l_1\cos\alpha_1-l_2\cos\alpha_2
\end{bmatrix}.
\]

为简化表达，记：

\[
\begin{aligned}
A&=l_1\cos\alpha_1+l_2\cos\alpha_2, &B&=l_1\sin\alpha_1+l_2\sin\alpha_2,\\
C&=l_2\cos\alpha_2, &D&=l_2\sin\alpha_2.
\end{aligned}
\]

平动接触 Jacobian：

\[
J=\frac{\partial p_f}{\partial q}=
\begin{bmatrix}
1&0&A&A&C\\
0&1&B&B&D
\end{bmatrix}.
\]

其速度相关偏置项：

\[
\dot Jv=\begin{bmatrix}
-l_1\sin\alpha_1\,\dot\alpha_1^2-l_2\sin\alpha_2\,\dot\alpha_2^2\\
l_1\cos\alpha_1\,\dot\alpha_1^2+l_2\cos\alpha_2\,\dot\alpha_2^2
\end{bmatrix}.
\]

采用本次简化的质量矩阵和重力项：

\[
M=\operatorname{diag}(m,m,I,0,0),\quad
h=[0,mg,0,0,0]^\top,\quad
S^\top=\begin{bmatrix}0_{3\times2}\\I_2\end{bmatrix}.
\]

**模型局限：** 腿部无质量使 \(M\) 奇异；这是一种便于理解浮动基座与接触力分配的教学模型，而非完整的人形机器人动力学模型。

## 3. 从虚功到接触广义力

同一个地面接触作用既可以在脚底接触坐标中表示，也可以在广义坐标中表示。对于任意虚位移：

\[
\delta p_w=J_w\delta q.
\]

接触 wrench 的虚功：

\[
\delta W_c=w_c^\top\delta p_w
=w_c^\top J_w\delta q
=(J_w^\top w_c)^\top\delta q.
\]

因此，接触力在广义坐标中的等效表达为：

\[
\boxed{Q_c=J_w^\top w_c}.
\]

这不是额外产生的一组物理力，而是**同一个物理作用的广义力表示**，不能与原接触力再重复相加。\(p_w\) 是 \(q\) 的函数，并不要求二者可逆或自由度相等。固定接触时，实际允许的虚位移满足 \(J_w\delta q=0\)，理想接触反力对这些允许虚位移做零虚功；但广义力映射 \(J_w^\top w_c\) 仍成立。

统一动力学为：

\[
\boxed{Ma+h=S^\top\tau+J_w^\top w_c}.
\]

这里 \(\tau\) 是电机施加的关节力矩，\(w_c\) 是地面施加在脚上的接触 wrench，二者是不同的物理作用。

## 4. Point contact 与 Flat foot

### 4.1 Point contact

只约束脚底中心的平动：

\[
Ja+\dot Jv=0.
\]

接触变量：

\[
w_c=[f_x,f_z]^\top.
\]

摩擦与单边接触约束：

\[
f_z\ge0,\qquad |f_x|\le\mu f_z.
\]

**特点：** 接触点无法提供独立的接触力矩，躯干转动必须通过水平/竖直接触力的力臂和系统运动来平衡。

### 4.2 Flat foot

本次使用脚底角度：

\[
\phi_f=\theta+q_1+q_2=\alpha_2,
\]

它是相对水平脚底朝向的简化角度定义，已吸收几何固定偏置。对应：

\[
J_\phi=[0,0,1,1,1],\qquad \dot J_\phi v=0.
\]

脚底完整平面位姿的 Jacobian 与固定接触约束：

\[
J_w=\begin{bmatrix}J\\J_\phi\end{bmatrix},\qquad
\boxed{J_wa+\begin{bmatrix}\dot Jv\\0\end{bmatrix}=0}.
\]

接触 wrench：

\[
w_c=[f_x,f_z,m_c]^\top.
\]

其中 \(m_c\) 是**地面对脚底中心的接触力矩**，不是髋/膝电机力矩。广义接触力：

\[
J_w^\top w_c=
\begin{bmatrix}
f_x\\f_z\\Af_x+Bf_z+m_c\\Af_x+Bf_z+m_c\\Cf_x+Df_z+m_c
\end{bmatrix}.
\]

增加足底半长 \(d=0.12\,\mathrm m\) 的 CoP 约束：

\[
\boxed{|m_c|\le df_z}.
\]

## 5. CoP 公式与物理来源

脚底沿长度方向坐标 \(s\in[-d,d]\)，法向压力分布为 \(p(s)\ge0\)。总竖直压力和脚底中心力矩为：

\[
f_z=\int_{-d}^d p(s)\,ds,\qquad
m_c=-\int_{-d}^d s\,p(s)\,ds.
\]

这里力矩符号沿用**本次代码的约定**。令分布式压力等效为作用在 CoP 的集中力：

\[
\boxed{x_{\mathrm{CoP}}=-\frac{m_c}{f_z}},\qquad f_z>0.
\]

于是：

\[
x_{\mathrm{CoP}}=
\frac{\int_{-d}^d s\,p(s)\,ds}{\int_{-d}^d p(s)\,ds}.
\]

因为它是非负压力的加权平均，必然位于脚底内：

\[
|x_{\mathrm{CoP}}|\le d
\quad\Longleftrightarrow\quad
\boxed{-df_z\le m_c\le df_z}.
\]

当 \(f_z=0\) 时 CoP 比值未定义；QP 中的力矩约束仍有意义。

## 6. 为什么不同加速度任务会耦合？

由浮动基座的前三行动力学：

\[
ma_x=f_x,\qquad m(a_z+g)=f_z,
\]

\[
Ia_\theta=Af_x+Bf_z+m_c.
\]

消去接触力：

\[
\boxed{Ia_\theta=mAa_x+mB(a_z+g)+m_c}.
\]

再用 \(m_c=-x_{\mathrm{CoP}}f_z\)：

\[
\boxed{
Ia_\theta=mAa_x+m(B-x_{\mathrm{CoP}})(a_z+g)
}.
\]

等价地：

\[
\boxed{
x_{\mathrm{CoP}}=B+\frac{Aa_x-(I/m)a_\theta}{a_z+g}
},\quad a_z+g>0.
\]

所以，\(a_x,a_z,a_\theta\) 并不能在固定接触下任意独立指定；它们还必须满足 CoP、摩擦等可行性约束。QP 的权重只能改变折中点，不能突破物理可行域。

## 7. 静态站姿实验

### 7.1 原站姿：静态平衡不可行

取 \(\theta=0,q_1=0.4,q_2=-0.4\)，则 \(\alpha_2=0\)，

\[
B=0.5\sin0.4\approx0.19471\,\mathrm m.
\]

静止时 \(a=0\)，要求 \(f_z=mg\) 且 \(m_c=-mgB\)，即

\[
m_c\approx-19.10\,\mathrm{N\,m},
\]

但足底最大允许力矩为

\[
dmg=0.12\times98.1=11.772\,\mathrm{N\,m}.
\]

故 **原站姿不能在保持脚底完整固定接触的情况下静态平衡**。即使目标水平加速度为零，QP 也会给出非零加速度以折中各项任务。

### 7.2 新站姿：Flat foot 可以静态平衡

改为 \(\theta=0,q_1=0.1,q_2=-0.1\)：

\[
B=0.5\sin0.1\approx0.0499167\,\mathrm m<d.
\]

本次程序输出（\(a_{x,\mathrm{cmd}}=0\)）：

```text
POINT CONTACT
  a = [-0.47262 -0.02375  0.28429  0.19148 -0.00392]
  objective = 23.297683

FLAT FOOT
  a = [-1.e-05 -0.e+00  0.e+00  2.e-05 -2.e-05]
  wrench = [-1.200000e-04  9.809999e+01 -4.896710e+00]
  tau = [-0.       4.89677]
  objective = 0.096476
  CoP = 0.04991545178732397 m
```

理论预测：

\[
x_{\mathrm{CoP}}=B\approx0.0499167\,\mathrm m,
\qquad m_c=-mgB\approx-4.8968\,\mathrm{N\,m},
\]

与数值结果一致。Point contact 仍无法完全静止，是因为它缺少抵消重力力矩所需的独立接触力矩。

### 7.3 静态平衡与支撑区域

在腿部无质量、无其他外力、关节力矩足够且脚底保持水平接触的条件下，整机 COM 就是躯干质心。静态平衡要求 COM 的水平投影落入脚底支撑区间：

\[
\boxed{|x_{\mathrm{COM}}-x_{\mathrm{foot}}|=|B|\le d}.
\]

这一点也解释了大型背包高达模型的站姿问题：背包会移动整机 COM，双脚前后错开可扩大前后支撑范围，但若左右间距过小，支撑多边形会狭长，容易向两脚连线的近似垂直方向倾倒。真实模型还受关节保持力矩、摩擦、零件弹性影响。

## 8. 水平加速度跟踪实验

新站姿下，保持 QP 原权重，对 \(a_{x,\mathrm{cmd}}\) 扫描。

| \(a_{x,\mathrm{cmd}}\) (m/s²) | Point \(a_x\) | Flat \(a_x\) | Flat \(a_\theta\) (rad/s²) | Flat CoP (m) |
|---:|---:|---:|---:|---:|
| 0 | -0.47262 | ≈0 | ≈0 | 0.049915 |
| 1 | -0.43538 | 0.70471 | 0.17457 | 0.120000 |
| 3 | -0.36089 | 0.77548 | 1.34282 | 0.120000 |
| 8 | -0.17467 | 0.95240 | 4.26343 | 0.120000 |

主要观察：

1. Flat foot 明显改善了静态平衡和水平加速度跟踪。
2. 从命令 \(a_x=1\) 开始，CoP 已经到达 \(+d\) 边缘。
3. CoP 饱和以后，进一步增大 \(a_x\) 必须付出更大的躯干角加速度代价。
4. 求解器状态为 `optimal`、动力学与接触残差接近机器精度，只能说明求得了满足当前数学模型约束的最优解，**不意味着运动目标被精确实现**。

### 8.1 预测 CoP 饱和临界值

在 \(a_z=a_\theta=0\) 时：

\[
x_{\mathrm{CoP}}=B+\frac{A}{g}a_x.
\]

当前 \(A=0.5\cos0.1+0.5\approx0.99750\)，\(B\approx0.04992\)。令 \(x_{\mathrm{CoP}}=d=0.12\)：

\[
\boxed{
a_{x,\max}=\frac{g(d-B)}{A}\approx0.689\,\mathrm{m/s^2}
}.
\]

**注意：** 这是在保持 \(a_z=a_\theta=0\) 条件下的临界加速度，不是允许躯干旋转时的绝对最大水平加速度。实际 QP 的切换点还受权重和正则项影响。

### 8.2 密集扫描验证

使用 `np.linspace(0, 1.2, 25)` 扫描命令加速度，绘制三张曲线：

- \(a_x\) 对命令值：前段近似线性跟踪，后段增长明显变慢。
- \(a_\theta\) 对命令值：前段约为零，CoP 饱和后明显增大。
- CoP 对命令值：从约 \(0.05\,\mathrm m\) 上升，随后停留在 \(0.12\,\mathrm m\)。

观察到分界点约在 \(0.6\sim0.7\,\mathrm{m/s^2}\)，与理论预测 \(0.689\,\mathrm{m/s^2}\) 相符。这是 QP **活跃约束集发生切换**的直观例子。

## 9. WBC-QP 的标准结构（本次模型）

决策变量：

\[
a\in\mathbb R^5,\quad\tau\in\mathbb R^2,\quad
w_c\in\mathbb R^2\text{（点接触）或}\mathbb R^3\text{（平足）}.
\]

代表性目标函数：

\[
\min_{a,\tau,w_c}
\ w_x(a_x-a_{x,\mathrm{cmd}})^2+w_z(a_z-a_{z,\mathrm{cmd}})^2
+w_\theta(a_\theta-a_{\theta,\mathrm{cmd}})^2
+\epsilon_w\|w_c\|^2+\epsilon_\tau\|\tau\|^2.
\]

约束：

\[
\begin{cases}
Ma+h=S^\top\tau+J_w^\top w_c,\\
J_wa+\dot J_wv=0,\\
f_z\ge0,\quad-\mu f_z\le f_x\le\mu f_z,\\
-df_z\le m_c\le df_z\quad\text{（Flat foot）}.
\end{cases}
\]

这是一个凸二次规划（在目标权重非负且约束保持线性的条件下）。目标函数体现**希望实现的运动**，等式与不等式约束体现**动力学和接触物理允许的运动**。

## 10. 本日最重要的认识

1. **虚功是广义力映射的基础：** \(Q_c=J_w^\top w_c\) 描述的是同一接触作用在不同变量下的表达。
2. **浮动基座欠驱动：** 基座坐标没有直接驱动器，基座运动必须依赖接触反力与关节驱动共同满足动力学。
3. **接触模型决定能力：** Point contact 没有独立接触力矩；Flat foot 允许接触力矩，但受 CoP 约束。
4. **站姿决定静态可行性：** 在本模型中，\(|B|\le d\) 是单脚静态平衡的核心几何条件。
5. **QP 不是万能的跟踪器：** 当任务与物理约束冲突时，它求解的是最优折中，而不是违反约束来精确跟踪。
6. **CoP 饱和产生分段行为：** 当 CoP 到达边缘，继续提高水平加速度通常需要更大的躯干角加速度。
7. **瞬时最优不等于闭环稳定：** 单步求解验证了数学关系，但还没有证明机器人在连续仿真中能长期站立或抗扰动。

## 11. 下一步：Day 10

**主题：MuJoCo 中的闭环 WBC 控制。**

从“给定一个 \(q,v\)，求一次 \(a,\tau,w_c\)”推进到：

1. 在 MuJoCo 中建立具有质量、惯量和脚底碰撞几何的机器人。
2. 每个控制周期读取 \(q,v\) 和接触状态。
3. 计算模型量 \(M,h,J,\dot Jv\)，构造并求解 QP。
4. 将关节力矩施加到仿真，循环更新状态。
5. 观察站立保持、姿态偏差、接触变化和数值稳定性。

需要特别注意：当前 Day 9 的无质量腿教学模型与 MuJoCo 中的真实正定惯量模型不同；进入 Day 10 时应检查动力学一致性、接触模式和状态反馈，而不是直接假设单步 QP 必然产生稳定闭环。
