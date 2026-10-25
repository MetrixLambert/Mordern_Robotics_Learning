# Day 26｜力控制、混合运动/力控制与阻抗控制

> ROBO_Learn · Month 1 · Day 26  
> 状态：**已完成**（理论问答与 MuJoCo 实验均已验证）  
> 核心主题：从自由空间位置控制过渡到机器人与外界交互时的柔顺控制。

## 1. 本日学习成果

- 理解力控制、混合运动/力控制、阻抗控制的目标与适用情形。
- 理解阻抗模型中的期望惯量、阻尼、刚度及其对稳态和瞬态响应的影响。
- 推导关节空间和任务空间的完整惯量整形阻抗控制律。
- 从储能函数理解无源性及其在机器人—环境交互中的作用。
- 在 MuJoCo 2R 机械臂上完成重力补偿 PD、完整关节阻抗、完整任务空间阻抗的对照实验。

## 2. 为什么单纯的位置控制不够

自由空间中的计算力矩控制（CTC）可通过逆动力学补偿得到近似线性的误差动力学：

$$
\tau=M(q)(\ddot q_d-K_D\dot e-K_Pe)+h(q,\dot q),\quad
\ddot e+K_D\dot e+K_Pe=0,
$$

其中 $e=q-q_d$、$h=C(q,\dot q)\dot q+g(q)$。当末端受到刚性环境约束时，期望位置可能无法到达，持续的位置误差会驱动控制器继续施力，因此接触任务需要显式考虑力和柔顺性。

### 2.1 力控制（Force Control）

目标是调节接触力 $F$，例如维持打磨或擦拭时的法向压力。力误差为 $e_F=F_d-F$。末端力与关节广义力矩的映射为：

$$\tau_{\rm ext}=J(q)^TF_{\rm ext}.$$

这一映射来自虚功关系 $F^T\delta x=\tau^T\delta q$。**力反馈控制只有在建立接触后才有明确的接触力调节对象**；实际控制器还需要考虑动力学补偿、接触方向与传感器信息。

### 2.2 混合运动/力控制（Hybrid Motion/Force Control）

典型场景：沿墙面切向跟踪轨迹，同时在法向维持接触压力。若法向为 $z$，一个简化的选择矩阵例子为：

$$S_m=\operatorname{diag}(1,1,0),\qquad S_f=\operatorname{diag}(0,0,1).$$

真实任务的运动/力子空间由接触约束几何决定，不能在任意坐标系中机械地套用上述矩阵。

### 2.3 阻抗控制（Impedance Control）

不直接要求“位置误差始终为零”或“接触力始终为某值”，而是规定**外力与运动响应之间的动态关系**：

$$\boxed{M_d\ddot e+D_d\dot e+K_de=F_{\rm ext}}$$

其中 $e=x-x_d$（关节空间时将 $x,F$ 换成 $q,\tau$）。$M_d$ 为期望惯量，$D_d$ 为期望阻尼，$K_d$ 为期望刚度。注意此处 $K_d$ 是**刚度**，不要与 CTC 中常用的微分增益 $K_D$ 混淆。

恒定外力下，若系统收敛到静态平衡：

$$\boxed{e_\infty=K_d^{-1}F_0}.$$

一维情形的自然频率和阻尼比：

$$\omega_n=\sqrt{K_d/M_d},\qquad \zeta=\frac{D_d}{2\sqrt{M_dK_d}}.$$

减小刚度会增大同一外力下的静态位移；改变阻尼主要改变瞬态响应；改变惯量影响加速度与自然频率。阻尼过大不一定收敛更快。

## 3. 关节空间阻抗控制

机器人受外部广义力矩时：

$$M(q)\ddot q+h(q,\dot q)=\tau+\tau_{\rm ext}.$$

希望闭环满足：

$$M_d\ddot e+D_d\dot e+K_de=\tau_{\rm ext},\qquad e=q-q_d.$$

理想完整阻抗控制律：

$$\boxed{\tau=M(q)\left[\ddot q_d+M_d^{-1}(\tau_{\rm ext}-D_d\dot e-K_de)\right]+h(q,\dot q)-\tau_{\rm ext}}.$$

该控制律在模型精确、外力已知、执行器可准确施加力矩等条件下实现指定的闭环惯量、阻尼与刚度。

### 与重力补偿 PD 的区别

固定目标下：

$$\tau=g(q)-K_Pe-K_D\dot e.$$

代入动力学得到：

$$M(q)\ddot e+[C(q,\dot q)+K_D]\dot e+K_Pe=\tau_{\rm ext}.$$

这种控制已经具有弹簧—阻尼式柔顺行为，但**真实惯量 $M(q)$ 并未被整形成指定的 $M_d$**。完整阻抗控制与 PD 可能有相同的稳态位移，但瞬态不同。

课堂算例：$\tau_{\rm ext}=[1,0]^T\,\mathrm{N\,m}$，$K_P=100I$，静态平衡误差为 $[0.01,0]^T\,\mathrm{rad}$。

## 4. 无源性（Passivity）

一维阻抗系统：

$$m\ddot e+d\dot e+ke=f_{\rm ext},\qquad m,d,k>0.$$

选择储能函数（动能与弹性势能）：

$$V=\tfrac12m\dot e^2+\tfrac12ke^2.$$

则：

$$\boxed{\dot V=f_{\rm ext}\dot e-d\dot e^2\leq f_{\rm ext}\dot e}.$$

将输入定义为 $u=f_{\rm ext}$、输出定义为 $y=\dot e$，则 $uy$ 是外界输入的瞬时功率。存在非负储能函数满足 $\dot V\leq u^Ty$，意味着系统不会在超过外界输入和初始储能的意义上凭空产生能量，这就是无源性的核心。

- **无源性不等于渐近稳定性**：无阻尼质量—弹簧系统可以无源，但零输入时持续振荡。
- 对 $f_{\rm ext}=0$，$\dot V=-d\dot e^2\leq0$；结合 LaSalle 不变性原理可证明原点渐近稳定。
- **负阻尼**系统 $m\ddot x-d\dot x=f$ 对同一动能函数满足 $\dot V=f\dot x+d\dot x^2$；零输入时速度可指数增长，反映内部能量注入。
- 机器人和环境若以功率守恒的方式互联且各自无源，则可将储能相加进行整体稳定性分析；渐近稳定还需额外条件。

## 5. 任务空间动力学与阻抗

末端位置运动学：

$$x=f(q),\quad \dot x=J(q)\dot q,\quad \ddot x=J(q)\ddot q+\dot J(q,\dot q)\dot q.$$

对于非冗余、非奇异的机械臂，定义任务空间惯量：

$$\boxed{\Lambda(q)=J^{-T}M(q)J^{-1}}.$$

对于满行秩的任务 Jacobian，更一般地：

$$\boxed{\Lambda(q)=(JM^{-1}J^T)^{-1}}.$$

非冗余情形下，任务空间动力学可写为：

$$\Lambda\ddot x+\eta=F_{\rm cmd}+F_{\rm ext},$$

其中：

$$\eta=J^{-T}\bigl(h-MJ^{-1}\dot J\dot q\bigr),\qquad \tau=J^TF_{\rm cmd}.$$

期望阻抗：

$$M_d\ddot e_x+D_d\dot e_x+K_de_x=F_{\rm ext},\qquad e_x=x-x_d.$$

因此：

$$\ddot x_{\rm cmd}=\ddot x_d+M_d^{-1}(F_{\rm ext}-D_d\dot e_x-K_de_x),$$

$$\boxed{F_{\rm cmd}=\Lambda\ddot x_{\rm cmd}+\eta-F_{\rm ext},\qquad \tau=J^TF_{\rm cmd}}.$$

对于本次 2R 平面机械臂，直接在关节加速度层实现等价的控制更方便：

$$\ddot q_{\rm cmd}=J^{-1}(\ddot x_{\rm cmd}-\dot J\dot q),$$

$$\tau=M\ddot q_{\rm cmd}+h-J^TF_{\rm ext}.$$

**重要限制**：需要避开 Jacobian 奇异构型；冗余机械臂还涉及零空间力矩设计；模型误差、力测量误差、时延、饱和和实际接触会影响理想阻抗的实现。

### 任务空间惯量的物理意义

课堂例子：

$$M=\begin{bmatrix}2&0\\0&1\end{bmatrix},\qquad J=\begin{bmatrix}1&0\\0&2\end{bmatrix}.$$

于是：

$$\Lambda=\begin{bmatrix}2&0\\0&1/4\end{bmatrix}.$$

在忽略偏置和其他输入的条件下，沿 $x$ 方向施加 $1\,\mathrm N$ 产生 $0.5\,\mathrm{m/s^2}$ 的加速度，沿 $y$ 方向施加同样大小的力产生 $4\,\mathrm{m/s^2}$。**任务空间惯量越小，相同方向外力引起的加速度越大**。它综合体现了质量分布与构型，不等同于速度可操作度。

## 6. MuJoCo 实验记录

本日已完成三个脚本，并已在本地修改参数、观察并验证结果：

| 文件 | 实验内容 | 主要结论 |
|---|---|---|
| `day26_joint_impedance.py` | 重力补偿 PD + 外部关节力矩 | 稳态偏移随刚度倒数变化；阻尼主要影响瞬态 |
| `day26_full_impedance.py` | 完整关节阻抗与 PD 对照 | 期望惯量 $M_d$ 可以独立设计；相同静态刚度不意味着相同动态响应 |
| `day26_task_space_impedance.py` | 2R 末端任务空间阻抗与 PD 对照 | 在笛卡尔方向直接指定刚度、阻尼与期望惯量 |

第三个脚本中，机械臂运动平面是世界坐标的 **$x$–$z$ 平面**，而非 $x$–$y$ 平面。默认设置：

$$K=\operatorname{diag}(200,50)\ \mathrm{N/m},$$

$$D=\operatorname{diag}(30,14)\ \mathrm{N\,s/m},\quad M_d=\operatorname{diag}(1.5,0.8)\ \mathrm{kg},$$

$$F_{\rm ext}=[0,5]^T\ \mathrm N.$$

理论上，$z$ 方向静态偏移为 $5/50=0.1\,\mathrm m$（在目标位置及可达工作空间等条件允许的情况下）。

运行示例：

```powershell
python day26_joint_impedance.py
python day26_full_impedance.py
python day26_task_space_impedance.py
python day26_task_space_impedance.py --mode pd
python day26_task_space_impedance.py --viewer --mode full
```

### 关键 MuJoCo API

```python
M = np.zeros((model.nv, model.nv), dtype=np.float64)
mujoco.mj_fullM(model, data, M)
h = data.qfrc_bias.copy()

jacp = np.zeros((3, model.nv))
jacr = np.zeros((3, model.nv))
mujoco.mj_jacSite(model, data, jacp, jacr, site_id)
J = jacp[[0, 2], :]  # x-z 平面

data.qfrc_applied[:] = J.T @ F_ext
```

- `qfrc`：与广义坐标/速度对应的 generalized force（广义力），旋转关节对应力矩。
- `qfrc_bias`：重力、科里奥利和离心等动力学偏置，**一般不能直接等同于重力项**。
- `qfrc_actuator`：执行器产生的广义力；`qfrc_applied`：外部施加的广义力；`qfrc_constraint`：约束力。
- 练习中通过 `J.T @ F_ext` 注入**已知的理想外力扰动**，不代表已经建立了真实的物理接触或使用了力传感器。

## 7. 本日最重要的认识

1. 位置控制关注跟踪误差；阻抗控制关注**外力—运动之间的动态关系**。
2. 重力补偿 PD 已能表现出一定柔顺性，但完整阻抗控制进一步支持**期望惯量整形**。
3. 任务空间刚度、阻尼和惯量不能直接由关节空间矩阵照搬；必须考虑 Jacobian 和动力学映射。
4. 无源性从能量交换的角度提供分析机器人—环境交互的重要工具，但不能与渐近稳定性画等号。
5. 目前实验是**无实际接触的已知扰动力矩/外力仿真**；接触约束、摩擦、接触力测量、控制时延与接触稳定性仍是后续内容。

## 8. 下一阶段预告

**Day 27：浮动基座动力学（Floating-Base Dynamics）入门**。

从固定基座机械臂转向腿足/人形机器人：

$$M(q)\ddot q+h(q,\dot q)=S^T\tau+J_c(q)^T\lambda.$$

重点理解：自由基座的 6 个非驱动自由度、驱动选择矩阵 $S$、接触 Jacobian $J_c$、接触力 $\lambda$，以及为何足端接触能改变机器人可实现的运动。之后再进入实际接触与约束动力学实验。

---

**完成情况**：Day 26 理论学习与三组 MuJoCo 实验均已完成。力控制和混合运动/力控制已建立概念基础，真实接触实验留待后续。
