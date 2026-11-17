# ROBO_Learn · Month 2 Day 19 学习总结

**主题**：二维落脚点规划、步态时序、周期 DCM 与反馈修正  
**状态**：已完成；依据用户提交的代码运行结果及图像（2026-11-17）

## 1. 与 Day 18 的衔接

Day 18 讨论单次捕获：触地时若令 $p_{next}=\xi(T)$，理想 LIPM 可渐近停止。连续行走则需要使每步触地后 DCM 相对支撑脚保持设计的偏移量 $b$，而不是每步都停下。

固定 ZMP 的 DCM 动力学及解析解：

$$\dot{\boldsymbol\xi}=\omega(\boldsymbol\xi-\boldsymbol p),\quad\omega=\sqrt{g/z_c},\quad \boldsymbol\xi(t)=\boldsymbol p+[\boldsymbol\xi(0)-\boldsymbol p]e^{\omega t}.$$

## 2. 几何足印与时序

- 世界坐标系：$x$ 朝前，$y$ 朝左。
- 左脚 $y=+W/2$，右脚 $y=-W/2$；相邻落脚点前向间隔 $L$。
- 简化步态：固定单支撑时间 $T$、双支撑时间 $0$，支撑脚交替切换。
- 步态索引 $k=\lfloor t/T\rfloor$，局部归一化相位 $s=(t-kT)/T\in[0,1)$。
- 这一阶段的时间驱动切换是理想化的：真实触地应检测接触事件，不能只依赖计时。

## 3. DCM 反馈落脚点修正

若当前实际 DCM 与参考 DCM 的误差为 $e_\xi=\xi_{actual}-\xi_{nom}$，剩余支撑时间为 $T_r$，且比较的两条轨迹在触地前使用相同 ZMP，则：

$$\Delta\xi(T_r)=e^{\omega T_r}e_\xi,\qquad \boxed{p_{corrected}=p_{nom}+e^{\omega T_r}e_\xi}.$$

实验中 $z_c=1$ m、$T_r=0.5$ s、$e_\xi=(0.05,0)$ m，则前向修正约 $0.239$ m；名义落脚点 $(0.30,-0.10)$ m 修正为约 $(0.539,-0.10)$ m。

该方法尚未考虑腿的可达范围、摆动速度、脚掌面积或接触约束。过大的修正不能直接执行。

## 4. 周期 DCM 的两轴偏移

令 $a=e^{\omega T}$。一维等步长前进的周期偏移：

$$\boxed{b_x=\frac{L}{a-1}}.$$

横向原地交替支撑，$p_L=+W/2$、$p_R=-W/2$，由周期条件可得：

$$\boxed{b_L=-\frac{W}{1+a}},\qquad\boxed{b_R=+\frac{W}{1+a}}.$$

即第 $k$ 步起始 DCM 参考可设为：

$$\boldsymbol\xi_k=\boldsymbol p_k+\begin{bmatrix}L/(a-1)\\(-1)^{k+1}W/(a+1)\end{bmatrix},\quad \boldsymbol p_k=\begin{bmatrix}kL\\(-1)^kW/2\end{bmatrix}.$$

在每个单支撑区间，使用固定当前支撑点 $\boldsymbol p_k$ 的 DCM 解析解生成连续参考。理想切换时 DCM 连续、ZMP 可跳变。

## 5. 已完成实验和验证

1. 编写 `generate_footsteps`：输出 8 个左右交替的二维足印。
2. 编写 `generate_gait_timing`：输出步态时间节点。
3. 编写 `predict_dcm`：起初将解析解首项误写成 `xi0`，更正为 `p` 后测试通过。额外可用 $T=0$ 的初值条件检查公式。
4. 编写 `correct_footstep`：实现 DCM 反馈修正，通过两条轨迹相对落脚点偏移相等的断言。
5. 进一步生成连续二维 DCM 轨迹；图示表现为沿前向前进、横向左右周期摆动，切换处无明显跳变。

## 6. 对后续控制架构的认识

- 周期 DCM 轨迹提供**动态平衡参考**，不是实际质心轨迹：$\xi=c+\dot c/\omega$。
- 高层 DCM 控制可根据误差调整可行 ZMP 或下一步落脚点；WBC 通常进一步落实质心/接触力/关节任务，并非必须直接跟踪 DCM。
- ZMP/CoP 对水平平面接触通常受当前支撑多边形约束；ZMP 落在支撑域内并不自动保证动态稳定。
- 实践瓶颈：支撑切换、脚掌接触和运动学可达性。Day 20 转向摆动足轨迹与 IK。

**验收结论**：Day 19 通过。尚未进行真实 MuJoCo 全身运动或支撑接触仿真，不应把当前足印/DCM 图解释为成功行走。
