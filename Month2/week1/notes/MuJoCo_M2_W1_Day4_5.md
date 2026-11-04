# ROBO_Learn｜MuJoCo 练习总结（Day 4.5）

> 阶段定位：机器人动力学基础之后、Day 5 机器人平衡与 ZMP 之前。  
> 本文记录已练习的 MuJoCo 建模、控制、运动学和 Python 接口，以及已发现的问题；接触力、外力、Sensor 等留待接触动力学仿真前集中学习。

## 1. 已完成的练习

### 1.1 单摆：MJCF + PD + 重力补偿

- 使用 MJCF 定义 `worldbody`、`body`、`joint`、`geom` 和 `actuator/motor`。
- 通过 `data.qpos`、`data.qvel` 读取状态，构造关节 PD 控制律。
- 观察到纯 PD 控制存在由重力引起的稳态误差。
- 在练习所用单摆参数和角度约定下加入重力补偿项：

  \[
  \tau_g(q)=9.81\times0.5\sin q
  \]

  实验中消除了观察到的稳态偏差。

**关键认识**：PD 不一定能在恒定外扰（包括目标位形下的重力矩）存在时实现零稳态误差；补偿项的符号必须与模型坐标、控制力矩方向一致。

### 1.2 双连杆：独立关节 PD 控制

- 定义嵌套 Body、两个 hinge Joint 和 Capsule Geom。
- 为两个关节设置独立 PD 控制器。
- 记录关节角、角速度和力矩，使用 NumPy/Matplotlib 绘图。
- 明确区分“定义初始值变量”与“真正设置仿真初始状态”：

```python
data.qpos[:] = qpos_init
data.qvel[:] = qvel_init
mujoco.mj_forward(model, data)
```

- 记录历史数据时使用 `.copy()`，避免保存可随仿真更新的数组视图：

```python
q_history.append(data.qpos.copy())
q_history = np.array(q_history)  # (T, nq)
```

### 1.3 Viewer：可视化和运行管理

```python
import mujoco.viewer

with mujoco.viewer.launch_passive(model, data) as viewer:
    while viewer.is_running():
        # 在此处设置控制输入
        mujoco.mj_step(model, data)
        viewer.sync()
```

- `launch_passive`：打开被动 Viewer。
- `viewer.sync()`：同步仿真数据与界面。
- `viewer.is_running()`：判断窗口是否仍在运行。
- `with ... as ...`：Python 上下文管理器，退出时由对象的上下文管理协议处理资源清理；不是简单地自动调用任意对象的 `close()`。

### 1.4 独立运动学练习：Body / Site / Jacobian

模型：两个绕局部 Y 轴转动的 hinge 关节，连杆沿局部负 Z 轴延伸；长度分别为 `1.0 m` 和 `0.8 m`；末端 Site 名为 `ee`。

测试位形和速度：

\[
q=[0.5,-0.3]^T,\qquad \dot q=[0.2,-0.4]^T.
\]

读取结果：

```text
link1 xpos: [ 0.          0.          0.        ]
link2 xpos: [-0.47942554  0.         -0.87758256]
ee xpos:    [-0.638361    0.         -1.66163582]

link2 xquat: [0.99500417 0. 0.09983342 0.]
```

Jacobian 计算：

```python
jacp = np.zeros((3, model.nv))
jacr = np.zeros((3, model.nv))
mujoco.mj_jacSite(model, data, jacp, jacr, model.site("ee").id)

v_ee = jacp @ data.qvel
omega_ee = jacr @ data.qvel
```

对应结果：

```text
jacp:
[[-1.66163582 -0.78405326]
 [ 0.          0.        ]
 [ 0.638361    0.15893546]]

jacr:
[[0. 0.]
 [1. 1.]
 [0. 0.]]

v_ee: [-0.01870586  0.  0.06409801]
```

**物理解释**：末端运动被约束在 XZ 平面，因此线速度 Jacobian 的 Y 行为零；两个关节均绕 Y 轴，角速度 Jacobian 仅 Y 行非零。

## 2. MJCF 的关键概念

| MJCF 元素 | 用途 | 重要区别 |
|---|---|---|
| `<worldbody>` | 世界中的刚体和几何体容器 | 世界坐标系为根 |
| `<body pos="...">` | 定义 Body 相对父 Body 的坐标原点 | `pos` 是局部位置 |
| `<joint type="hinge" axis="...">` | 定义旋转自由度 | 轴方向按局部坐标表达 |
| `<geom type="capsule" fromto="..." size="...">` | 定义几何体 | `fromto` 决定轴线端点，Capsule 的 `size` 首项为半径 |
| `<site pos="...">` | 定义附着在 Body 上的参考点/坐标系 | Site 随所属 Body 运动 |
| `<actuator><motor .../></actuator>` | 定义驱动输入 | 与关节、传动关系关联 |

**容易混淆的地方**：本练习 Capsule 的 `size="0.04"` 是半径，不是连杆长度。连杆长度来自 Body 层级、Site 位置及相关几何定义。

## 3. Model 与 Data：常用读取接口

`model` 保存模型定义；`data` 保存当前状态和计算结果。

| 接口 | 含义 |
|---|---|
| `model.nq`, `model.nv` | 广义位置、广义速度维数 |
| `data.qpos`, `data.qvel` | 当前广义位置、速度 |
| `model.body("link2").pos` | Body 相对父 Body 的局部平移 |
| `data.body("link2").xpos` | Body 坐标系原点的世界位置 |
| `data.body("link2").xipos` | Body 质心的世界位置 |
| `data.body("link2").xmat` | Body 局部系到世界系的旋转矩阵 |
| `data.body("link2").xquat` | Body 世界姿态四元数，顺序 `[w,x,y,z]` |
| `model.site("ee").pos` | Site 在所属 Body 局部系的位置 |
| `data.site("ee").xpos` | Site 的世界位置 |
| `data.site("ee").xmat` | Site 的世界姿态旋转矩阵 |
| `model.site("ee").id` | Site ID |
| `data.qfrc_bias` | 动力学偏置广义力 |
| `data.qfrc_actuator` | 驱动产生的广义力 |
| `data.qfrc_applied` | 用户施加的广义力 |

还可使用 `mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "link2")` 查询 ID，并以 `data.xpos[body_id]` 访问相应数据。

### 局部坐标与世界坐标

对于 `link2` 上的 Site：

\[
{}^Wp_{ee}={}^Wp_2+{}^WR_2\,{}^2p_{ee}.
\]

在本练习中：

\[
{}^2p_{ee}=[0,0,-0.8]^T,\qquad {}^WR_2=R_y(q_1+q_2).
\]

Site 的局部 `pos` 固定不变，但世界 `xpos` 会随关节运动变化。Site 的偏移由 MJCF 指定，不是由 Geom 形状自动决定。

## 4. 正运动学与 Jacobian

两杆长度 \(l_1=1,\ l_2=0.8\)，定义：

\[
p_{ee}=
\begin{bmatrix}
-l_1\sin q_1-l_2\sin(q_1+q_2)\\
0\\
-l_1\cos q_1-l_2\cos(q_1+q_2)
\end{bmatrix}.
\]

位置 Jacobian：

\[
J_v=
\begin{bmatrix}
-l_1\cos q_1-l_2\cos(q_1+q_2)&-l_2\cos(q_1+q_2)\\
0&0\\
l_1\sin q_1+l_2\sin(q_1+q_2)&l_2\sin(q_1+q_2)
\end{bmatrix}.
\]

角速度 Jacobian：

\[
J_\omega=
\begin{bmatrix}0&0\\1&1\\0&0\end{bmatrix}.
\]

于是：

\[
v_{ee}=J_v\dot q,\qquad \omega_{ee}=J_\omega\dot q.
\]

`mj_jacSite` 的 `jacp` 与 `jacr` 分别对应世界坐标系下 Site 原点的线速度 Jacobian、角速度 Jacobian。它们不能不加区分地直接视为《Modern Robotics》中按 `[ω;v]` 排列的 Space Twist Jacobian。

## 5. `mj_forward`、`mj_step` 与数值验证

- `mujoco.mj_forward(model, data)`：在当前状态下重新计算运动学、动力学等派生量，**不推进时间**。
- `mujoco.mj_step(model, data)`：推进仿真一步。
- 仅验证位置对关节角的导数时，使用 `mj_forward` 即可，无须进行动力学积分。

### 有限差分验证末端线速度

\[
v_{fd}\approx\frac{p(q+\dot q\Delta t)-p(q)}{\Delta t}.
\]

```python
data.qpos[:] = q
data.qvel[:] = dq
mujoco.mj_forward(model, data)
p_now = data.site("ee").xpos.copy()

dt = 1e-6
data.qpos[:] = q + dq * dt
mujoco.mj_forward(model, data)
p_next = data.site("ee").xpos.copy()

v_fd = (p_next - p_now) / dt
print("finite difference:", v_fd)
print("Jacobian:", jacp @ dq)  # jacp 应在原始 q 处计算
```

注意：如果在更新状态之后才读取 `p_now`，就会错误地计算 `(p_next-p_next)/dt=0`。

## 6. 练习中发现的两个问题

**问题 A：正运动学长度取错。**

曾将 `size="0.04"` 误用为两杆长度，导致理论末端位置与仿真值相差约 `1.70095 m`。应使用 `l1=1.0, l2=0.8`。

**问题 B：有限差分旧位置未保存。**

曾在 `mj_forward` 更新状态后，用 `p_next - data.site("ee").xpos` 计算速度，结果必然为零。应在修改 `qpos` **之前**执行 `p_now = ...xpos.copy()`。

这两处已定位并给出修正方法；本总结不声称已再次运行修正后的程序。

## 7. 下一阶段安排

**立即进入 Day 5（先理论、后仿真）：**

1. 静态平衡、支撑多边形。
2. CoP 与 ZMP 的定义及其区别。
3. LIPM 的假设、推导、稳定性与控制。
4. NumPy/Python 中的 LIPM 控制实验。

**暂缓至接触动力学仿真前集中学习的 MuJoCo 内容：**

- `data.contact`、`data.ncon` 与 `mj_contactForce()`。
- `data.qfrc_applied`、`data.xfrc_applied` 的外力施加。
- Sensor 的 MJCF 定义及 `data.sensordata`。
- Actuator 的传动、控制输入和驱动力接口（在此前电机基础上扩展）。
- 接触模型、地面反作用力、CoP/ZMP 数据提取与调试。

> 学习原则：MuJoCo 接口按实际任务需求学习；先建立理论和控制方法，再用仿真验证。
