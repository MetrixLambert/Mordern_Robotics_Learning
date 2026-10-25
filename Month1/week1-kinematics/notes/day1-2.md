# ROBO_Learn 学习笔记：Day 1–2

**教材**：Lynch & Park, *Modern Robotics: Mechanics, Planning, and Control*，第 2 章（构型空间）与第 3 章（刚体旋转）  
**学习日期**：2026-10-11  
**学习范围**：Day 1 构型空间与旋转矩阵；Day 2 旋转复合与角速度。  
**下一单元**：Day 3，SO(3) 指数映射、Rodrigues 公式及对数映射。

> 约定：坐标系 `{s}` 为世界系，`{b}` 为刚体系；列向量记法。`R_sb` 将向量的 **b 系坐标**变为 **s 系坐标**：`v_s = R_sb v_b`。旋转矩阵均属于 `SO(3)`。

---

## Day 1｜构型空间与旋转矩阵

### 1. 构型、自由度与构型空间

- **构型（configuration）**：完整描述机器人或刚体几何状态所需的信息。
- **自由度（DOF）**：描述构型所需的独立参数数目。
- 三维空间中自由刚体有 **6 个自由度**：3 个平移 + 3 个旋转。
- 自由刚体构型属于 `SE(3)`；其姿态属于 `SO(3)`。直观上，构型由位置 `p ∈ ℝ³` 与姿态 `R ∈ SO(3)` 共同决定，但 **`SO(3)` 不是普通的 `ℝ³` 向量空间**。

### 2. 旋转矩阵的几何意义

设 `R_ab` 表示坐标系 `{b}` 相对于 `{a}` 的姿态，则

\[
\boxed{v_a=R_{ab}v_b},\qquad \boxed{v_b=R_{ab}^{T}v_a}.
\]

`R_ab` 的三列分别是 **b 系 x、y、z 单位轴在 a 系中的坐标**：

\[
R_{ab}=\begin{bmatrix}\hat x_b^a&\hat y_b^a&\hat z_b^a\end{bmatrix}.
\]

**旋转矩阵的性质**：

\[
\boxed{R^TR=RR^T=I},\qquad \boxed{\det R=+1},\qquad R^{-1}=R^T.
\]

`det = +1` 表示保持右手系；`det = -1` 的正交矩阵还包含反射，不属于 `SO(3)`。

### 3. 主动旋转与被动坐标变换

- **主动旋转**：几何向量本身发生旋转，坐标系保持不变。例如同一坐标系内 `v' = Qv`。
- **被动坐标变换**：几何向量不变，只改变其坐标表达。例如 `v_a = R_ab v_b`。
- 两者可以出现相同数值的矩阵，但**物理含义、下标和方向必须明确**；不能只看一个 `R` 就判断是在主动旋转还是换坐标。

绕固定 z 轴的正向旋转矩阵：

\[
\boxed{R_z(\theta)=\begin{bmatrix}
\cos\theta&-\sin\theta&0\\
\sin\theta&\cos\theta&0\\
0&0&1
\end{bmatrix}}.
\]

右手定则：绕 +z 旋转 +90°，+x 方向转向 +y。

### 4. Day 1 易错点

1. `R_ab` 的**列**是 b 轴在 a 中的坐标；不是 a 轴在 b 中的坐标。
2. `v_a=R_ab v_b` 与 `v_b=R_ab^T v_a` 不可混用。
3. 不要把旋转矩阵的 9 个元素当成 9 个独立自由度：正交性与行列式约束使姿态只有 3 个自由度。

---

## Day 2｜旋转复合、左右乘与角速度

### 1. 连续坐标变换

\[
\boxed{R_{ac}=R_{ab}R_{bc}}
\]

因为 `v_a = R_ab v_b = R_ab R_bc v_c`。矩阵乘法从**右向左**作用，顺序不能随意交换。

### 2. 世界系旋转与刚体系旋转

令 `R = R_sb` 表示当前刚体姿态。对姿态施加旋转增量：

| 旋转增量描述在哪个坐标系 | 姿态更新公式 | 直观理解 |
|---|---|---|
| 世界系 `{s}` | `R_new = Q_s R_old` | **左乘**：在世界系中旋转当前刚体各轴 |
| 当前刚体系 `{b}` | `R_new = R_old Q_b` | **右乘**：绕刚体当前自身轴旋转 |

同一次物理旋转在两坐标系下的算子满足：

\[
\boxed{Q_s=RQ_bR^T},\qquad \boxed{Q_sR=RQ_b}.
\]

**推导思路**：对于任意向量 `v_s`，先变换到 b 系 (`R^T v_s`)，在 b 系执行旋转 (`Q_b R^T v_s`)，再变回 s 系 (`R Q_b R^T v_s`)；与直接在 s 系旋转 (`Q_s v_s`) 相同。

**具体例子**：当 `R = R_z(90°)` 时，刚体 +y 轴对应世界 −x 轴，因此

\[
R_x(-30^\circ)R_z(90^\circ)=R_z(90^\circ)R_y(30^\circ).
\]

注意：`Q_s` 与 `Q_b` 是**同一物理旋转增量的不同坐标表达**，不能在不改变旋转轴定义的情况下直接互换。

### 3. 三维角速度的物理意义

\[
\boxed{\omega=\dot\theta\,\hat\omega}
\]

- `\hat ω`：瞬时旋转轴方向（右手定则）；`\dot θ`：有符号旋转速率（rad/s）。
- 对绕固定点 `O` 纯旋转的刚体，点 `P` 相对 `O` 的位置向量为 `r_OP`，则其速度为

\[
\boxed{v_P=\omega\times r_{OP}}.
\]

- 若参考点 `O` 同时平移，则应写为 `v_P=v_O+ω×r_OP`。
- 角速度是**几何向量**，一般**不等于三个欧拉角导数直接组成的向量**。

### 4. 反对称矩阵：skew 与 vee

定义叉乘矩阵（hat/skew 算子）：

\[
\boxed{[\omega]_\times=\begin{bmatrix}
0&-\omega_z&\omega_y\\
\omega_z&0&-\omega_x\\
-\omega_y&\omega_x&0
\end{bmatrix}},\qquad
\boxed{[\omega]_\times v=\omega\times v}.
\]

反向映射 `vee`：若 `W=[ω]_×`，则

\[
\boxed{\operatorname{vee}(W)=\begin{bmatrix}W_{32}&W_{13}&W_{21}\end{bmatrix}^T}.
\]

其中 `W_ij` 用数学上的 **1-based** 下标；NumPy 中对应 `W[2,1]`、`W[0,2]`、`W[1,0]`。

### 5. 角速度与姿态矩阵导数：本节核心

世界系角速度（spatial）与刚体系角速度（body）分别满足：

\[
\boxed{\dot R=[\omega_s]_\times R},\qquad
\boxed{\dot R=R[\omega_b]_\times}.
\]

反向提取：

\[
\boxed{[\omega_s]_\times=\dot R R^T},\qquad
\boxed{[\omega_b]_\times=R^T\dot R}.
\]

**推导主线**：

1. 取固定在刚体上的向量 `r_b`，所以 `\dot r_b=0`，但世界坐标 `r_s=Rr_b` 会变化。
2. 对 `r_s=Rr_b` 求导：`\dot r_s=\dot R r_b`。
3. 由瞬时旋转运动：`\dot r_s=ω_s×r_s=[ω_s]_×Rr_b`。
4. 对任意 `r_b` 成立，故 `\dot R=[ω_s]_×R`。

**正交性带来的结构**：对 `R^TR=I` 求导，

\[
\dot R^TR+R^T\dot R=0
\quad\Rightarrow\quad
(R^T\dot R)^T=-R^T\dot R.
\]

因此 `R^T\dot R` 是反对称矩阵，恰好可以用一个三维角速度向量表示。

### 6. 同一个角速度在两个坐标系下的表示

\[
\boxed{\omega_s=R\omega_b},\qquad
\boxed{\omega_b=R^T\omega_s}.
\]

角速度的模长不变：`||ω_s|| = ||ω_b||`。叉乘矩阵对应地满足：

\[
\boxed{[\omega_s]_\times=R[\omega_b]_\times R^T}.
\]

这是旋转保持叉乘的结果：`R(a×b)=(Ra)×(Rb)`，其中 `R∈SO(3)`。

**何时 `ω_s = ω_b`？**

\[
(R-I)\omega_b=0.
\]

- `ω_b=0` 时总成立。
- `R=I` 时对所有角速度都成立。
- 对 `R≠I` 的三维旋转，非零角速度须沿 `R` 的固定旋转轴（即 `ker(R−I)`）方向，才有相同坐标。任意 `SO(3)` 旋转都有特征值 1。
- 例：`R=R_z(θ)` 且 `ω` 沿 z 轴时，`ω_s=ω_b=[0,0,ω]^T`。

### 7. Day 2 已完成的 NumPy 练习

实现了 `rot_z(theta)`、`skew(w)`、`vee(W)`；用中心差分近似姿态导数：

\[
\dot R(t)\approx\frac{R(t+h)-R(t-h)}{2h},\qquad h=10^{-5}.
\]

对 `R(t)=R_z(2t)`、`t=0.5 s`，分别计算

```python
Omega_s = R_dot @ R.T
Omega_b = R.T @ R_dot
omega_s = vee(Omega_s)
omega_b = vee(Omega_b)
```

数值结果均约为 `[0, 0, 2]` rad/s，与理论一致。注意：有限差分得到的矩阵不一定**严格**反对称（可能有很小数值误差）；需要时可取 `(Omega - Omega.T)/2` 后再 `vee`。

---

## 两天速记卡

| 问题 | 必须能快速写出的答案 |
|---|---|
| `R_ab` 的列表示什么？ | b 系三个单位轴在 a 系中的坐标 |
| 坐标如何变换？ | `v_a=R_ab v_b`；反向用转置 |
| 旋转矩阵满足什么约束？ | `R^TR=I`，`det R=1` |
| 连续坐标变换？ | `R_ac=R_ab R_bc` |
| 世界系增量与刚体系增量？ | `Q_s R`（左乘）；`R Q_b`（右乘） |
| 同一旋转算子换坐标系？ | `Q_s=R Q_b R^T` |
| 角速度对应的点速度？ | `v=ω×r`（绕固定点纯旋转） |
| skew 的作用？ | `[ω]_×v=ω×v` |
| 从 `R, R_dot` 提取世界系角速度？ | `[ω_s]_×=R_dot R^T` |
| 提取刚体系角速度？ | `[ω_b]_×=R^T R_dot` |
| 两系角速度如何转换？ | `ω_s=R ω_b` |

## 自测题（先闭卷，再展开答案）

1. 若 `R_sb=R_z(90°)`，b 系的 +x、+y 轴在 s 系分别指向哪里？
2. 世界系绕 +x 轴旋转 `+30°`，当当前姿态为 `R_z(90°)` 时，等效的刚体系旋转轴和角度是什么？
3. 为什么 `R^T R_dot` 必定是反对称矩阵？
4. 若 `R=R_z(90°)`、`ω_b=[3,0,0]^T`，求 `ω_s`。
5. 已知 `R(t)` 与 `R_dot(t)`，写出求 `ω_s` 和 `ω_b` 的完整步骤。

<details>
<summary>参考答案（完成后查看）</summary>

1. b 系 +x 指向 s 系 +y；b 系 +y 指向 s 系 −x。
2. 世界 +x 等于刚体 −y，因此等效为绕刚体 −y 轴 `+30°`，或绕刚体 +y 轴 `−30°`。
3. 对 `R^TR=I` 求导：`R_dot^T R + R^T R_dot=0`，前一项正是后一项的转置。
4. `ω_s=R ω_b=[0,3,0]^T`。
5. 计算 `Omega_s=R_dot R^T`、`Omega_b=R^T R_dot`，再分别应用 `vee`；理想连续情形两者严格反对称。

</details>

## 学习进度记录

- **Day 1**：完成。构型空间、旋转矩阵几何意义、正交性、主动/被动旋转及 NumPy 验证。
- **Day 2**：完成。旋转复合、世界/刚体系左右乘、角速度、skew/vee、`R_dot` 与 `ω` 的关系及 NumPy 数值验证。
- **当前掌握程度**：Day 1 约 2–3/3；Day 2 约 2.5/3。后续通过一般三维旋转的例子巩固。
- **下次学习**：Day 3，先闭卷复习 3 道 Day 2 问题，再学习 `SO(3)` 指数与对数映射。
