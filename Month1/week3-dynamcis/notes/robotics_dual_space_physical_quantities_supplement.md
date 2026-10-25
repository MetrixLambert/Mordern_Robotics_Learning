# ROBO_Learn 补充笔记：物理量、参考系与对偶空间

> 适用范围：Modern Robotics 第 3、5、8 章；对应 Day 8（静力学）、Day 11（刚体动力学）、Day 12（RNEA）中的讨论。  
> 目的：把 **Twist、Wrench、广义速度/广义力、对偶空间、Adjoint 及其转置** 放在统一框架中理解。  
> 全文采用 Modern Robotics 的 **角速度/力矩在前** 约定。

---

## 1. 先分清：物理对象、坐标分量、参考点

机器人学中经常出现三个不同层次：

1. **物理对象**：同一个刚体的运动、同一个力系，不因更换坐标系而改变。
2. **坐标表示**：选定坐标系后得到的数字向量；更换坐标系，数字通常改变。
3. **参考点**：角速度与合力不依赖所选原点，但线速度和力矩往往与参考点有关。

例如，作用于点 $P$ 的力 $\boldsymbol f$，相对于原点 $O$ 的力矩为

$$
\boldsymbol m_O=(\boldsymbol p_P-\boldsymbol p_O)\times\boldsymbol f.
$$

如果把原点换到 $O'$，同一个力对应的力矩会改变：

$$
\boldsymbol m_{O'}=\boldsymbol m_O-(\boldsymbol p_{O'}-\boldsymbol p_O)\times\boldsymbol f.
$$

**重要**：改变力矩参考点，不意味着实际施加的力发生了变化。

---

## 2. Twist：刚体的瞬时运动

### 2.1 定义与排列顺序

$$
V=\begin{bmatrix}\boldsymbol\omega\\\boldsymbol v\end{bmatrix}\in\mathbb R^6,
$$

其中 $\boldsymbol\omega$ 是角速度分量，$\boldsymbol v$ 是 Twist 的线性分量。

给定 $T_{sb}=\begin{bmatrix}R&p\\0&1\end{bmatrix}$，其中 $p$ 为 Body 原点在 Space 系中的位置：

$$
[V_s]=\dot T_{sb}T_{sb}^{-1},
\qquad
[V_b]=T_{sb}^{-1}\dot T_{sb}.
$$

展开：

$$
\boldsymbol\omega_s=R\boldsymbol\omega_b,\qquad
\boldsymbol v_s=\dot p-\boldsymbol\omega_s\times p,
$$

$$
\boldsymbol v_b=R^T\dot p.
$$

**易错点**：一般而言，Space Twist 的 $\boldsymbol v_s$ **不是** Body 原点的实际线速度 $\dot p$；两者相差 $\boldsymbol\omega_s\times p$。

对 Space 系中位置为 $x$ 的刚体上某点，其瞬时速度为

$$
\dot x=\boldsymbol\omega_s\times x+\boldsymbol v_s.
$$

取 $x=p$，便得到 $\dot p=\boldsymbol\omega_s\times p+\boldsymbol v_s$。

### 2.2 Adjoint：同一个 Twist 的两种表示

$$
V_s=\operatorname{Ad}_{T_{sb}}V_b,
$$

$$
\operatorname{Ad}_{T}=
\begin{bmatrix}
R&0\\
[p]_\times R&R
\end{bmatrix}.
$$

逆变换：

$$
V_b=\operatorname{Ad}_{T^{-1}}V_s,\qquad
\operatorname{Ad}_{T^{-1}}=(\operatorname{Ad}_{T})^{-1}.
$$

这里的 Adjoint 同时处理旋转和原点变化，不能仅用 $\operatorname{diag}(R,R)$ 替代。

---

## 3. Wrench：与 Twist 配对的力学量

### 3.1 定义

$$
F=\begin{bmatrix}\boldsymbol m\\\boldsymbol f\end{bmatrix}\in\mathbb R^6.
$$

其中 $\boldsymbol m$ 是关于**所用坐标系原点**的合力矩，$\boldsymbol f$ 是合力。

Twist 与 Wrench 之间最重要的关系是**瞬时功率**：

$$
\boxed{P=F^TV=\boldsymbol m^T\boldsymbol\omega+\boldsymbol f^T\boldsymbol v.}
$$

这要求 $F$ 和 $V$ 是**匹配的坐标表示**：例如都是 Space 表示，或者都是 Body 表示。

注意 Twist 的 $\boldsymbol v_s$ 不是末端原点速度时，仍可以使用 $F_s^TV_s$，但必须把 $\boldsymbol m_s$ 定义为关于 Space 原点的力矩。

### 3.2 为什么 Wrench 不能像 Twist 一样变换？

令 $A=\operatorname{Ad}_T$。Twist 满足

$$
V_s=AV_b.
$$

同一物理运动与力系产生的功率不应依赖坐标系：

$$
F_s^TV_s=F_b^TV_b.
$$

代入 $V_s=AV_b$：

$$
F_s^TAV_b=F_b^TV_b.
$$

由于对任意 $V_b$ 成立，

$$
A^TF_s=F_b,
$$

所以

$$
\boxed{F_s=A^{-T}F_b=\operatorname{Ad}_T^{-T}F_b.}
$$

反向写为

$$
\boxed{F_b=\operatorname{Ad}_T^TF_s.}
$$

**记忆方法**：速度用 $A$，与其做功率配对的力用 $A^{-T}$。不要把两个方向的公式混用。

展开可得

$$
\boldsymbol f_s=R\boldsymbol f_b,\qquad
\boldsymbol m_s=R\boldsymbol m_b+p\times(R\boldsymbol f_b).
$$

后一项反映了力矩参考点的变化。

---

## 4. 什么是“对偶空间”？

### 4.1 向量与线性泛函

设 $X$ 为向量空间。对偶空间 $X^*$ 是从 $X$ 到实数的**线性映射**的集合：

$$
X^*=\{\ell:X\to\mathbb R\mid \ell\text{ 是线性的}\}.
$$

在有限维情况下，选定一组基以后，可以把 $X$ 的向量表示为列向量 $x$，把 $X^*$ 中的线性泛函表示为行向量 $\ell^T$，于是

$$
\ell(x)=\ell^Tx.
$$

虽然欧氏内积可以帮助我们将向量与对偶向量相联系，但二者在**坐标变换规则**上扮演不同角色。

- Twist 是“运动方向”，可视作速度空间中的向量。
- Wrench 是“给定一个运动方向，输出对应功率”的线性泛函，因此属于速度空间的对偶空间。
- 关节速度 $\dot q$ 与广义力 $\tau$ 也是一对相应的量：$\tau^T\dot q$ 为功率。

更严谨地说，在构型 $q$ 处，$\dot q\in T_qQ$，广义力 $\tau\in T_q^*Q$；机器人运动速度与力的配对也可按切空间与余切空间理解。

### 4.2 二维例子：为什么一定是逆转置？

取

$$
x=\begin{bmatrix}1\\2\end{bmatrix},\qquad
f=\begin{bmatrix}3\\4\end{bmatrix},\qquad
A=\begin{bmatrix}2&0\\0&1\end{bmatrix}.
$$

原始配对为

$$
f^Tx=3+8=11.
$$

若向量换表示为

$$
x'=Ax=\begin{bmatrix}2\\2\end{bmatrix},
$$

为了让同一个线性泛函的数值不变，必须取

$$
f'=A^{-T}f=\begin{bmatrix}1.5\\4\end{bmatrix}.
$$

验证：

$$
(f')^Tx'=1.5\times2+4\times2=11.
$$

若误用 $f'=Af$，就无法保持这个配对不变。

> 关键直觉：**对偶向量不是“另一种普通速度向量”，而是对速度向量进行线性测量的对象。**

### 4.3 一个需要澄清的细节

有时我们把 Wrench 也存储成 NumPy 的 `(6,)` 数组，看起来与 Twist 没有区别。这是**存储方式相同**，并不代表它们的变换规律相同。

同样，`F @ V` 在 NumPy 中是内积，数学意义仍然是 $F^TV$ 这一对偶配对。

---

## 5. Jacobian 转置：对偶映射如何出现？

关节速度通过 Jacobian 映射到末端 Twist：

$$
V=J(q)\dot q.
$$

如果末端受到外部 Wrench $F_{\mathrm{ext}}$，它对机械臂做功的瞬时功率为

$$
P=F_{\mathrm{ext}}^TV
=F_{\mathrm{ext}}^TJ\dot q
=(J^TF_{\mathrm{ext}})^T\dot q.
$$

另一方面，外力对应的关节广义力 $\tau_{\mathrm{ext}}$ 满足

$$
P=\tau_{\mathrm{ext}}^T\dot q.
$$

因此

$$
\boxed{\tau_{\mathrm{ext}}=J^TF_{\mathrm{ext}}.}
$$

这不是人为规定，而是由**功率守恒/虚功关系**得到的对偶映射。

注意符号：

- $J^TF_{\mathrm{ext}}$：**外力对关节产生的广义力**。
- 若只考虑该外力、且要求静力平衡，电机需要施加的补偿力矩为 $-J^TF_{\mathrm{ext}}$。

对于位置 Jacobian，若 $\dot p=J_p\dot q$ 且只考虑作用于末端点的力 $\boldsymbol f$，则

$$
\tau_{\mathrm{ext}}=J_p^T\boldsymbol f.
$$

若用 Space Jacobian，需把该力转换成关于 Space 原点的完整 Wrench：

$$
F_s=\begin{bmatrix}p\times\boldsymbol f\\\boldsymbol f\end{bmatrix},
\qquad
J_s^TF_s=J_p^T\boldsymbol f.
$$

**不要**把“末端点的纯力”简单写为 $[0;\boldsymbol f]$ 后直接与任意 Space Jacobian 配对；其力矩部分取决于参考原点。

---

## 6. RNEA 中的两种映射：为什么一个是 Adjoint，另一个是转置？

令 $T_{i,i+1}={}^{i}T_{i+1}$ 表示子连杆坐标系 $\{i+1\}$ 相对于父连杆坐标系 $\{i\}$ 的位姿，并定义

$$
H_{i+1}=\operatorname{Ad}_{T_{i,i+1}^{-1}}.
$$

### 6.1 前向：速度从父系传到子系

$$
\boxed{
V_{i+1}=H_{i+1}V_i+A_{i+1}\dot q_{i+1}.
}
$$

这里 $H_{i+1}V_i$ 把父连杆的运动表示到子连杆坐标系。

### 6.2 反向：载荷从子系传到父系

忽略外部载荷时：

$$
\boxed{
F_i=G_i\dot V_i-\operatorname{ad}_{V_i}^TG_iV_i
+H_{i+1}^TF_{i+1}.
}
$$

因为 $H_{i+1}$ 把**父系速度表示**映射到**子系速度表示**，所以其对偶映射 $H_{i+1}^T$ 恰好把**子系 Wrench**传回**父系**。

也可以直接由功率看出：

$$
F_{i+1}^T(H_{i+1}V_i)
=(H_{i+1}^TF_{i+1})^TV_i.
$$

### 6.3 今天遇到的数组索引错误

若 `T[k]` 存的是第 `k` 根连杆相对其父连杆的变换（Python 采用 0-based 索引），那么：

- 前向计算第 `i` 根连杆时，用 `T[i]`；
- 反向把**子连杆 `i+1`** 的 Wrench 传到当前连杆 `i` 时，用 `T[i+1]`。

```python
# forward: parent -> current
H = adjoint(inv_transform(T[i]))
V_i = H @ V_prev + robot.A[i] * qd[i]

# backward: child -> current
if i < n - 1:
    H_child = adjoint(inv_transform(T[i + 1]))
    F_i += H_child.T @ F_child

tau[i] = robot.A[i] @ F_i
```

这个索引错误是**参考系/连杆编号错误**，与转置本身是两件独立的事。

---

## 7. 空间惯性矩阵：将速度映射为动量的对偶量

刚体的动能可以写为

$$
K=\frac12 V_b^TG_bV_b.
$$

当 Body 坐标系原点位于质心时：

$$
G_b=
\begin{bmatrix}
I_C&0\\0&mI_3
\end{bmatrix}.
$$

六维动量定义为

$$
H_b=G_bV_b.
$$

因此可以把空间惯性矩阵 $G_b$ 看作一种将 Twist 映射到其对偶动量的线性算子；动能是二次配对 $\frac12V_b^TH_b$。

若 $V_s=AV_b$，为了保持动能不变：

$$
V_s^TG_sV_s=V_b^TG_bV_b.
$$

所以

$$
\boxed{G_s=A^{-T}G_bA^{-1}.}
$$

这与 Wrench 的逆转置变换一致，但惯性矩阵有**两个指标/两个向量位置**，因此左右各出现一个变换因子。

### 7.1 为什么不能直接写 $F=G\dot V$？

在旋转 Body 系中，刚体的 Newton–Euler 方程为

$$
\boxed{
F_b=G_b\dot V_b-\operatorname{ad}_{V_b}^TG_bV_b.
}
$$

额外项源于运动坐标系中的动量变化及其耦合。

其中

$$
\operatorname{ad}_{V}=
\begin{bmatrix}
[\omega]_\times&0\\
[v]_\times&[\omega]_\times
\end{bmatrix}.
$$

必须区分：

- $\operatorname{Ad}_T$：**有限坐标变换**，是 $6\times6$ 矩阵；
- $\operatorname{ad}_V$：**Lie bracket 的线性算子**，与速度相关；
- $\operatorname{ad}_V^T$：出现在动量/Wrench 的对偶动力学表达中。

不能因为名称相似就互换。

---

## 8. Twist 的导数不等于某点的物理加速度

对于 Body 原点（假设它固定在刚体上）：

$$
\boldsymbol v_b=R^T\dot p.
$$

对其求导：

$$
\dot{\boldsymbol v}_b
=R^T\ddot p-\boldsymbol\omega_b\times\boldsymbol v_b.
$$

因此

$$
\boxed{
R^T\ddot p
=\dot{\boldsymbol v}_b+\boldsymbol\omega_b\times\boldsymbol v_b.
}
$$

即使 $\dot{\boldsymbol v}_b=0$，物理加速度 $\ddot p$ 也可能不为零。

### 8.1 为什么第 1 根连杆匀速转动时 $\dot V_1=0$？

固定基座、单个固定螺旋轴：

$$
V_1=A_1\dot q_1.
$$

完整加速度递推为

$$
\dot V_1=A_1\ddot q_1+
\operatorname{ad}_{V_1}(A_1\dot q_1).
$$

由于 $V_1=A_1\dot q_1$，而 Lie bracket 具有反对称性，

$$
\operatorname{ad}_{V_1}V_1=[V_1,V_1]=0.
$$

故

$$
\boxed{\dot V_1=A_1\ddot q_1.}
$$

当 $\ddot q_1=0$，有 $\dot V_1=0$，但连杆上偏离旋转轴的点仍有向心加速度。

第二根连杆则不同：

$$
V_2=\operatorname{Ad}_{T_{1,2}^{-1}}V_1+A_2\dot q_2.
$$

因为 $V_2$ 通常不与 $A_2\dot q_2$ 平行，故

$$
\operatorname{ad}_{V_2}(A_2\dot q_2)
$$

不一定为零，即使 $\ddot q_1=\ddot q_2=0$。

---

## 9. 物理量与变换规则速查表

| 物理量 | 常见符号 | 类型/作用 | 典型变换 |
|---|---|---|---|
| 位姿 | $T\in SE(3)$ | 刚体坐标变换 | 矩阵复合 |
| Twist | $V=[\omega;v]$ | 速度向量 | $V_s=\operatorname{Ad}_T V_b$ |
| Wrench | $F=[m;f]$ | 速度的对偶量 | $F_s=\operatorname{Ad}_T^{-T}F_b$ |
| 关节速度 | $\dot q$ | 广义速度 | $V=J\dot q$ |
| 关节广义力 | $\tau$ | 关节速度的对偶量 | $\tau=J^TF$ |
| 空间动量 | $H=GV$ | 与 Twist 配对 | $H_s=\operatorname{Ad}_T^{-T}H_b$ |
| 空间惯性 | $G$ | Twist $\to$ 动量 | $G_s=\operatorname{Ad}_T^{-T}G_b\operatorname{Ad}_T^{-1}$ |

**通用检查法**：凡是遇到速度、力、力矩的坐标变换，先写下功率 $F^TV$ 是否保持不变。这个方法比机械记忆转置位置更可靠。

---

## 10. NumPy 中的行/列向量问题

数学上

$$
\tau_i=A_i^TF_i,\quad A_i,F_i\in\mathbb R^{6\times1}.
$$

但 NumPy 的一维数组形状通常为 `(6,)`，此时

```python
A.shape       # (6,)
A.T.shape     # (6,)  一维数组转置不改变形状
tau = A @ F   # 标量内积
```

因此 `A @ F` 与数学上的 $A^TF$ 对应。

如果显式使用二维列向量，则需要：

```python
A_col = A.reshape(6, 1)
F_col = F.reshape(6, 1)
tau = (A_col.T @ F_col).item()
```

**数学对象的“转置”与 NumPy 数组的 `.T` 操作不能脱离数组维度来理解。**

---

## 11. 自检题（建议隔几天闭卷回答）

1. 为什么 Space Twist 的线性分量通常不等于 Body 原点在 Space 系中的实际速度？
2. 若 $V_s=AV_b$，如何**仅通过功率不变性**推导 $F_s=A^{-T}F_b$？
3. 为什么 $J^TF$ 是关节广义力，而不是 $JF$？
4. RNEA 前向使用 `T[i]`，反向传递子连杆载荷为什么使用 `T[i+1]`？
5. 为什么单关节匀速运动可以有 $\dot V_1=0$，但连杆上的点仍有向心加速度？
6. 为什么空间惯性矩阵的变换同时含有左侧逆转置和右侧逆矩阵？

---

## 12. 与现有学习笔记的衔接

- **Day 8**：Wrench、静力学、Jacobian 转置与功率。
- **Day 11**：空间惯性、Body Newton–Euler、旋转坐标系中的导数。
- **Day 12**：RNEA 的前向/反向递推及实现中的参考系索引。

本笔记是针对上述讨论的**概念补充**，不是替代 Day 8–12 的主线笔记。
