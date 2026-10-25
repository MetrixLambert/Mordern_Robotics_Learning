# ROBO_Learn — Day 5 学习笔记

## Product of Exponentials (PoE) 与机器人正运动学

**参考教材**：Kevin M. Lynch & Frank C. Park, *Modern Robotics: Mechanics, Planning, and Control*，第 4 章（重点 §4.1）。

**学习成果**：理解旋转/移动关节的 Screw Axis，掌握 Space PoE 和 Body PoE 的定义及等价关系；独立实现正运动学和 Screw Axis 坐标转换，并通过固定案例与随机一致性测试。

---

## 1. 正运动学问题

对于串联机器人，给定关节变量 $q=(q_1,\ldots,q_n)$，求末端坐标系 $\{b\}$ 相对于固定空间坐标系 $\{s\}$ 的位姿：

$$
T_{sb}(q)=\begin{bmatrix}R_{sb}(q)&p_{sb}(q)\\0&1\end{bmatrix}\in SE(3).
$$

PoE 的核心思路是：**用每个关节对应的刚体运动指数映射表示整个机械臂的正运动学。**

约定六维向量采用**角速度在前、线速度在后**的排列：

$$
\mathcal S=\begin{bmatrix}\omega\\v\end{bmatrix},\qquad
[\mathcal S]=\begin{bmatrix}[\omega]_\times&v\\0&0\end{bmatrix}.
$$

## 2. Screw Axis：关节运动的六维表示

### 2.1 旋转关节（Revolute）

设旋转轴的单位方向为 $\omega$，轴线上任意一点为 $a$，两者均用同一个坐标系表示，则

$$
\boxed{\mathcal S=\begin{bmatrix}\omega\\-\omega\times a\end{bmatrix}},\qquad \|\omega\|=1.
$$

关节转角为 $\theta$ 时，产生的刚体变换为

$$
e^{[\mathcal S]\theta}
=\begin{bmatrix}R&(I-R)a\\0&1\end{bmatrix},
\qquad R=e^{[\omega]_\times\theta}.
$$

这里 $a$ 是**旋转轴上的一点**，不一定是关节处某个专门建立的坐标系原点。$v=-\omega\times a$ 是 Screw Axis 的线速度部分，不应直接误认为末端原点的实际速度。

### 2.2 移动关节（Prismatic）

设移动方向为单位向量 $u$，则

$$
\boxed{\mathcal S=\begin{bmatrix}0\\u\end{bmatrix}},\qquad \|u\|=1,
$$

关节位移为 $d$ 时：

$$
e^{[\mathcal S]d}=\begin{bmatrix}I&du\\0&1\end{bmatrix}.
$$

旋转关节的变量通常以弧度计，移动关节的变量通常以长度计。归一化使关节变量具有直接的物理意义。

## 3. Home Configuration 与 Space PoE

令

$$
M=T_{sb}(0)
$$

表示**所有关节处于零位时，所选末端坐标系相对于空间坐标系的位姿**。

在零位构型下，将每个关节轴都表示在同一个固定空间坐标系中，得到 $\mathcal S_1,\ldots,\mathcal S_n$。则

$$
\boxed{
T_{sb}(q)=e^{[\mathcal S_1]q_1}e^{[\mathcal S_2]q_2}\cdots e^{[\mathcal S_n]q_n}M
}.
$$

**理解要点**：

- 所有 $\mathcal S_i$ 都在**零位构型**下确定，计算过程中保持不变。
- 机器人运动后，下游关节的**实际空间轴线**会随上游关节运动；这已由指数矩阵的有序乘积自动体现，不需要在 Space PoE 中逐次重写 $\mathcal S_i$。
- 不必为每个关节单独建立坐标系才能使用 PoE。
- 矩阵乘法通常不交换，指数项的顺序不能任意调换。
- 物理上以不同先后顺序到达同一组关节变量，并不意味着 PoE 公式中的矩阵可以交换顺序。

## 4. 平面 2R 机械臂

两段连杆长度为 $L_1=2,L_2=1$，两个关节均绕 $z$ 轴转动。零位时机械臂沿 $+x$ 方向伸直，末端坐标系方向与空间坐标系一致。

$$
\mathcal S_1=\begin{bmatrix}0\\0\\1\\0\\0\\0\end{bmatrix},\qquad
\mathcal S_2=\begin{bmatrix}0\\0\\1\\0\\-2\\0\end{bmatrix},\qquad
M=\begin{bmatrix}I&(3,0,0)^T\\0&1\end{bmatrix}.
$$

由 Space PoE 可得到熟悉的几何正运动学：

$$
\begin{aligned}
p_x&=2\cos\theta_1+\cos(\theta_1+\theta_2),\\
p_y&=2\sin\theta_1+\sin(\theta_1+\theta_2),\\
p_z&=0,\\
R&=R_z(\theta_1+\theta_2).
\end{aligned}
$$

检验案例：

| 关节变量                              | 末端位置         |
| --------------------------------- | ------------ |
| $\theta_1=0,\ \theta_2=\pi/2$     | $(2,1,0)^T$  |
| $\theta_1=\pi/2,\ \theta_2=\pi/2$ | $(-1,2,0)^T$ |

**易错点**：第二段连杆相对空间系的绝对角度为 $\theta_1+\theta_2$，而不是单独的 $\theta_2$。

## 5. RP 机械臂：旋转关节 + 移动关节

第一关节绕空间 $z$ 轴旋转，第二关节在零位时沿空间 $x$ 方向平移；零位末端位于 $(2,0,0)^T$：

$$
\mathcal S_1=\begin{bmatrix}0\\0\\1\\0\\0\\0\end{bmatrix},\quad
\mathcal S_2=\begin{bmatrix}0\\0\\0\\1\\0\\0\end{bmatrix},\quad
M=\begin{bmatrix}I&(2,0,0)^T\\0&1\end{bmatrix}.
$$

因此

$$
T(\theta,d)=e^{[\mathcal S_1]\theta}e^{[\mathcal S_2]d}M
=\begin{bmatrix}R_z(\theta)&R_z(\theta)(2+d,0,0)^T\\0&1\end{bmatrix}.
$$

当 $\theta=\pi/2,d=1$ 时：

$$
\boxed{p=(0,3,0)^T},\qquad
\boxed{R=\begin{bmatrix}0&-1&0\\1&0&0\\0&0&1\end{bmatrix}}.
$$

**物理解释**：第二关节的平移方向会随着第一关节转动；零位时沿 $x$，此时沿空间 $y$。下游移动关节不改变上游旋转关节的轴线，但其伸长量会影响上游关节转动时末端的线速度大小，例如本例有 $\|\partial p/\partial\theta\,\dot\theta\|=|2+d|\,|\dot\theta|$。这为后续 Jacobian 学习埋下伏笔。

## 6. Body Screw Axis 与 Body PoE

### 6.1 定义

与 Space Screw Axis 不同，Body Screw Axis $\mathcal B_i$ 使用**末端在零位时的坐标系**来表达第 $i$ 个关节轴。

二者由 Day 4 学习的 Adjoint 联系：

$$
\boxed{\mathcal B_i=\operatorname{Ad}_{M^{-1}}\mathcal S_i},\qquad
\boxed{\mathcal S_i=\operatorname{Ad}_{M}\mathcal B_i}.
$$

注意这里用的是**零位变换 $M$**，而不是随关节变量变化的 $T_{sb}(q)$。$\mathcal B_i$ 本身也是固定的零位参数，不是关节轴在当前运动末端系下的实时坐标。

### 6.2 Body PoE 公式

$$
\boxed{
T_{sb}(q)=M e^{[\mathcal B_1]q_1}e^{[\mathcal B_2]q_2}\cdots e^{[\mathcal B_n]q_n}
}.
$$

与 Space PoE 的关键区别是 **$M$ 在左侧，而不是右侧**。

### 6.3 两种 PoE 为什么等价？

由于 $\mathcal S_i=\operatorname{Ad}_M\mathcal B_i$，有

$$
e^{[\mathcal S_i]q_i}=M e^{[\mathcal B_i]q_i}M^{-1}.
$$

于是

$$
\begin{aligned}
T_{sb}(q)
&=e^{[\mathcal S_1]q_1}\cdots e^{[\mathcal S_n]q_n}M\\
&=(Me^{[\mathcal B_1]q_1}M^{-1})\cdots(Me^{[\mathcal B_n]q_n}M^{-1})M\\
&=Me^{[\mathcal B_1]q_1}\cdots e^{[\mathcal B_n]q_n}.
\end{aligned}
$$

### 6.4 2R 机械臂的 Body Screw Axis

对于第 4 节的机械臂，零位末端原点在 $(3,0,0)$，两个旋转轴相对于该原点的位置分别为 $(-3,0,0)$、$(-1,0,0)$，因此

$$
\boxed{\mathcal B_1=(0,0,1,0,3,0)^T},\qquad
\boxed{\mathcal B_2=(0,0,1,0,1,0)^T}.
$$

## 7. 末端坐标系改变时，哪些参数会变？

假设 2R 机械臂的物理结构、空间坐标系不变，仅将所选末端坐标系的**零位原点**从 $(3,0,0)$ 移至 $(4,0,0)$，方向不变：

| 参数                          | 是否改变 | 新值或原因                                |
| --------------------------- | ---- | ------------------------------------ |
| $\mathcal S_1,\mathcal S_2$ | 不变   | 关节轴在固定空间坐标系下的表示未变                    |
| $\mathcal B_1$              | 改变   | $(0,0,1,0,4,0)^T$                    |
| $\mathcal B_2$              | 改变   | $(0,0,1,0,2,0)^T$                    |
| $M$                         | 改变   | $M_{\mathrm{new}}=[I,(4,0,0)^T;0,1]$ |

**关键认识**：$M$ 描述的是**所选末端坐标系**在零位时的位姿，不是独立于坐标系选择的“机器人实体末端位置”。所选坐标系原点甚至可以位于机械臂实体之外。

## 8. Python 实践与实现要点

本日练习复用了 Day 4 的 `so3_exp`、`se3_exp`、`adjoint`、`trans_inv` 等函数，并独立完成以下三个函数。

### 8.1 Space PoE

```python
def fk_space(M, Slist, q):
    T = M.copy()
    for i in range(len(q) - 1, -1, -1):
        T = se3_exp(Slist[:, i], q[i]) @ T
    return T
```

从 $M$ 开始，**倒序遍历、逐次左乘**，保证最终乘积为 $e^{[\mathcal S_1]q_1}\cdots e^{[\mathcal S_n]q_n}M$。

### 8.2 Space → Body Screw Axis

```python
def space_to_body_screws(M, Slist):
    M_inv = trans_inv(M)
    return adjoint(M_inv) @ Slist
```

约定 `Slist` 形状为 `(6, n)`，每一列是一个关节的 Screw Axis。

### 8.3 Body PoE

```python
def fk_body(M, Blist, q):
    T = M.copy()
    for i in range(len(q)):
        T = T @ se3_exp(Blist[:, i], q[i])
    return T
```

从 $M$ 开始，**正序遍历、逐次右乘**。

### 8.4 测试记录

以下测试均已通过：

- `test_fk_space_2r()`：平面 2R 的位姿计算。

- `test_fk_space_rp()`：RP 机械臂的旋转与平移耦合。

- `test_space_to_body()`：已知解析值的 Screw Axis 转换。

- `test_fk_body()`：同一 2R 机械臂的 Space/Body PoE 结果一致。

- `test_random_consistency()`：20 组随机末端零位位姿、旋转/移动混合关节轴与关节变量，验证
  
  $$
  \operatorname{FK}_{\mathrm{space}}(M,S,q)
\approx
\operatorname{FK}_{\mathrm{body}}(M,\operatorname{Ad}_{M^{-1}}S,q).
  $$

随机测试主要检验**两种实现的内部一致性**；与已知几何解对照的 2R、RP 测试则提供了独立的正确性验证。

## 9. 本日知识清单

- [x] 区分旋转关节与移动关节的 Screw Axis。
- [x] 由关节轴方向和轴上一点构造 $\mathcal S$。
- [x] 理解 $M$ 与零位构型的关系。
- [x] 理解 Space PoE 中关节轴固定表示与实际轴线运动的区别。
- [x] 独立计算平面 2R 与 RP 机械臂的正运动学。
- [x] 理解 Body Screw Axis 的零位坐标系含义。
- [x] 利用 Adjoint 转换 Space/Body Screw Axis。
- [x] 推导 Space PoE 与 Body PoE 的等价性。
- [x] 独立实现并测试 `fk_space`、`fk_body`、`space_to_body_screws`。

## 10. 后续衔接

Day 5 解决了**给定关节变量，末端在哪里**的问题。下一步自然进入 **Jacobian 与速度运动学**：给定关节速度 $\dot q$，求末端 Twist $V$，并理解为什么下游关节位置会影响上游关节对末端线速度的贡献。

---

*ROBO_Learn · Day 5 完成*
