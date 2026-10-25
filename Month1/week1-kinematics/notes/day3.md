# ROBO_Learn · Day 3 学习笔记

**主题：SO(3) 的指数映射与对数映射（Exponential and Logarithm Maps）**  
**日期：2026-10-11**  
**参考：** Lynch & Park, *Modern Robotics*, 第 3 章（旋转的指数坐标与相关表示）

## 1. 本日学习目标与成果

本日的主线是建立两种三维旋转描述之间的联系：

- **旋转矩阵**：$R\in SO(3)$，满足 $R^TR=I$、$\det R=1$。
- **轴角表示**：绕单位轴 $\hat\omega$ 旋转角度 $\theta$。
- **指数坐标（旋转向量）**：$\phi=\theta\hat\omega\in\mathbb R^3$。

核心对应关系为

$$
\boxed{R=\exp([\phi]_\times)},\qquad
\boxed{\phi=\operatorname{Log}(R)^\vee}\quad\text{（选取适当的对数分支）}.
$$

**完成情况：** 已理解矩阵指数、Rodrigues 公式、对数映射及特殊情况；独立实现 `so3_exp` 与 `so3_log`，并完成数值测试。

## 2. 从角速度到矩阵指数

Day 2 已得到旋转矩阵的运动学关系：

$$
\dot R=[\omega_s]_\times R=R[\omega_b]_\times,
$$

其中 $\omega_s$ 是世界坐标系表达的角速度，$\omega_b$ 是本体坐标系表达的角速度。

当 **$\omega_s$ 为常向量** 时：

$$
R(t)=\exp([\omega_s]_\times t)R(0).
$$

当 **$\omega_b$ 为常向量** 时：

$$
R(t)=R(0)\exp([\omega_b]_\times t).
$$

因此，世界坐标系中的旋转增量**左乘**，本体坐标系中的旋转增量**右乘**。这与 Day 2 的旋转复合规则一致。

矩阵指数的定义是

$$
e^A=I+A+\frac{A^2}{2!}+\frac{A^3}{3!}+\cdots.
$$

当 $A$ 为实反对称矩阵时，$e^A\in SO(3)$。

## 3. 轴角与指数坐标

若 $\|\hat\omega\|=1$，定义旋转向量

$$
\phi=\theta\hat\omega,\qquad \theta=\|\phi\|.
$$

旋转矩阵可写为

$$
R=\exp([\hat\omega]_\times\theta)=\exp([\phi]_\times).
$$

对任意 $a=(a_x,a_y,a_z)^T$，其叉乘矩阵为

$$
[a]_\times=\begin{bmatrix}
0&-a_z&a_y\\
a_z&0&-a_x\\
-a_y&a_x&0
\end{bmatrix},\qquad [a]_\times b=a\times b.
$$

**注意：** 旋转向量不是唯一的。例如，$(0,0,\pi/2)^T$ 与 $(0,0,5\pi/2)^T$ 产生相同的旋转矩阵。指数映射不是一一对应的。

## 4. Rodrigues 公式及推导

记 $K=[\hat\omega]_\times$。对于单位轴，利用向量三重积：

$$
K^2=\hat\omega\hat\omega^T-I,\qquad K^3=-K.
$$

在指数级数中，奇数次幂归并为 $\sin\theta$，偶数次幂归并为 $1-\cos\theta$，得到

$$
\boxed{R=I+\sin\theta K+(1-\cos\theta)K^2.}
$$

**示例：** 绕 z 轴旋转 $\theta$：

$$
K=\begin{bmatrix}0&-1&0\\1&0&0\\0&0&0\end{bmatrix},\qquad
K^2=\operatorname{diag}(-1,-1,0),
$$

因此

$$
R_z(\theta)=\begin{bmatrix}
\cos\theta&-\sin\theta&0\\
\sin\theta&\cos\theta&0\\
0&0&1
\end{bmatrix}.
$$

**编程注意：** 角度使用弧度；对于 $\phi=0$，直接返回 $I$。任意角度均可代入 Rodrigues 公式，不必预先将角度归一化到 $[0,\pi]$。

## 5. 从旋转矩阵恢复轴角：SO(3) 对数映射

由 Rodrigues 公式取迹，得到

$$
\operatorname{tr}(R)=1+2\cos\theta,
$$

所以主值角度为

$$
\boxed{\theta=\arccos\left(\frac{\operatorname{tr}(R)-1}{2}\right)\in[0,\pi].}
$$

又因为 $K^T=-K$、$(K^2)^T=K^2$，有

$$
R-R^T=2\sin\theta K.
$$

因此，当 $0<\theta<\pi$ 时，

$$
\boxed{\hat\omega=\operatorname{vee}\left(\frac{R-R^T}{2\sin\theta}\right)},
$$
$$
\boxed{\phi=\frac{\theta}{2\sin\theta}\operatorname{vee}(R-R^T).}
$$

其中 `vee` 是 `skew` 的逆操作：

$$
\operatorname{vee}\!\left(\begin{bmatrix}0&-c&b\\c&0&-a\\-b&a&0\end{bmatrix}\right)
=\begin{bmatrix}a\\b\\c\end{bmatrix}.
$$

**示例：**

$$
R=\begin{bmatrix}0&-1&0\\1&0&0\\0&0&1\end{bmatrix}
\implies \operatorname{tr}(R)=1\implies \theta=\pi/2,
$$

由 $R-R^T$ 得 $\hat\omega=(0,0,1)^T$，故 $\phi=(0,0,\pi/2)^T$。

## 6. 两个特殊情况

### 6.1 零旋转：$\theta=0$

此时 $R=I$，轴方向不唯一；主值旋转向量约定为

$$
\boxed{\phi=0.}
$$

一般公式包含 $\sin\theta$ 分母，不能直接使用。练习代码以阈值 `theta < 1e-8` 处理。

### 6.2 半周旋转：$\theta=\pi$

此时 $\sin\theta=0$，一般公式再次失效。Rodrigues 公式简化为

$$
R=I+2K^2=2\hat\omega\hat\omega^T-I,
$$

即

$$
\boxed{R+I=2\hat\omega\hat\omega^T.}
$$

**几何解释：** 沿旋转轴方向的向量在主动旋转下保持不变，所以

$$
\boxed{R\hat\omega=\hat\omega},\quad\text{即}\quad (R-I)\hat\omega=0.
$$

因此，旋转轴是旋转矩阵**特征值为 1 的特征向量**。这一性质对任意非零角度的三维旋转均成立，并非只适用于 $\pi$。

例如

$$
R=\operatorname{diag}(-1,-1,1)
$$

对应绕 z 轴旋转 $\pi$，可取 $\hat\omega=(0,0,1)^T$ 或 $(0,0,-1)^T$。两种轴方向在 $\theta=\pi$ 时对应同一个旋转矩阵。

数值上可以利用 `np.linalg.eigh((R + R.T)/2)` 找到特征值最接近 1 的特征向量，并归一化。注意 `eig_vals == 1` 的精确比较不稳健。

## 7. NumPy 实现（本日练习版）

```python
import numpy as np


def skew(w):
    wx, wy, wz = np.asarray(w, dtype=float).reshape(3)
    return np.array([
        [0.0, -wz, wy],
        [wz, 0.0, -wx],
        [-wy, wx, 0.0]
    ])


def vee(W):
    return np.array([W[2, 1], W[0, 2], W[1, 0]])


def so3_exp(phi):
    phi = np.asarray(phi, dtype=float).reshape(3)
    theta = np.linalg.norm(phi)
    if theta == 0:
        return np.eye(3)
    omega = phi / theta
    K = skew(omega)
    return np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * (K @ K)


def so3_log(R):
    """教学用实现：输入假设为有效旋转矩阵。"""
    cos_theta = np.clip((np.trace(R) - 1) / 2, -1.0, 1.0)
    theta = np.arccos(cos_theta)

    if theta < 1e-8:
        return np.zeros(3)

    if np.isclose(theta, np.pi, atol=1e-6):
        R_sym = (R + R.T) / 2
        eig_vals, eig_vecs = np.linalg.eigh(R_sym)
        idx = np.argmin(np.abs(eig_vals - 1))
        omega = eig_vecs[:, idx]
        omega = omega / np.linalg.norm(omega)
        return theta * omega

    return theta / (2 * np.sin(theta)) * vee(R - R.T)
```

### 测试与判断标准

```python
# Exp 基本测试
assert np.allclose(so3_exp([0, 0, np.pi / 2]),
                   [[0, -1, 0], [1, 0, 0], [0, 0, 1]])
assert np.allclose(so3_exp([0, 0, 0]), np.eye(3))

phi = np.array([0.3, -0.5, 0.8])
R = so3_exp(phi)
assert np.allclose(R.T @ R, np.eye(3))
assert np.isclose(np.linalg.det(R), 1)
assert np.allclose(so3_log(R), phi)

# 特殊角度：比较旋转矩阵，而非直接比较旋转向量
for phi in [
    np.zeros(3),
    np.array([0.3, -0.5, 0.8]),
    np.array([0, 0, np.pi]),
    np.array([np.pi / np.sqrt(2), np.pi / np.sqrt(2), 0]),
]:
    R = so3_exp(phi)
    R_recovered = so3_exp(so3_log(R))
    assert np.allclose(R, R_recovered)
```

**实际练习结果：**

- 绕 z 轴旋转 $90^\circ$：通过。
- 任意旋转矩阵正交性：通过。
- 行列式：约为 `1.0000000000000002`，属浮点误差。
- 零旋转：通过。
- 一般情况 $\phi=(0.3,-0.5,0.8)^T$ 的往返误差：`5.551115123125783e-17`。
- 加入 $\theta=0$、$\theta=\pi$ 特殊分支后，四组重构测试通过。

**关键判断标准：** 一般不能要求 `so3_log(so3_exp(phi)) == phi` 对任意 $\phi$ 成立；更可靠的检查是

$$
\boxed{\operatorname{Exp}(\operatorname{Log}(R))=R.}
$$

## 8. 数值实现的边界与后续改进

上面的代码是**本日学习用版本**，不是面向所有病态输入的工业级实现：

1. `arccos` 输入要裁剪到 $[-1,1]$，避免舍入误差。
2. 角度接近 0 时，直接返回零向量会忽略很小但非零的旋转；更稳健的实现需要小角度展开。
3. 角度接近 $\pi$ 时，迹公式与反对称部分可能数值敏感；特征向量法需要结合输入正交性、阈值与误差控制。近似 $\pi$ 时，轴方向符号还应与 $R-R^T$ 的信息保持一致，避免重构角度符号错误。
4. 输入若不是有效的 $SO(3)$ 矩阵，需要另行验证或投影。

这些问题可在后续机器人姿态计算时深入处理，本日以理解数学结构和完成基础实现为主。

## 9. 与其他旋转表示的关系

机器人学中还会遇到欧拉角、Roll–Pitch–Yaw (RPY)、单位四元数等表示。本阶段优先掌握旋转矩阵、轴角与指数坐标；后续学习姿态控制、状态估计和机器人仿真时，再系统补充四元数。

指数坐标的意义还在于可以自然推广到刚体运动：

$$
SO(3):\ R=e^{[\hat\omega]_\times\theta}
\quad\longrightarrow\quad
SE(3):\ T=e^{[\mathcal S]\theta}.
$$

这将连接后续的刚体位姿、twist、正运动学与 PoE（Product of Exponentials）公式。

## 10. 本日回顾问题

1. 为什么恒定世界系角速度对应指数矩阵**左乘**，而恒定本体系角速度对应**右乘**？
2. 如何从 $K^3=-K$ 推出 Rodrigues 公式？
3. 为什么 $\phi=(0,0,\pi/2)^T$ 与 $(0,0,5\pi/2)^T$ 对应相同的 $R$？
4. 为什么 $\operatorname{tr}(R)$ 可以用于恢复主值旋转角度？
5. 为什么 $\theta=0$ 与 $\theta=\pi$ 要特殊处理？
6. 为什么旋转轴满足 $R\hat\omega=\hat\omega$？
7. 为什么验证 `Exp(Log(R)) == R` 比总是比较恢复前后的旋转向量更合适？

---

**Day 3 状态：完成。** 下一步：复习 SO(3) 的关键公式，并进入 $SE(3)$ 的刚体位姿表示与齐次变换。
