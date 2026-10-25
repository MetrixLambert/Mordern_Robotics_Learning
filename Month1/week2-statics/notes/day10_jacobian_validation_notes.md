# Day 10 学习笔记：Jacobian 的数值验证与 Week 2 运动学总结

> 课程：Modern Robotics（Lynch & Park）｜主题：有限差分、Space/Body Jacobian 一致性与测试方法

## 1. 学习目标

1. 从正运动学独立地数值估计末端位置 Jacobian 与 Space Jacobian。
2. 理解有限差分（Finite Difference，FD）与解析 Jacobian 的关系。
3. 使用 Adjoint 关系验证 Space/Body Jacobian 的坐标一致性。
4. 区分“两个程序相互一致”与“程序实际正确”。

## 2. 位置 Jacobian：中心有限差分

对于末端位置 $p(q)\in\mathbb R^m$，位置 Jacobian 定义为

$$
J_p(q)=\frac{\partial p}{\partial q}\in\mathbb R^{m\times n}.
$$

对第 $i$ 个关节单独施加扰动，中心差分给出

$$
\boxed{J_{p,i}^{\mathrm{FD}}(q)=\frac{p(q+\varepsilon e_i)-p(q-\varepsilon e_i)}{2\varepsilon}.}
$$

其中 $e_i$ 是第 $i$ 个标准基向量，$\varepsilon>0$ 是差分步长。比较误差时使用 Frobenius 范数：

$$
e_{p,\mathrm{FD}}=\|J_p-J_p^{\mathrm{FD}}\|_F.
$$

**实现要点：** 用 `q.copy()` 创建 `q_plus` 和 `q_minus`，只扰动一个分量；Jacobian 的第 $i$ **列**对应第 $i$ 个关节。

## 3. SE(3) 上的 Space Jacobian 数值验证

正运动学写为 $T(q)\in SE(3)$。Space Twist 的定义是

$$
[V_s]=\dot T T^{-1},\qquad V_s=J_s(q)\dot q.
$$

因此，对第 $i$ 个关节进行扰动，定义

$$
T_+=T(q+\varepsilon e_i),\qquad T_-=T(q-\varepsilon e_i).
$$

利用 $SE(3)$ 的矩阵对数，可以构造

$$
\boxed{J_{s,i}^{\mathrm{FD}}\approx\frac{\operatorname{Log}_6(T_+T_-^{-1})}{2\varepsilon}.}
$$

这里 $\operatorname{Log}_6$ 表示先计算 $SE(3)$ 矩阵对数，再通过 vee 映射转成角速度在前的六维向量 $(\omega,v)$。乘法顺序 **$T_+T_-^{-1}$** 对应 Space 表达。

类似地，Body Jacobian 可通过

$$
\boxed{J_{b,i}^{\mathrm{FD}}\approx\frac{\operatorname{Log}_6(T_-^{-1}T_+)}{2\varepsilon}}
$$

数值估计。乘法顺序 **$T_-^{-1}T_+$** 对应 Body 表达。

在 PoE 模型中，固定其他关节、只改变一个关节时，上述相对位姿沿同一单参数子群变化；在矩阵对数的适当分支内，该差分公式可具有比一般光滑轨迹更强的精确性，但实际计算仍受浮点和对数实现误差影响。

解析与数值 Space Jacobian 的误差为

$$
e_{s,\mathrm{FD}}=\|J_s-J_s^{\mathrm{FD}}\|_F.
$$

## 4. Space/Body Jacobian 的 Adjoint 一致性

对于同一个末端位姿 $T(q)$，

$$
\boxed{J_s(q)=\operatorname{Ad}_{T(q)}J_b(q)}
$$

其中

$$
\operatorname{Ad}_T=
\begin{bmatrix}
R&0\\
[p]_\times R&R
\end{bmatrix},\qquad
T=\begin{bmatrix}R&p\\0&1\end{bmatrix}.
$$

对应检验误差为

$$
e_{\mathrm{adj}}=\|J_s-\operatorname{Ad}_T J_b\|_F.
$$

**重要认识：** 即使 $e_{\mathrm{adj}}\approx0$，也不能单独证明 $J_s$ 和 $J_b$ 都正确。两个实现可能共享同一错误。例如任取非零矩阵 $E$，令 $\widetilde J_b=J_b+E$、$\widetilde J_s=\operatorname{Ad}_T\widetilde J_b$，一致性检验仍然通过。因此应当组合解析特例、正运动学差分与 Adjoint 一致性测试。

## 5. 三项测试及实际结果

Day 10 的 `verify_jacobian_consistency()` 返回三个**标量误差**，而不是三个 Jacobian 矩阵。

| 字典键 | 数学表达式 | 用户本次运行结果 |
|---|---|---:|
| `space_fd` | $\|J_s-J_s^{\mathrm{FD}}\|_F$ | $4.424\times10^{-11}$ |
| `body_adjoint` | $\|J_s-\operatorname{Ad}_T J_b\|_F$ | $7.964\times10^{-16}$ |
| `position_fd` | $\|J_p-J_p^{\mathrm{FD}}\|_F$ | $4.854\times10^{-11}$ |

另外两个独立练习的测试结果为：位置 Jacobian FD 误差约 $9.037\times10^{-11}$，Space Jacobian FD 误差约 $4.424\times10^{-11}$。三个练习均通过测试，并生成 `day10_fd_errors.png`。

核心 Python 表达式：

```python
Jp_fd[:, i] = (position_fn(q_plus) - position_fn(q_minus)) / (2 * eps)
Js_fd[:, i] = se3_log_vector(T_plus @ trans_inv(T_minus)) / (2 * eps)
space_fd_error = np.linalg.norm(Js - Js_fd, ord="fro")
body_adjoint_error = np.linalg.norm(Js - adjoint(T) @ Jb, ord="fro")
position_fd_error = np.linalg.norm(Jp - Jp_fd, ord="fro")
```

## 6. 差分步长与数值误差

对于一般光滑函数的中心差分，误差可概括为

$$
\text{误差}\sim O(\varepsilon^2)+O\left(\frac{\epsilon_{\mathrm{mach}}}{\varepsilon}\right).
$$

- $O(\varepsilon^2)$：中心差分截断误差，通常随步长缩小而下降。
- $O(\epsilon_{\mathrm{mach}}/\varepsilon)$：浮点舍入误差及相近量相减带来的放大效应，可能在步长过小时上升。

所以 **$\varepsilon$ 不是越小越好**。实际应在多个步长（例如 $10^{-2}$ 至 $10^{-8}$）下比较误差；具体误差曲线依赖模型、实现和数值精度。对于单关节 PoE 的 Lie 群对数差分，上面的普通中心差分截断误差估计可能偏保守。

## 7. 与 Week 2 知识的联系

- **PoE 正运动学**：$T(q)=e^{[S_1]q_1}\cdots e^{[S_n]q_n}M$，为数值差分提供 $T(q\pm\varepsilon e_i)$。
- **Space/Body Jacobian**：分别对应 $[V_s]=\dot T T^{-1}$ 与 $[V_b]=T^{-1}\dot T$。
- **位置 Jacobian**：末端原点速度 $\dot p=J_p\dot q$，且 $J_p=[-[p]_\times\ \ I]J_s$（角速度在前的约定）。
- **逆运动学**：Jacobian 将关节增量映射为任务空间运动，误差会直接影响迭代更新。
- **静力学**：$\tau=J^TF$；Jacobian 的正确性也关系到关节力矩计算。
- **可操作度**：奇异值、速度椭球、条件数都依赖正确的 Jacobian。

## 8. 实践建议

1. 先用已知解析几何模型验证 FK，例如平面 2R 的末端位置。
2. 在多个随机、非奇异构型下运行有限差分测试。
3. 同时检查 $J_s=\operatorname{Ad}_T J_b$。
4. 对位置任务，额外检查 $J_p=[-[p]_\times\ \ I]J_s$ 与位置差分结果是否一致。
5. 为测试设置合理容差，避免把机器精度附近的波动误判为程序错误。

## 9. 关键结论

> 一致性检验说明不同计算路径遵循某个理论关系；独立数值验证进一步检查结果是否描述了真实运动。二者互补，而非互相替代。

**Day 10 完成，接下来进入 Week 2 综合复习。**
