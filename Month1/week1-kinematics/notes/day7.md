# Day 7｜逆运动学（Inverse Kinematics）学习总结

> 教材：Lynch & Park, *Modern Robotics*, Chapter 6  
> 学习内容：解析逆运动学、基于 Body Jacobian 的数值逆运动学、SE(3) 对数映射、阻尼最小二乘（DLS）与奇异性。

## 1. 正运动学与逆运动学

正运动学：给定关节变量 \(q\)，计算末端位姿




$$

$$

T(q)\in SE(3).
\]

逆运动学：给定目标位姿 \(T_d\)，寻找满足

\[
T(q)=T_d
\]

的关节变量 \(q\)。逆运动学可能存在多个解、无解或奇异构型附近的数值困难。

## 2. 平面 2R 机械臂的解析逆运动学

连杆长度 \(L_1,L_2\)，末端位置为

\[
\begin{aligned}
x&=L_1\cos q_1+L_2\cos(q_1+q_2),\\
y&=L_1\sin q_1+L_2\sin(q_1+q_2).
\end{aligned}
\]

目标位置 \((x_d,y_d)\) 满足可达条件

\[
|L_1-L_2|\leq \sqrt{x_d^2+y_d^2}\leq L_1+L_2.
\]

由余弦定理：

\[
\cos q_2=\frac{x_d^2+y_d^2-L_1^2-L_2^2}{2L_1L_2},
\]

\[
q_2=\pm\arccos\left(\frac{x_d^2+y_d^2-L_1^2-L_2^2}{2L_1L_2}\right),
\]

\[
q_1=\operatorname{atan2}(y_d,x_d)-\operatorname{atan2}(L_2\sin q_2,L_1+L_2\cos q_2).
\]

正负分支通常对应肘上、肘下两个位置解。**位置可达不代表任意指定的完整末端位姿可达**：2R 机械臂只有两个自由度。

**课堂例子**：\(L_1=2,L_2=1\)，\((x_d,y_d)=(2,1)\)，有 \(\cos q_2=0\)。两组解为

\[
(q_1,q_2)=(0,\pi/2),\qquad
(q_1,q_2)=(2\arctan(1/2),-\pi/2).
\]

## 3. 位置逆运动学：Newton–Raphson

给定位置误差 \(e_p=p_d-p(q)\)，线性化

$$
p(q+\Delta q)\approx p(q)+J_p(q)\Delta q.
$$

求解 \(J_p\Delta q\approx e_p\)，得到

\[
$\boxed{q^{(k+1)}=q^{(k)}+J_p(q^{(k)})^\dagger[p_d-p(q^{(k)})].}$
\]

\(J_p^\dagger\) 为 Moore–Penrose 伪逆，最小化线性化残差；在残差最小的解中选取最小范数的关节增量。

**2R 例子**：\(L_1=2,L_2=1\)，\(q=(0,0)\)，有

\[
J_p=\begin{bmatrix}0&0\\3&1\end{bmatrix},\qquad
J_p^\dagger=\begin{bmatrix}0&0.3\\0&0.1\end{bmatrix}.
\]

对目标 \((2,1)\)，当前位置为 \((3,0)\)，误差 \((-1,1)\)，因此首步 \(\Delta q=(0.3,0.1)\)。该构型的 \(J_p\) 秩亏，无法在一阶近似下产生 x 方向位移。

## 4. 完整位姿误差：为什么使用 SE(3) 对数

不能直接把 \(T_d-T(q)\) 当作六维误差：旋转矩阵处于非线性流形 \(SO(3)\)，矩阵相减不是 Twist。

定义 Body Frame 中的相对位姿

\[
$\boxed{T_{\mathrm{err}}=T(q)^{-1}T_d.}$
\]

目标达到时 \(T_{\mathrm{err}}=I_4\)。通过矩阵对数转换成六维误差坐标：

\[
$\boxed{V_b^{\mathrm{err}}=\operatorname{Log}_6(T(q)^{-1}T_d).}$
\]

这里的 \(V_b^{\mathrm{err}}\) 是**有限位姿差的对数坐标**，不是机器人当前的真实瞬时速度。

### 4.1 SE(3) 对数的计算

令



$$T=\begin{bmatrix}R&p\\0&1\end{bmatrix},\qquad
\log T=\begin{bmatrix}[\phi]_\times&\rho\\0&0\end{bmatrix}.$

$$

对一般旋转角 \(0<\theta<\pi\)，


$$\theta=\arccos\left(\frac{\operatorname{tr}R-1}{2}\right),\qquad
[\phi]_\times=\frac{\theta}{2\sin\theta}(R-R^T).$

旋转和平移通过 \(SO(3)\) 左 Jacobian 耦合：

\[
p=J_l(\phi)\rho,\qquad \rho=J_l(\phi)^{-1}p,
\]

\[
J_l(\phi)^{-1}=I-\frac12[\phi]_\times+
\left(\frac1{\theta^2}-\frac{1}{2\theta\tan(\theta/2)}\right)[\phi]_\times^2.
\]

所以

\[
\operatorname{Log}_6(T)=\begin{bmatrix}\phi\\\rho\end{bmatrix}.
\]

当 \(R=I\)（纯平移）时，\(\phi=0,\rho=p\)。对于非零旋转，一般 \(\rho\neq p\)。

**练习例子**：绕 z 轴旋转 \(\pi/2\) 且 \(p=(1,0,0)^T\)，有

\[
\operatorname{Log}_6(T)=
\begin{bmatrix}0&0&\pi/2&\pi/4&-\pi/4&0\end{bmatrix}^T.
\]

原因是刚体旋转与平移相互耦合；这也与瞬时关系 \(\dot p=\omega_s\times p+v_s\) 相呼应，但要区分有限变换的对数坐标与瞬时 Twist。

**数值注意**：直接使用 \(\theta/(2\sin\theta)\) 不适合 \(\theta\to0\) 或恰好 \(\theta=\pi\) 的情形；工程实现需针对小角度采用级数展开、对 \(\pi\) 附近采用稳健的旋转对数分支。

## 5. Body Jacobian 数值逆运动学

Body Jacobian 满足

\[
V_b=J_b(q)\dot q.
\]

对小的关节变化，一阶近似为

\[
$T(q+\Delta q)\approx T(q)\exp([J_b(q)\Delta q]).$
\]

与目标表达式 \(T_d=T(q)\exp([V_b^{\mathrm{err}}])\) 对比，得到局部近似

\[
$J_b(q)\Delta q\approx V_b^{\mathrm{err}}.$
\]

因此基本迭代为

\[
$\boxed{q^{(k+1)}=q^{(k)}+J_b(q^{(k)})^\dagger
\operatorname{Log}_6(T(q^{(k)})^{-1}T_d).}$
\]

这属于**局部线性化**，并不保证从任意初始值全局收敛。

停止条件可分别检查

\[
\|\phi_{\mathrm{err}}\|<\varepsilon_\omega,\qquad
\|\rho_{\mathrm{err}}\|<\varepsilon_v.
\]

旋转与平移的物理单位不同，分别设置阈值更合适。

### 5.1 伪逆公式回顾

- 列满秩：\(A^\dagger=(A^TA)^{-1}A^T\)。
- 行满秩：\(A^\dagger=A^T(AA^T)^{-1}\)。
- 一般情况：若 \(A=U\Sigma V^T\)，则 \(A^\dagger=V\Sigma^\dagger U^T\)。
- NumPy：`np.linalg.pinv(A)`。

### 5.2 Space 与 Body 的区别

\[
V_s^{\mathrm{err}}=\operatorname{Ad}_{T(q)}V_b^{\mathrm{err}},\qquad
J_s=\operatorname{Ad}_{T(q)}J_b.
\]

若用欧氏最小二乘伪逆，在非方阵或欠驱动情况下，Body 与 Space 形式的数值关节增量**不必相同**，因为 Adjoint 通常不是正交变换。

## 6. 奇异性与阻尼最小二乘（DLS）

在奇异点附近，某些奇异值 \(\sigma_i\) 非常小，伪逆对应增益 \(1/\sigma_i\) 可能极大，导致关节更新过大、线性化失效。

DLS 将关节增量的范数作为正则化惩罚：

\[
\boxed{\min_{\Delta q}\ \|J\Delta q-e\|^2+\lambda^2\|\Delta q\|^2,\qquad \lambda>0.}
\]

对 \(\Delta q\) 求导并令梯度为零：

\[
(J^TJ+\lambda^2I)\Delta q=J^Te.
\]

故

\[
\boxed{\Delta q=(J^TJ+\lambda^2I)^{-1}J^Te.}
\]

实现时优先使用

```python
A = J.T @ J + damping**2 * np.eye(J.shape[1])
b = J.T @ error
dq = np.linalg.solve(A, b)
```

而不是显式计算逆矩阵。

### 6.1 奇异值视角

对于 \(J=U\Sigma V^T\)，DLS 在第 i 个奇异方向的增益为

\[
\boxed{g_i(\lambda)=\frac{\sigma_i}{\sigma_i^2+\lambda^2}.}
\]

相比伪逆的 \(1/\sigma_i\)，DLS 对小奇异值方向进行抑制；\(\lambda\) 越大，通常关节增量越小，但线性化任务误差可能更大。

**课堂数值例子**：

\[
J=\operatorname{diag}(1,0.001),\quad e=(0,1)^T.
\]

普通伪逆给出 \(\Delta q=(0,1000)^T\)；取 \(\lambda=0.1\)，DLS 给出

\[
\Delta q=\left(0,\frac{0.001}{0.001^2+0.1^2}\right)^T\approx(0,0.09999)^T.
\]

### 6.2 与 Body IK 结合

将原先的伪逆步长替换为

\[
\boxed{\Delta q_k=(J_b(q_k)^TJ_b(q_k)+\lambda^2I)^{-1}J_b(q_k)^TV_b^{\mathrm{err}}.}
\]

固定阻尼实现简洁，但未处理关节限位、步长约束和全局收敛性。后续可研究自适应阻尼、线搜索或信赖域方法。

## 7. Python 练习与实测结果

本日已完成两个练习文件中的核心 TODO：

1. `day7_inverse_kinematics_exercises.py`：实现 `se3_log()` 和 `ik_body()`；单位测试包含纯平移、纯旋转、耦合变换和 100 组随机位姿。
2. `day7_damped_least_squares_exercises.py`：实现 `dls_step()` 和 `ik_body_dls()`；测试近奇异矩阵及可达目标的 Body IK。

**普通 Body IK**（目标 \(q=(0.5,0.7)\)）：约 3 次迭代后得到

\[
q\approx(0.49999814,\ 0.70000672).
\]

**DLS Body IK**（相同目标）：约 3 次迭代后得到

\[
q\approx(0.50004774,\ 0.69988569),
\]

角度误差范数约 \(6.66\times10^{-5}\)，平移对数误差范数约 \(6.18\times10^{-5}\)。

**近奇异平面 2R 位置 Jacobian**：奇异值约为

\[
\sigma=(3.16227728,\ 0.0006324555).
\]

| 方法                   | \(\Delta q_1\) | \(\Delta q_2\) |
| -------------------- | --------------:| --------------:|
| 伪逆                   | 5.0000         | -15.0000       |
| DLS \(\lambda=0.01\) | 0.019919       | -0.059761      |
| DLS \(\lambda=0.1\)  | 0.00019879     | -0.00060038    |
| DLS \(\lambda=0.5\)  | 0.000006829    | -0.000024390   |

**重要区别**：测试中的 `achieved_linear_dp = J @ dq` 是线性化预测位移，不是执行有限关节更新后由正运动学计算的真实位移。尤其对 \(\|\Delta q\|\) 很大的伪逆步长，二者可能相差显著。

另外，平面 2R 的**位置 Jacobian** 奇异，不等于六维 **Twist Jacobian** 必然奇异；奇异性必须结合所选任务和 Jacobian 理解。

## 8. 本日核心结论与后续衔接

- 解析 IK 可以揭示多解性与可达条件，但一般机器人常需数值 IK。
- 完整位姿误差适合通过 \(\operatorname{Log}_6(T^{-1}T_d)\) 表示；平移对数分量一般不等于简单的位置差。
- Body IK 使用 \(J_b^\dagger V_b^{\mathrm{err}}\) 进行局部 Newton 更新。
- DLS 是带 \(\ell_2\) 正则化的最小二乘问题，通过抑制小奇异值方向的增益改善数值稳定性。
- 阻尼与跟踪精度之间存在权衡；固定阻尼不是所有场景的最优选择。

**下一主题（Day 8）**：机器人静力学、Wrench、虚功原理，以及 \(\tau=J^TF\) 的 Jacobian 对偶关系。
