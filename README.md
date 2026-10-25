# Modern Robotics Learning

本仓库（`Mordern_Robotics_Learning`）用于记录个人在**机器人运动学、动力学、控制与机器人学习**方向的学习过程，保存课程笔记、公式推导、Python 实现及 MuJoCo 仿真实验。

仓库主要是**学习与实验记录**，而非成熟的机器人控制软件库。代码以理解原理、复现结果和便于修改参数为目标，部分内容仍在持续完善。

## 学习路线

整体学习沿着以下技术路线推进：

**机器人运动学与动力学 → 经典控制与仿真 → 连续控制强化学习 → 足式机器人控制与学习 → 控制理论与学习方法融合**

## Month 1：机器人基础与控制

### 参考教材与工具

**主要教材**

- Kevin M. Lynch and Frank C. Park, *Modern Robotics: Mechanics, Planning, and Control*（2017）。重点涉及第 2–6、8、9、11 章的相关内容。
- Russ Tedrake, *Underactuated Robotics*（[在线教材](https://underactuated.mit.edu/)）。作为后续欠驱动系统、浮动基座与接触动力学的补充参考。

**编程与仿真**

- Python、NumPy：数学计算、运动学与动力学算法实现。
- Matplotlib：轨迹、误差、力矩及系统响应可视化。
- MuJoCo：机器人动力学与控制仿真。

### 已学习内容

| 模块 | 主要内容 |
| --- | --- |
| 刚体运动与坐标变换 | SO(3)、SE(3)、旋转矩阵、齐次变换、指数映射、twist、screw axis、Adjoint、wrench |
| 正逆运动学 | 指数积公式（PoE）、Space/Body Jacobian、奇异性、可操作度、数值逆运动学 |
| 机器人动力学 | Euler–Lagrange 方程、质量矩阵、科里奥利与重力项、Newton–Euler 递推、正/逆动力学、能量分析 |
| 轨迹生成与经典控制 | 五次多项式轨迹、PD、重力补偿、Computed Torque Control（CTC）、模型失配、稳定性及输入饱和 |
| MuJoCo 仿真 | 2R 机械臂建模、状态与动力学量读取、控制器实现、轨迹跟踪及 Viewer 可视化 |
| 力与阻抗控制（Day 26） | 力控制、混合运动/力控制、关节空间与任务空间阻抗、期望惯量整形、无源性与能量分析 |

### 主要实践

- **2R 机械臂运动学与动力学**：实现并验证正运动学、Jacobian、逆运动学和动力学计算。
- **控制器对比**：实现 PD、重力补偿 PD 和 CTC，分析跟踪误差、模型不确定性和执行器饱和。
- **MuJoCo 2R 仿真**：使用仿真动力学量实现控制与可视化，并验证不同目标位置下的跟踪表现。
- **阻抗控制实验（Day 26）**：比较重力补偿 PD 与完整关节空间阻抗控制；进一步实现任务空间阻抗，验证刚度、阻尼和期望惯量对系统响应的影响。

### 后续学习

Month 1 后续将初步学习**浮动基座动力学、欠驱动特性和接触约束**，为之后的腿足机器人控制、全身控制（WBC）以及强化学习打下基础。

## 仓库说明

- 笔记与代码主要用于个人学习、复习和实验复现。
- 不同阶段的示例可能使用不同的模型参数与控制器约定，请以对应练习文件中的说明为准。
- 内容随学习进度持续更新，尚未涵盖完整的机器人学习与足式控制技术栈。


## 环境配置与依赖安装

本仓库主要使用 **Python** 实现机器人运动学、动力学、控制算法及仿真实验，推荐使用 Anaconda 管理 Python 环境。

### 1. 创建 Python 环境

建议使用 Python 3.11：

```bash
conda create -n robotics python=3.11
conda activate robotics
```

### 2. 安装依赖

```bash
pip install numpy scipy matplotlib sympy modern_robotics mujoco pytest
```

主要依赖及其用途如下：

| Python 包 | 用途 |
|---|---|
| `numpy` | 矩阵运算、线性代数、运动学与动力学计算 |
| `scipy` | 数值积分、优化、轨迹计算 |
| `matplotlib` | 数据可视化、控制响应曲线及动画 |
| `sympy` | 符号运算与动力学公式推导 |
| `modern_robotics` | *Modern Robotics* 教材配套算法库 |
| `mujoco` | 机器人动力学与控制仿真 |
| `pytest` | 算法实现与数值计算的单元测试 |

### 3. 验证安装

```bash
python -c "import numpy, scipy, matplotlib, sympy, modern_robotics, mujoco, pytest; print('All dependencies installed successfully!')"
```

部分实验使用 MuJoCo Viewer 进行可视化，需要支持 OpenGL 的图形环境。

**说明：** 本仓库中的部分算法（如正逆运动学、递归牛顿—欧拉动力学等）采用 NumPy 独立实现，以加深对算法原理的理解，并不完全依赖第三方机器人算法库。


## References

1. Lynch, K. M., & Park, F. C. (2017). *Modern Robotics: Mechanics, Planning, and Control*. Cambridge University Press. [Modern Robotics 官网](https://modernrobotics.northwestern.edu/)
2. Tedrake, R. *Underactuated Robotics*. [在线教材](https://underactuated.mit.edu/)
3. [MuJoCo Documentation](https://mujoco.readthedocs.io/)


