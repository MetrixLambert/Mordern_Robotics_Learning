
"""
ROBO_Learn - Month 2 Day 19
Footstep Planning and Gait Timing

目标：
1. 生成左右脚交替落脚的二维足印
2. 生成单支撑步态时序
3. 预测未来 DCM
4. 实现 DCM 反馈落脚点修正

依赖：numpy, matplotlib

完成 TODO 1-4 后运行。
"""

from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt


@dataclass
class Config:
    g: float = 9.81
    zc: float = 1.0

    step_length: float = 0.30
    step_width: float = 0.20
    step_time: float = 0.50

    n_steps: int = 8

    disturbance: np.ndarray = None

    def __post_init__(self):
        if self.disturbance is None:
            self.disturbance = np.array([0.05, 0.0])

    @property
    def omega(self):
        return np.sqrt(self.g / self.zc)


def generate_footsteps(cfg):
    """
    TODO 1：生成左右脚交替落脚的二维足印。

    返回：
        footsteps: ndarray, shape (n_steps, 2)

    约定：
        x：前进方向
        y：横向方向

        第 0 步：左脚
        第 1 步：右脚
        第 2 步：左脚
        ...

    要求：
        每一步 x 增加 step_length
        左脚 y = +step_width/2
        右脚 y = -step_width/2

    注意：
        此处的 step_length 指相邻足印的前向距离。
        只生成几何足印，不处理 DCM。
    """
    footsteps = np.zeros((cfg.n_steps, 2))

    # TODO 1
    for i in range(cfg.n_steps): 
        x = i * cfg.step_length
        y = cfg.step_width / 2 if i % 2 == 0 else -cfg.step_width / 2
        footsteps[i] = np.array([x, y])
    return footsteps


def generate_gait_timing(cfg):
    """
    TODO 2：生成步态时间序列。

    返回：
        times: ndarray, shape (n_steps + 1,)

    第 k 步开始时间：
        t_k = k * step_time

    暂时假设：
        双支撑时间 = 0
        每个单支撑周期时长相同
    """
    # TODO 2
    times = np.arange(0, cfg.n_steps * cfg.step_time, cfg.step_time) 
    return times

def predict_dcm(xi0, p, T, omega):
    """
    TODO 3：预测 T 秒后的二维 DCM。

    输入：
        xi0: shape (2,)
        p:   shape (2,)
        T:   float
        omega: float

    返回：
        xi_T: shape (2,)

    提示：
        x、y 两个方向独立满足 LIPM 动力学。

    注意：
        此处假设整个预测区间 ZMP 固定。
    """
    xi_T = p + (xi0 - p) * np.exp(omega * T) 
    return xi_T


def correct_footstep(p_nom, xi_actual, xi_nom, T, omega):
    """
    TODO 4：DCM 反馈落脚点修正。

    输入：
        p_nom: 原计划落脚点
        xi_actual: 当前实际 DCM
        xi_nom: 当前期望 DCM
        T: 距离下一次触地的剩余时间

    返回：
        p_corrected: 修正后的落脚点

    提示：
        先计算当前 DCM 误差。
        然后预测触地时刻的误差。

    暂时不考虑：
        关节限位
        最大步长
        最大步宽
        摆动脚速度约束
    """
    # TODO 4
    p_corrected = p_nom + (xi_actual - xi_nom) * np.exp(omega * T) 
    return p_corrected



def periodic_dcm_offset(step_index, cfg):
    """
    TODO 5:
    计算第 step_index 步开始时的周期 DCM 偏移。

    返回：
        ndarray, shape (2,)

    提示：
        bx = L / (exp(omega*T) - 1)

        by 的符号随左右支撑脚交替改变。
    """
    bx = cfg.step_length / (np.exp(cfg.omega * cfg.step_time) - 1) 
    by = - cfg.step_width / (np.exp(cfg.omega * cfg.step_time) + 1) if step_index % 2 == 0 else cfg.step_width / (np.exp(cfg.omega * cfg.step_time) + 1)

    # dcm_offset = np.array([bx, by]) 
    dcm_offset = np.array([step_index * cfg.step_length, cfg.step_width / 2 if step_index % 2 == 0 else -cfg.step_width / 2])
    dcm_offset = dcm_offset + np.array([bx, by])
    return dcm_offset


def generate_dcm_trajectory(footsteps, cfg, dt=0.01):
    """
    TODO 6:
    根据足印序列生成连续二维 DCM 轨迹。

    输入：
        footsteps: (N, 2)
        cfg: Config
        dt: 采样时间间隔

    返回：
        times: (M,)
        xi_ref: (M, 2)

    约定：
        每个单支撑周期的 ZMP 固定在当前支撑脚。
        不包含双支撑阶段。
    """
    times = []
    xi_ref = []

    for k in range(len(footsteps)):
        p = footsteps[k]

        # TODO 6.1: 计算当前步态开始时的周期 DCM
        xi_start = periodic_dcm_offset(k, cfg)

        # 当前单支撑阶段内的采样时刻
        local_times = np.arange(0, cfg.step_time, dt)

        for tau in local_times:
            # TODO 6.2:
            # 使用 DCM 解析解计算 xi(t)
            xi_t = p + (xi_start - p) * np.exp(cfg.omega * tau) 

            times.append(k * cfg.step_time + tau)
            xi_ref.append(xi_t)

    return np.asarray(times), np.asarray(xi_ref)


def plot_footsteps(footsteps, corrected=None):
    fig, ax = plt.subplots(figsize=(10, 5))

    for i, p in enumerate(footsteps):
        color = "tab:blue" if i % 2 == 0 else "tab:orange"

        ax.scatter(p[0], p[1], c=color, s=100)
        ax.text(p[0], p[1] + 0.025, f"{i}", ha="center")

        if i > 0:
            ax.plot(
                footsteps[i-1:i+1, 0],
                footsteps[i-1:i+1, 1],
                "k--",
                alpha=0.3
            )

    if corrected is not None:
        ax.scatter(
            corrected[0], corrected[1],
            c="red", marker="x", s=150,
            label="Corrected Footstep"
        )
        ax.legend()

    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.set_title("Day 19 - Footstep Planning")
    ax.axis("equal")
    ax.grid(True)


def main():
    cfg = Config()

    footsteps = generate_footsteps(cfg)
    times = generate_gait_timing(cfg)

    print("Footsteps:")
    print(footsteps)

    print("\nGait timing:")
    print(times)

    # DCM 预测测试
    xi0 = np.array([0.1, 0.05])
    p0 = np.array([0.0, 0.1])

    xi_T = predict_dcm(
        xi0, p0, cfg.step_time, cfg.omega
    )

    print("\nPredicted DCM:")
    print(xi_T)

    # DCM 扰动反馈测试
    xi_nom = xi0.copy()
    xi_actual = xi_nom + cfg.disturbance

    p_nom = footsteps[1]

    p_corrected = correct_footstep(
        p_nom,
        xi_actual,
        xi_nom,
        cfg.step_time,
        cfg.omega
    )

    print("\nNominal footstep:", p_nom)
    print("Corrected footstep:", p_corrected)
    print("Correction:", p_corrected - p_nom)

    # 验证修正后的触地 DCM 相对足印的偏移不变
    xi_nom_T = predict_dcm(
        xi_nom, p0, cfg.step_time, cfg.omega
    )
    xi_actual_T = predict_dcm(
        xi_actual, p0, cfg.step_time, cfg.omega
    )

    np.testing.assert_allclose(
        xi_actual_T - p_corrected,
        xi_nom_T - p_nom,
        atol=1e-10
    )

    plot_footsteps(footsteps, p_corrected)

    times, xi_ref = generate_dcm_trajectory(footsteps, cfg)

    plt.figure(figsize=(10, 5))
    plt.plot(xi_ref[:, 0], xi_ref[:, 1], label="DCM")
    plt.scatter(footsteps[:, 0], footsteps[:, 1],
                c="red", label="Footsteps")
    plt.axis("equal")
    plt.grid(True)
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
