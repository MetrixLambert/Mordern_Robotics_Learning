"""ROBO_Learn / Month 2 / Day 18: LIPM, DCM and Capture Point.

运行：python day18_lipm_practice.py

依赖：numpy, scipy, matplotlib

练习约定：
  x 为水平质心位置 [m]，v 为水平质心速度 [m/s]；
  p 为水平 ZMP 位置 [m]，zc 为恒定质心高度 [m]。
  假设 p 可瞬间切换：这是理想化模型，不模拟实际迈步过程。

完成 TODO 1–4 后运行。不要修改底部的检查与绘图代码。
"""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Config:
    g: float = 9.81
    zc: float = 1.0
    x0: float = 0.0
    v0: float = 0.5
    touch_time: float = 0.5
    end_time: float = 2.0
    dt: float = 0.005


def omega(cfg: Config) -> float:
    """TODO 1: 实现 LIPM 的 omega，单位 [1/s]。"""
    omega = np.sqrt(cfg.g / cfg.zc)
    return omega


def dcm(x: float, v: float, cfg: Config) -> float:
    """TODO 2: 根据 (x, v) 计算 DCM / instantaneous capture point [m]。"""
    xi = x + v / omega(cfg)
    return xi


def lipm_rhs(t: float, state: np.ndarray, p: float, cfg: Config) -> np.ndarray:
    """TODO 3: 返回 [dx/dt, dv/dt]，其中 state=[x, v]。

    提示：p 是 ZMP，不是驱动力；公式为 x_ddot = omega^2 * (x-p)。
    """
    x, v = state
    x_ddot = omega(cfg)**2 * (x - p) 
    return np.array([v, x_ddot])


def predict_touchdown_dcm(x0: float, v0: float, p0: float,
                          T: float, cfg: Config) -> float:
    """TODO 4: 假设接触前 ZMP=p0 恒定，解析预测 T 秒后的 DCM。

    提示：先计算 xi(0)，再求解 xi_dot=omega*(xi-p0)。
    不要通过数值积分来实现此函数；用于与积分结果交叉验证。
    """
    xi0 = dcm(x0, v0, cfg) 
    w = omega(cfg)
    return p0 + (xi0 - p0) * np.exp(w * T)


def simulate_case(cfg: Config, p_after: float, times: np.ndarray) -> np.ndarray:
    """接触前 p=0；在 touch_time 时理想切换为 p_after。"""
    split = int(np.searchsorted(times, cfg.touch_time))
    t1 = times[:split + 1]
    t2 = times[split:]
    initial = np.array([cfg.x0, cfg.v0], dtype=float)

    sol1 = solve_ivp(
        lipm_rhs, (t1[0], t1[-1]), initial,
        args=(0.0, cfg), t_eval=t1, rtol=1e-10, atol=1e-12,
    )
    if not sol1.success:
        raise RuntimeError(sol1.message)
    sol2 = solve_ivp(
        lipm_rhs, (t2[0], t2[-1]), sol1.y[:, -1],
        args=(p_after, cfg), t_eval=t2, rtol=1e-10, atol=1e-12,
    )
    if not sol2.success:
        raise RuntimeError(sol2.message)
    return np.concatenate([sol1.y[:, :-1], sol2.y], axis=1)


def main() -> None:
    cfg = Config()
    times = np.arange(0.0, cfg.end_time + cfg.dt / 2, cfg.dt)
    assert np.isclose(times, np.round(times / cfg.dt) * cfg.dt).all()
    assert np.any(np.isclose(times, cfg.touch_time))

    xi0 = dcm(cfg.x0, cfg.v0, cfg)
    p_capture = predict_touchdown_dcm(cfg.x0, cfg.v0, 0.0, cfg.touch_time, cfg)
    print(f"omega={omega(cfg):.5f}  xi(0)={xi0:.5f} m")
    print(f"predicted touchdown capture point={p_capture:.5f} m")

    # 数值验证：触地前的 DCM 必须等于解析预测。
    before = solve_ivp(
        lipm_rhs, (0, cfg.touch_time), [cfg.x0, cfg.v0],
        args=(0.0, cfg), rtol=1e-10, atol=1e-12,
    )
    xi_numeric = dcm(*before.y[:, -1], cfg)
    np.testing.assert_allclose(xi_numeric, p_capture, atol=1e-7)

    # 三种触地后的理想情形：精确捕获、落脚偏后、落脚偏前。
    cases = {
        "capture": p_capture,
        "behind CP (-0.08 m)": p_capture - 0.08,
        "ahead of CP (+0.08 m)": p_capture + 0.08,
    }
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    for name, p_next in cases.items():
        path = simulate_case(cfg, p_next, times)
        axes[0].plot(times, path[0], label=name)
        axes[1].plot(times, path[1], label=name)
        print(f"{name:25s} p={p_next:.4f} m; "
              f"x(end)={path[0,-1]:+.4f} m; v(end)={path[1,-1]:+.4f} m/s")

    for ax in axes:
        ax.axvline(cfg.touch_time, color="gray", ls="--", label="touchdown")
        ax.grid(True, alpha=0.3)
        ax.legend()
    axes[0].set_ylabel("CoM x [m]")
    axes[1].set_ylabel("CoM v [m/s]")
    axes[1].set_xlabel("time [s]")
    fig.suptitle("Day 18: LIPM with ideal instantaneous ZMP switch")
    fig.tight_layout()
    plt.show()

    print("\n思考题：")
    print("1. 为什么三种方案在 touchdown 前的轨迹完全一致？")
    print("2. 精确捕获时，触地后的 x 与 v 最终趋向哪里？")
    print("3. 若触地延迟 0.1 s，而仍采用旧捕获点，会出现什么趋势？")
    print("4. 解释为什么此练习不能直接代表真实机器人的迈步过程。")


if __name__ == "__main__":
    main()
