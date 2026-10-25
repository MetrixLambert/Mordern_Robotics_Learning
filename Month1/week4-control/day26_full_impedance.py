"""Day 26 — Full joint-space impedance with desired inertia shaping.

Run:
    python day26_full_impedance.py
    python day26_full_impedance.py --viewer --mode full
    python day26_full_impedance.py --mode pd

Compare 'full' and 'pd', then change MD, DD, KD_STIFF below.
Requires: mujoco 3.x, numpy, matplotlib.
"""
from __future__ import annotations

import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
import mujoco

from day26_joint_impedance import MJCF

# -------- Parameters for the experiment --------
Q_DES = np.array([0.35, -0.50])
QD_DES = np.zeros(2)
QDD_DES = np.zeros(2)
MD = np.diag([0.50, 0.30])       # Desired inertia (kg m^2), NOT actual robot inertia
DD = np.diag([8.0, 6.0])        # Desired damping (N m s/rad)
KD_STIFF = np.diag([100.0, 100.0])  # Desired stiffness (N m/rad)
TAU_EXT = np.array([1.0, 0.0])
T_DISTURB = 2.0
T_END = 8.0
SAMPLE_EVERY = 5
# -----------------------------------------------


def mass_matrix(model: mujoco.MjModel, data: mujoco.MjData) -> np.ndarray:
    """Extract full generalized inertia matrix using MuJoCo 3.x API."""
    M = np.zeros((model.nv, model.nv), dtype=np.float64)
    mujoco.mj_fullM(model, data, M)
    return M


def full_impedance_control(model: mujoco.MjModel, data: mujoco.MjData,
                           tau_ext: np.ndarray) -> np.ndarray:
    """Exact joint-space inertia-shaping impedance, assuming measured tau_ext.

    Desired: MD e_ddot + DD e_dot + KD_STIFF e = tau_ext.
    Dynamics: M q_ddot + h = tau + tau_ext.
    """
    q = data.qpos[:2].copy()
    qdot = data.qvel[:2].copy()
    e = q - Q_DES
    edot = qdot - QD_DES
    M = mass_matrix(model, data)
    h = data.qfrc_bias[:2].copy()  # Coriolis/centrifugal + gravity
    e_ddot_cmd = np.linalg.solve(MD, tau_ext - DD @ edot - KD_STIFF @ e)
    q_ddot_cmd = QDD_DES + e_ddot_cmd
    return M @ q_ddot_cmd + h - tau_ext


def gravity_compensation_pd(model: mujoco.MjModel, data: mujoco.MjData,
                            gravity_data: mujoco.MjData) -> np.ndarray:
    """Baseline: gravity compensation plus PD; physical inertia is unchanged."""
    gravity_data.qpos[:] = data.qpos
    gravity_data.qvel[:] = 0.0
    mujoco.mj_forward(model, gravity_data)
    g = gravity_data.qfrc_bias[:2].copy()
    e = data.qpos[:2] - Q_DES
    edot = data.qvel[:2] - QD_DES
    return g - KD_STIFF @ e - DD @ edot


def control_step(model: mujoco.MjModel, data: mujoco.MjData,
                 gravity_data: mujoco.MjData, mode: str) -> None:
    # mj_forward updates qM and qfrc_bias for the current q and qdot.
    mujoco.mj_forward(model, data)
    tau_ext = TAU_EXT if data.time >= T_DISTURB else np.zeros(2)
    data.qfrc_applied[:] = 0.0
    data.qfrc_applied[:2] = tau_ext
    if mode == "full":
        data.ctrl[:] = full_impedance_control(model, data, tau_ext)
    else:
        data.ctrl[:] = gravity_compensation_pd(model, data, gravity_data)


def simulate(model: mujoco.MjModel, mode: str) -> dict[str, np.ndarray]:
    data = mujoco.MjData(model)
    gravity_data = mujoco.MjData(model)
    data.qpos[:2] = Q_DES
    mujoco.mj_forward(model, data)
    samples = {name: [] for name in ("t", "e", "edot", "tau", "tau_ext", "M11")}
    steps = int(round(T_END / model.opt.timestep))
    for i in range(steps):
        control_step(model, data, gravity_data, mode)
        if i % SAMPLE_EVERY == 0:
            samples["t"].append(float(data.time))
            samples["e"].append((data.qpos[:2] - Q_DES).copy())
            samples["edot"].append(data.qvel[:2].copy())
            samples["tau"].append(data.ctrl.copy())
            samples["tau_ext"].append(data.qfrc_applied[:2].copy())
            samples["M11"].append(mass_matrix(model, data)[0, 0])
        mujoco.mj_step(model, data)
    return {key: np.asarray(value) for key, value in samples.items()}


def plot_comparison(model: mujoco.MjModel, output: str) -> None:
    logs = {mode: simulate(model, mode) for mode in ("pd", "full")}
    theoretical_error = np.linalg.solve(KD_STIFF, TAU_EXT)
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for mode, log in logs.items():
        label = "Full impedance (inertia shaping)" if mode == "full" else "Gravity + PD"
        axes[0].plot(log["t"], log["e"][:, 0], label=label)
        axes[1].plot(log["t"], log["edot"][:, 0], label=label)
        axes[2].plot(log["t"], log["tau"][:, 0], label=label)
        print(f"{mode:4s}: final joint errors = {log['e'][-1]} rad; "
              f"M11 initial/final = {log['M11'][0]:.4f}/{log['M11'][-1]:.4f}")
    axes[0].axhline(theoretical_error[0], linestyle="--", color="black",
                    label=f"Theoretical steady state = {theoretical_error[0]:.4f} rad")
    for ax in axes:
        ax.axvline(T_DISTURB, color="gray", linestyle=":", label="Disturbance on")
        ax.grid(alpha=0.3)
        ax.legend(loc="best")
    axes[0].set_ylabel("Joint 1 position error (rad)")
    axes[1].set_ylabel("Joint 1 velocity (rad/s)")
    axes[2].set_ylabel("Motor 1 torque (N m)")
    axes[2].set_xlabel("Time (s)")
    fig.suptitle("Day 26: Gravity+PD vs Full Inertia-Shaping Impedance")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    print("Theoretical steady-state error:", theoretical_error, "rad")
    print("Saved plot:", output)
    plt.show()


def run_viewer(model: mujoco.MjModel, mode: str) -> None:
    import mujoco.viewer
    data = mujoco.MjData(model)
    gravity_data = mujoco.MjData(model)
    data.qpos[:2] = Q_DES
    mujoco.mj_forward(model, data)
    with mujoco.viewer.launch_passive(model, data) as viewer:
        while viewer.is_running():
            start = time.perf_counter()
            control_step(model, data, gravity_data, mode)
            mujoco.mj_step(model, data)
            viewer.sync()
            remaining = model.opt.timestep - (time.perf_counter() - start)
            if remaining > 0:
                time.sleep(remaining)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["pd", "full"], default="full")
    parser.add_argument("--viewer", action="store_true")
    parser.add_argument("--output", default="day26_full_impedance.png")
    args = parser.parse_args()
    model = mujoco.MjModel.from_xml_string(MJCF)
    if model.nv != 2 or model.nu != 2:
        raise RuntimeError("Expected a 2R arm with two direct-drive motors")
    if True:
        run_viewer(model, args.mode)
    else:
        plot_comparison(model, args.output)


if __name__ == "__main__":
    main()
