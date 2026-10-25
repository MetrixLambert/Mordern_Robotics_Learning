"""Day 26: MuJoCo 2R joint-space impedance exercise.

Dependencies: numpy, matplotlib, mujoco (tested API pattern for MuJoCo 3.x).
Run: python day26_joint_impedance.py
Optional interactive viewer: python day26_joint_impedance.py --viewer

Exercise: modify KP, KD, disturbance and compare steady-state errors.
"""
from __future__ import annotations

import argparse
import numpy as np
import matplotlib.pyplot as plt
import mujoco

MJCF = r"""
<mujoco model="day26_2r_impedance">
  <compiler angle="radian"/>
  <option timestep="0.002" gravity="0 0 -9.81" integrator="RK4"/>
  <default>
    <joint damping="0" armature="0"/>
    <geom type="capsule" size="0.045" density="850" rgba="0.22 0.48 0.85 1"/>
  </default>
  <worldbody>
    <light pos="0 -2 3"/>
    <geom type="plane" size="2 2 0.1" rgba="0.9 0.9 0.9 1"/>
    <body name="link1" pos="0 0 1.1">
      <joint name="joint1" type="hinge" axis="0 1 0" range="-2.8 2.8" limited="true"/>
      <geom name="upper" fromto="0 0 0 0.65 0 0"/>
      <body name="link2" pos="0.65 0 0">
        <joint name="joint2" type="hinge" axis="0 1 0" range="-2.8 2.8" limited="true"/>
        <geom name="lower" fromto="0 0 0 0.55 0 0" rgba="0.20 0.70 0.52 1"/>
        <site name="tip" pos="0.55 0 0" size="0.035" rgba="1 0.2 0.2 1"/>
      </body>
    </body>
  </worldbody>
  <actuator>
    <motor name="motor1" joint="joint1" gear="1"/>
    <motor name="motor2" joint="joint2" gear="1"/>
  </actuator>
</mujoco>
"""

# -------- Experiment settings (change these first) --------
Q_DES = np.array([0.35, -0.50], dtype=float)
KP = np.diag([100.0, 100.0])  # N*m/rad
KD = np.diag([10.0, 10.0])    # N*m*s/rad
TAU_EXT = np.array([1.0, 0.0])  # N*m; applied after T_DISTURB
T_DISTURB = 2.0
T_END = 7.0
# ---------------------------------------------------------


def gravity_compensation_pd(model: mujoco.MjModel, data: mujoco.MjData) -> np.ndarray:
    """tau = g(q) - KP (q-qd) - KD qdot.

    Gravity g(q) is obtained by zeroing qvel and qacc in a separate data object.
    This avoids incorrectly using qfrc_bias as pure gravity when the arm moves.
    """
    gravity_data = mujoco.MjData(model)
    gravity_data.qpos[:] = data.qpos
    gravity_data.qvel[:] = 0.0
    mujoco.mj_forward(model, gravity_data)
    gravity = gravity_data.qfrc_bias.copy()
    e = data.qpos[:2] - Q_DES
    edot = data.qvel[:2]
    return gravity - KP @ e - KD @ edot


def apply_control(model: mujoco.MjModel, data: mujoco.MjData) -> None:
    data.ctrl[:] = gravity_compensation_pd(model, data)
    data.qfrc_applied[:] = 0.0
    if data.time >= T_DISTURB:
        data.qfrc_applied[:2] = TAU_EXT


def simulate(model: mujoco.MjModel) -> dict[str, np.ndarray]:
    data = mujoco.MjData(model)
    data.qpos[:2] = Q_DES
    mujoco.mj_forward(model, data)
    samples = {key: [] for key in ("t", "q", "qd", "tau", "tau_ext")}
    n_steps = int(round(T_END / model.opt.timestep))
    for step in range(n_steps):
        apply_control(model, data)
        if step % 5 == 0:
            samples["t"].append(float(data.time))
            samples["q"].append(data.qpos[:2].copy())
            samples["qd"].append(data.qvel[:2].copy())
            samples["tau"].append(data.ctrl.copy())
            samples["tau_ext"].append(data.qfrc_applied[:2].copy())
        mujoco.mj_step(model, data)
    return {k: np.asarray(v) for k, v in samples.items()}


def plot_results(log: dict[str, np.ndarray], output: str) -> None:
    t = log["t"]
    e = log["q"] - Q_DES
    e_theory = np.linalg.solve(KP, TAU_EXT)
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for i in range(2):
        axes[0].plot(t, e[:, i], label=f"joint {i+1} error")
        axes[0].hlines(e_theory[i], T_DISTURB, T_END,
                       linestyles="--", colors=f"C{i}", alpha=0.65,
                       label=f"joint {i+1} theory: {e_theory[i]:.4f} rad")
        axes[1].plot(t, log["qd"][:, i], label=f"joint {i+1}")
        axes[2].plot(t, log["tau_ext"][:, i], label=f"joint {i+1}")
    axes[0].set_ylabel("Position error (rad)")
    axes[1].set_ylabel("Joint velocity (rad/s)")
    axes[2].set_ylabel("External torque (N m)")
    axes[2].set_xlabel("Time (s)")
    for ax in axes:
        ax.axvline(T_DISTURB, color="gray", ls=":", alpha=0.8)
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best")
    fig.suptitle("Day 26 - Joint-space impedance under external torque")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    print(f"Saved plot: {output}")
    print(f"Theoretical steady-state error: {e_theory} rad")
    print(f"Simulated final error:         {e[-1]} rad")
    plt.show()


def run_viewer(model: mujoco.MjModel) -> None:
    import time
    import mujoco.viewer
    data = mujoco.MjData(model)
    data.qpos[:2] = Q_DES
    mujoco.mj_forward(model, data)
    with mujoco.viewer.launch_passive(model, data) as viewer:
        while viewer.is_running():
            start = time.perf_counter()
            apply_control(model, data)
            mujoco.mj_step(model, data)
            viewer.sync()
            remaining = model.opt.timestep - (time.perf_counter() - start)
            if remaining > 0:
                time.sleep(remaining)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--viewer", action="store_true", help="Show interactive MuJoCo viewer")
    parser.add_argument("--output", default="day26_joint_impedance.png")
    args = parser.parse_args()
    model = mujoco.MjModel.from_xml_string(MJCF)
    if model.nv != 2 or model.nu != 2:
        raise RuntimeError("Expected two hinge DOFs and two motors")
    if args.viewer == "True":
        run_viewer(model)
    else:
        plot_results(simulate(model), args.output)


if __name__ == "__main__":
    main()
