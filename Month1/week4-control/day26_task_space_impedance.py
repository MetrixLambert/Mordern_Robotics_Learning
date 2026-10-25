"""Day 26, stage 2: Cartesian impedance of a planar 2R arm in MuJoCo.

Place next to day26_joint_impedance.py. Run:
  python day26_task_space_impedance.py
  python day26_task_space_impedance.py --viewer --mode full
  python day26_task_space_impedance.py --mode pd

The 2D task coordinates are world (x,z), not (x,y): hinge axes are world y.
The disturbance is an ideal, known Cartesian force mapped by J.T to qfrc_applied.
No physical contact, force sensor, actuator saturation or joint-limit contact is modeled.
"""
from __future__ import annotations

import argparse
import time
import numpy as np
import matplotlib.pyplot as plt
import mujoco

from day26_joint_impedance import MJCF

Q_HOME = np.array([0.55, -1.05])
K = np.diag([200.0, 50.0])   # Cartesian stiffness, N/m: x and z
D = np.diag([30.0, 14.0])   # Cartesian damping, N*s/m
MD = np.diag([1.5, 0.8])   # desired Cartesian inertia, kg
F_EXT = np.array([0.0, 5.0])  # world x,z force, N
T_DISTURB = 2.0
T_END = 9.0


def tip_state(model, data, site_id):
    """Return world (x,z) tip position, velocity and 2x2 Jacobian."""
    jp = np.zeros((3, model.nv))
    jr = np.zeros((3, model.nv))
    mujoco.mj_jacSite(model, data, jp, jr, site_id)
    J = jp[[0, 2], :]
    x = np.asarray(data.site_xpos[site_id])[[0, 2]].copy()
    xd = J @ data.qvel
    return x, xd, J


def jacobian_dot_qdot(model, data, site_id):
    """Analytical Jdot*qdot for this specific 2R x-z arm.

    With q1,q2 rotating about +y, link lengths l1=0.65,l2=0.55:
      x=l1*cos(q1)+l2*cos(q1+q2)
      z=z0-l1*sin(q1)-l2*sin(q1+q2)
    """
    q1, q2 = data.qpos[:2]
    v1, v2 = data.qvel[:2]
    l1, l2 = 0.65, 0.55
    a = v1 + v2
    return np.array([
        -l1*np.cos(q1)*v1*v1 - l2*np.cos(q1+q2)*a*a,
         l1*np.sin(q1)*v1*v1 + l2*np.sin(q1+q2)*a*a,
    ])


def make_target(model, site_id):
    reference = mujoco.MjData(model)
    reference.qpos[:2] = Q_HOME
    mujoco.mj_forward(model, reference)
    return np.asarray(reference.site_xpos[site_id])[[0, 2]].copy()


def controller(model, data, site_id, x_des, mode, force):
    """Full inertia-shaped impedance, or simpler Cartesian PD + gravity compensation.

    Dynamics: M qdd + h = tau + J.T F_ext.
    Full control enforces MD e_dd + D e_d + K e = F_ext.
    """
    mujoco.mj_forward(model, data)
    x, xd, J = tip_state(model, data, site_id)
    e = x - x_des
    if abs(np.linalg.det(J)) < 1e-4:
        raise RuntimeError("Jacobian near singular; choose a different Q_HOME")
    if mode == "full":
        M = np.zeros((model.nv, model.nv), dtype=np.float64)
        mujoco.mj_fullM(model, data, M)
        # qdd = J^-1 (xdd_cmd - Jdot*qdot)
        xdd_cmd = np.linalg.solve(MD, force - D @ xd - K @ e)
        qdd_cmd = np.linalg.solve(J, xdd_cmd - jacobian_dot_qdot(model, data, site_id))
        tau = M @ qdd_cmd + data.qfrc_bias.copy() - J.T @ force
    elif mode == "pd":
        # Cartesian spring-damper force plus gravity compensation (not inertia shaping).
        gravity_data = mujoco.MjData(model)
        gravity_data.qpos[:] = data.qpos
        gravity_data.qvel[:] = 0
        mujoco.mj_forward(model, gravity_data)
        tau = gravity_data.qfrc_bias.copy() + J.T @ (-K @ e - D @ xd)
    else:
        raise ValueError(mode)
    return tau, x, xd, e, J


def step_control(model, data, site_id, x_des, mode):
    force = F_EXT if data.time >= T_DISTURB else np.zeros(2)
    tau, x, xd, e, J = controller(model, data, site_id, x_des, mode, force)
    data.ctrl[:] = tau
    data.qfrc_applied[:] = J.T @ force
    return x, xd, e, force


def simulate(model, mode):
    data = mujoco.MjData(model)
    data.qpos[:2] = Q_HOME
    mujoco.mj_forward(model, data)
    site_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, "tip")
    x_des = make_target(model, site_id)
    log = {key: [] for key in ("t", "e", "v", "f", "tau")}
    for i in range(round(T_END / model.opt.timestep)):
        _, xd, e, force = step_control(model, data, site_id, x_des, mode)
        if i % 5 == 0:
            for key, value in (("t", data.time), ("e", e), ("v", xd),
                               ("f", force), ("tau", data.ctrl.copy())):
                log[key].append(np.array(value, copy=True))
        mujoco.mj_step(model, data)
    return {key: np.asarray(value) for key, value in log.items()}


def plot(log, output, mode):
    t = log["t"]
    fig, axes = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
    for i, axis in enumerate(("x", "z")):
        axes[0].plot(t, log["e"][:, i], label=f"{axis} error")
        axes[0].hlines(np.linalg.solve(K, F_EXT)[i], T_DISTURB, T_END,
                       color=f"C{i}", ls="--", alpha=.7,
                       label=f"{axis} theoretical equilibrium")
        axes[1].plot(t, log["v"][:, i], label=f"{axis} velocity")
        axes[2].plot(t, log["f"][:, i], label=f"{axis} external force")
    for ax in axes:
        ax.axvline(T_DISTURB, color="gray", ls=":")
        ax.grid(alpha=.3)
        ax.legend()
    axes[0].set_ylabel("Cartesian error (m)")
    axes[1].set_ylabel("Tip velocity (m/s)")
    axes[2].set_ylabel("External force (N)")
    axes[2].set_xlabel("Time (s)")
    fig.suptitle(f"Day 26 - Task-space impedance ({mode})")
    fig.tight_layout()
    fig.savefig(output, dpi=160)
    print("Theoretical equilibrium error [x,z] (m):", np.linalg.solve(K, F_EXT))
    print("Final simulated error [x,z] (m):       ", log["e"][-1])
    print("Saved:", output)
    plt.show()


def viewer(model, mode):
    import mujoco.viewer
    data = mujoco.MjData(model)
    data.qpos[:2] = Q_HOME
    mujoco.mj_forward(model, data)
    site_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, "tip")
    x_des = make_target(model, site_id)
    with mujoco.viewer.launch_passive(model, data) as v:
        while v.is_running():
            start = time.perf_counter()
            step_control(model, data, site_id, x_des, mode)
            mujoco.mj_step(model, data)
            v.sync()
            time.sleep(max(0, model.opt.timestep - (time.perf_counter()-start)))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("full", "pd"), default="full")
    parser.add_argument("--viewer", action="store_true")
    parser.add_argument("--output", default="day26_task_space_impedance.png")
    args = parser.parse_args()
    model = mujoco.MjModel.from_xml_string(MJCF)
    if model.nv != 2 or model.nu != 2:
        raise RuntimeError("Expected a two-joint, two-motor model")
    if args.viewer:
        viewer(model, args.mode)
    else:
        plot(simulate(model, args.mode), args.output, args.mode)


if __name__ == "__main__":
    main()
