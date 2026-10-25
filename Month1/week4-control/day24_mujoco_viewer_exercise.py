"""Day 24 · MuJoCo Viewer + Computed Torque Control (CTC)

Requirements:
    python -m pip install mujoco numpy

Run:
    python day24_mujoco_viewer_exercise.py

Fill in TODO 1–4, then run the script. Designed for MuJoCo 3.17.0.
"""

import time
import numpy as np
import mujoco
import mujoco.viewer

XML = """
<mujoco model="two_link_arm">
    <option timestep="0.002" gravity="0 0 -9.81"/>
    <worldbody>
        <light pos="0 -3 3" dir="0 1 -1"/>
        <geom name="floor" type="plane" pos="0 0 -2.4" size="3 3 0.1" rgba="0.85 0.85 0.85 1"/>
        <body name="link1" pos="0 0 0">
            <joint name="joint1" type="hinge" axis="0 1 0"/>
            <geom name="arm1" type="capsule" fromto="0 0 0 0 0 -1" size="0.05" mass="2" rgba="0.25 0.5 0.85 1"/>
            <body name="link2" pos="0 0 -1">
                <joint name="joint2" type="hinge" axis="0 1 0"/>
                <geom name="arm2" type="capsule" fromto="0 0 0 0 0 -1" size="0.05" mass="1" rgba="0.95 0.55 0.2 1"/>
            </body>
        </body>
    </worldbody>
    <actuator>
        <motor joint="joint1" gear="1"/>
        <motor joint="joint2" gear="1"/>
    </actuator>
</mujoco>
"""


def computed_torque_control(model, data, q_des, qd_des, qdd_des, kp, kd):
    """Exact-model CTC for a fixed-base, directly torque-actuated 2R robot."""
    M = np.zeros((model.nv, model.nv), dtype=np.float64)
    mujoco.mj_fullM(model, data, M)  # MuJoCo 3.17.0 signature

    e = data.qpos - q_des
    ed = data.qvel - qd_des
    v = qdd_des - kp * e - kd * ed
    return M @ v + data.qfrc_bias


def visualize_ctc(model, data, kp=25.0, kd=10.0):
    mujoco.mj_resetData(model, data)
    data.qpos[:] = [0.0, 0.0]
    data.qvel[:] = [0.0, 0.0]

    targets = [
        np.array([0.5, -0.3]),
        np.array([-0.5, 0.6]),
        np.array([0.2, 0.2]),
    ]

    dt = model.opt.timestep
    T = 10.0
    steps = int(T / dt)
    qd_des = np.zeros(2)
    qdd_des = np.zeros(2)

    with mujoco.viewer.launch_passive(model, data) as viewer:
        for k in range(steps):
            if not viewer.is_running():
                break

            wall_start = time.perf_counter()
            t = data.time
            target_index = min(int(t // 3.0), len(targets) - 1)
            q_des = targets[target_index]

            mujoco.mj_forward(model, data)

            # TODO 1: Call computed_torque_control to calculate tau.
            tau = computed_torque_control(model, data, q_des, qd_des, qdd_des, kp, kd)

            # TODO 2: Apply tau to the actuators using data.ctrl.
            data.ctrl[:] = tau
            # TODO 3: Advance the simulation by one timestep.
            mujoco.mj_step(model, data)

            # TODO 4: Synchronize the viewer with the simulation state.
            viewer.sync()

            # Approximate real-time playback, compensating for computation time.
            elapsed = time.perf_counter() - wall_start
            if elapsed < dt:
                time.sleep(dt - elapsed)

            if k % 500 == 0:
                print(f"t={data.time:5.2f}s | target={q_des} | q={data.qpos.copy()}")

    print("Visualization finished.")


def main():
    model = mujoco.MjModel.from_xml_string(XML)
    data = mujoco.MjData(model)
    print(f"MuJoCo {mujoco.__version__} | nq={model.nq}, nv={model.nv}, nu={model.nu}")
    visualize_ctc(model, data)


if __name__ == "__main__":
    main()
