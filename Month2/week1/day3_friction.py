
import mujoco
import numpy as np
import matplotlib.pyplot as plt

def run_simulation(mu, force_x=4.0, duration=2.0):
    xml = f"""
    <mujoco model="friction_test">
        <option timestep="0.002" gravity="0 0 -9.81"/>

        <worldbody>
            <geom name="ground"
                  type="plane"
                  size="3 3 0.1"
                  friction="{mu} 0.005 0.0001"/>

            <body name="slider" pos="0 0 0.1">
                <joint name="slide_x" type="slide" axis="1 0 0"/>
                <joint name="slide_z" type="slide" axis="0 0 1"/>

                <geom name="box"
                      type="box"
                      size="0.1 0.1 0.1"
                      mass="1"
                      friction="{mu} 0.005 0.0001"/>
            </body>
        </worldbody>
    </mujoco>
    """

    model = mujoco.MjModel.from_xml_string(xml)
    data = mujoco.MjData(model)

    body_id = mujoco.mj_name2id(
        model, mujoco.mjtObj.mjOBJ_BODY, "slider"
    )

    times = []
    positions = []
    velocities = []

    while data.time < duration:
        # 施加水平外力
        data.xfrc_applied[body_id, 0] = force_x

        mujoco.mj_step(model, data)

        times.append(data.time)
        positions.append(data.qpos[0])
        velocities.append(data.qvel[0])

    return np.array(times), np.array(positions), np.array(velocities)


results = {}

for mu in [0.8, 0.2]:
    results[mu] = run_simulation(mu)

fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

for mu, (t, x, vx) in results.items():
    axes[0].plot(t, x, label=f"mu = {mu}")
    axes[1].plot(t, vx, label=f"mu = {mu}")

axes[0].set_ylabel("Horizontal position (m)")
axes[1].set_ylabel("Horizontal velocity (m/s)")
axes[1].set_xlabel("Time (s)")

for ax in axes:
    ax.grid(True)
    ax.legend()

plt.tight_layout()
plt.show()
