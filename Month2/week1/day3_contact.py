
import mujoco
import numpy as np
import matplotlib.pyplot as plt

XML = """
<mujoco model="contact_day3">
    <option timestep="0.002" gravity="0 0 -9.81"/>

    <worldbody>
        <geom name="ground"
              type="plane"
              size="2 2 0.1"
              friction="0.8 0.005 0.0001"/>

        <body name="box" pos="0 0 0.5">
            <freejoint/>
            <geom name="box_geom"
                  type="box"
                  size="0.1 0.1 0.1"
                  mass="1"
                  friction="0.8 0.005 0.0001"/>
        </body>
    </worldbody>
</mujoco>
"""

model = mujoco.MjModel.from_xml_string(XML)
data = mujoco.MjData(model)

times = []
heights = []
velocities = []
contact_counts = []
normal_forces = []

while data.time < 2.0:
    mujoco.mj_step(model, data)

    total_normal_force = 0.0

    for i in range(data.ncon):
        contact = data.contact[i]

        # Contact frame 中第一个方向为接触法向
        force = np.zeros(6)
        mujoco.mj_contactForce(model, data, i, force)

        total_normal_force += force[0]

    times.append(data.time)
    heights.append(data.qpos[2])
    velocities.append(data.qvel[2])
    contact_counts.append(data.ncon)
    normal_forces.append(total_normal_force)

fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

axes[0].plot(times, heights)
axes[0].set_ylabel("Center height (m)")

axes[1].plot(times, velocities)
axes[1].set_ylabel("Vertical velocity (m/s)")

axes[2].plot(times, normal_forces, label="Contact normal force")
axes[2].axhline(9.81, linestyle="--", label="mg = 9.81 N")
axes[2].set_ylabel("Normal force (N)")
axes[2].set_xlabel("Time (s)")
axes[2].legend()

for ax in axes:
    ax.grid(True)

plt.tight_layout()
plt.show()

print("Final center height:", heights[-1])
print("Final vertical velocity:", velocities[-1])
print("Final contact count:", contact_counts[-1])
print("Final total normal force:", normal_forces[-1])
