"""ROBO_Learn - Month 2 / Week 1 / Day 3: friction contact forces.

Requirements: pip install mujoco numpy matplotlib
"""
import mujoco
import numpy as np
import matplotlib.pyplot as plt


def run_simulation(mu: float, force_x: float = 4.0, duration: float = 2.0):
    xml = f"""
    <mujoco model="friction_test">
      <option timestep="0.002" gravity="0 0 -9.81"/>
      <worldbody>
        <geom name="ground" type="plane" size="10 10 0.1"
              friction="{mu} 0.005 0.0001"/>
        <body name="slider" pos="0 0 0.1">
          <joint name="slide_x" type="slide" axis="1 0 0"/>
          <joint name="slide_z" type="slide" axis="0 0 1"/>
          <geom name="box" type="box" size="0.1 0.1 0.1" mass="1"
                friction="{mu} 0.005 0.0001"/>
        </body>
      </worldbody>
    </mujoco>
    """
    model = mujoco.MjModel.from_xml_string(xml)
    data = mujoco.MjData(model)
    body_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_BODY, "slider")
    box_geom_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM, "box")
    ground_geom_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM, "ground")

    times, positions, velocities = [], [], []
    normal_forces, tangential_forces, contact_counts = [], [], []

    while data.time < duration - 1e-12:
        data.xfrc_applied[body_id, :] = 0.0
        data.xfrc_applied[body_id, 0] = force_x
        mujoco.mj_step(model, data)

        force_on_box_world = np.zeros(3)
        count = 0
        for i in range(data.ncon):
            contact = data.contact[i]
            if {contact.geom1, contact.geom2} != {ground_geom_id, box_geom_id}:
                continue
            force_local = np.zeros(6)
            mujoco.mj_contactForce(model, data, i, force_local)
            # MuJoCo contact.frame rows are contact-frame basis vectors in world coordinates.
            force_world = contact.frame.reshape(3, 3).T @ force_local[:3]
            # mj_contactForce returns force on geom2 by geom1.
            if contact.geom2 == box_geom_id:
                force_on_box_world += force_world
            else:
                force_on_box_world -= force_world
            count += 1

        times.append(data.time)
        positions.append(data.qpos[0])
        velocities.append(data.qvel[0])
        normal_forces.append(force_on_box_world[2])
        tangential_forces.append(force_on_box_world[0])
        contact_counts.append(count)

    return {key: np.asarray(value) for key, value in {
        "t": times, "x": positions, "vx": velocities,
        "fn": normal_forces, "ft": tangential_forces, "ncon": contact_counts,
    }.items()}


def main():
    force_x = 4.0
    results = {mu: run_simulation(mu, force_x=force_x) for mu in (0.8, 0.2)}
    fig, axes = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    for mu, r in results.items():
        label = f"mu = {mu}"
        axes[0].plot(r["t"], r["x"], label=label)
        axes[1].plot(r["t"], r["vx"], label=label)
        axes[2].plot(r["t"], r["fn"], label=label)
        axes[3].plot(r["t"], r["ft"], label=label)

    axes[0].set_ylabel("Horizontal position (m)")
    axes[1].set_ylabel("Horizontal velocity (m/s)")
    axes[2].set_ylabel("Normal force (N)")
    axes[3].set_ylabel("Tangential force (N)")
    axes[3].set_xlabel("Time (s)")
    axes[2].axhline(9.81, color="gray", linestyle="--", label="mg")
    axes[3].axhline(-4.0, color="gray", linestyle="--", label="static balance: -4 N")
    axes[3].axhline(-0.2 * 9.81, color="black", linestyle=":", label="sliding: -mu*mg (mu=0.2)")
    for ax in axes:
        ax.grid(True)
        ax.legend()
    fig.tight_layout()
    plt.show()

    print("\nFinal simulation values (mean over final 0.5 s):")
    for mu, r in results.items():
        mask = r["t"] >= 1.5
        fn = float(np.mean(r["fn"][mask]))
        ft = float(np.mean(r["ft"][mask]))
        ratio = abs(ft) / fn if fn > 1e-8 else float("nan")
        print(f"mu={mu:.1f}: x_end={r['x'][-1]:.4f} m, "
              f"vx_end={r['vx'][-1]:.4f} m/s, "
              f"mean_fn={fn:.4f} N, mean_ft={ft:.4f} N, "
              f"|ft|/fn={ratio:.4f}, last_ncon={r['ncon'][-1]}")


if __name__ == "__main__":
    main()
