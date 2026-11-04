"""ROBO_Learn | Day 4.5 — MuJoCo kinematics practice

Independent exercise: body/site poses, Jacobians, and end-effector velocity.
Run: python day46_mujoco_kinematics.py

Fill in TODO 1–5 yourself. No PD controller is needed.
"""

import numpy as np
import mujoco


MODEL_XML = """
<mujoco model="day46_kinematics">
    <option timestep="0.002" gravity="0 0 -9.81"/>
    <worldbody>
        <body name="link1" pos="0 0 0">
            <joint name="joint1" type="hinge" axis="0 1 0"/>
            <geom name="link1_geom" type="capsule"
                  fromto="0 0 0 0 0 -1" size="0.04" mass="1"/>
            <body name="link2" pos="0 0 -1">
                <joint name="joint2" type="hinge" axis="0 1 0"/>
                <geom name="link2_geom" type="capsule"
                      fromto="0 0 0 0 0 -0.8" size="0.04" mass="0.8"/>
                <site name="ee" pos="0 0 -0.8" size="0.025" rgba="1 0 0 1"/>
            </body>
        </body>
    </worldbody>
</mujoco>
"""


def main():
    model = mujoco.MjModel.from_xml_string(MODEL_XML)
    data = mujoco.MjData(model)

    # Fixed configuration for kinematics experiments.
    q = np.array([0.5, -0.3])
    dq = np.array([0.2, -0.4])
    data.qpos[:] = q
    data.qvel[:] = dq
    mujoco.mj_forward(model, data)

    print("nq, nv:", model.nq, model.nv)
    print("q:", data.qpos)
    print("dq:", data.qvel)

    # TODO 1: Read the world-frame positions of link1, link2 and ee.
    # Hint: data.body("link1").xpos, data.site("ee").xpos
    # Print all three and compare with the geometry in MODEL_XML.
    print("link1 xpos:", data.body("link1").xpos)
    print("link2 xpos:", data.body("link2").xpos)
    print("ee xpos:", data.site("ee").xpos)

    # TODO 2: Read the world-frame orientation of link2 and ee.
    # Hint: data.body("link2").xmat (flat or 3x3 depending on API view)
    #       data.body("link2").xquat
    # Question: What is the relationship between the two orientations?
    print("link2 xmat:", data.body("link2").xmat)
    print("link2 xquat:", data.body("link2").xquat)

    # TODO 3: Calculate the theoretical ee position using planar FK.
    # x = -l1*sin(q1) - l2*sin(q1+q2)
    # y = 0
    # z = -l1*cos(q1) - l2*cos(q1+q2)
    # Compare to the ee world position; print norm of the difference.
    x = -0.04 * np.sin(q[0]) - 0.04 * np.sin(q[0] + q[1])
    y = 0
    z = -0.04 * np.cos(q[0]) - 0.04 * np.cos(q[0] + q[1])
    print("Theoretical ee position:", np.array([x, y, z]))
    print("Difference:", np.linalg.norm(data.site("ee").xpos - np.array([x, y, z])))


    # TODO 4: Obtain translational and rotational site Jacobians.
    # Hint: jacp = np.zeros((3, model.nv))
    #       jacr = np.zeros((3, model.nv))
    #       mujoco.mj_jacSite(model, data, jacp, jacr, site_id)
    #       site_id = model.site("ee").id
    # Print jacp and jacr, and calculate v_ee = jacp @ dq.
    jacp = np.zeros((3, model.nv))
    jacr = np.zeros((3, model.nv))
    mujoco.mj_jacSite(model, data, jacp, jacr, model.site("ee").id)
    print("jacp:", jacp)
    print("jacr:", jacr)
    print("v_ee:", jacp @ dq)

    # TODO 5: Independently verify linear velocity using finite differences.
    # Choose a small dt (e.g. 1e-6), set q_next = q + dq*dt,
    # call mj_forward at q_next, and calculate:
    # v_fd = (p_next - p_now) / dt
    # IMPORTANT: save p_now with .copy() before modifying qpos.
    # Print v_fd, jacp @ dq, and their difference.
    # This is a *kinematic* finite difference; do not call mj_step.
    dt = 1e-6
    q_next = q + dq * dt
    data.qpos[:] = q_next
    data.qvel[:] = dq
    mujoco.mj_forward(model, data)
    p_next = data.site("ee").xpos.copy()
    v_fd = (p_next - data.site("ee").xpos) / dt
    print("v_fd:", v_fd)
    print("Difference:", np.linalg.norm(v_fd - jacp @ dq))
    print("\nComplete TODO 1–5, then compare your outputs.")


if __name__ == "__main__":
    main()
