"""Day 10 / Experiment 10.2: validate MuJoCo dynamics and foot Jacobian.

Place this file next to the existing planar_leg.xml, then run:
    python run_10_2.py

No controller and no viewer are needed for this diagnostic experiment.
"""
from pathlib import Path
import numpy as np
import mujoco

np.set_printoptions(precision=5, suppress=True, linewidth=125)

MODEL_PATH = Path(__file__).with_name("planar_leg.xml")
model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
data = mujoco.MjData(model)

foot_site_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, "foot_center")
key_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_KEY, "initial")
assert foot_site_id >= 0 and key_id >= 0
assert (model.nq, model.nv, model.nu) == (5, 5, 2), "Unexpected robot dimensions"


def reset_initial():
    mujoco.mj_resetDataKeyframe(model, data, key_id)
    data.ctrl[:] = 0.0
    mujoco.mj_forward(model, data)


def foot_pose_and_jacobian(d):
    """Planar foot pose [x,z,phi] and Jacobian in Day-9 sign convention.

    MuJoCo rotational hinge axes in planar_leg.xml point along negative Y,
    so phi_dot = -omega_y (not +omega_y).
    """
    jacp = np.zeros((3, model.nv))
    jacr = np.zeros((3, model.nv))
    mujoco.mj_jacSite(model, d, jacp, jacr, foot_site_id)
    pos = d.site_xpos[foot_site_id]
    # All 3 hinges use negative-Y axes; their accumulated angle is q[2]+q[3]+q[4].
    phi = float(d.qpos[2] + d.qpos[3] + d.qpos[4])
    pose = np.array([pos[0], pos[2], phi])
    J = np.vstack((jacp[0], jacp[2], -jacr[1]))
    return pose, J


def numerical_jacobian(qpos, eps=1e-6):
    """Central differences of the foot pose (valid here: all joints are 1-DoF)."""
    scratch = mujoco.MjData(model)
    J_fd = np.zeros((3, model.nv))
    for j in range(model.nv):
        poses = []
        for sign in (+1, -1):
            scratch.qpos[:] = qpos
            scratch.qpos[j] += sign * eps
            mujoco.mj_forward(model, scratch)
            poses.append(foot_pose_and_jacobian(scratch)[0])
        J_fd[:, j] = (poses[0] - poses[1]) / (2 * eps)
    return J_fd


def get_mass_matrix(d):
    """Reconstruct M from matrix-vector products; does not access data.qM."""
    nv = model.nv
    M = np.empty((nv, nv))
    for j in range(nv):
        e = np.zeros(nv)
        e[j] = 1.0
        col = np.zeros(nv)
        mujoco.mj_mulM(model, d, col, e)
        M[:, j] = col
    return M


def diagnose(label):
    # mj_forward calculates qM, qfrc_bias, qacc, contact constraints, etc.
    mujoco.mj_forward(model, data)
    M = get_mass_matrix(data)
    h = data.qfrc_bias.copy()
    Q_act = data.qfrc_actuator.copy()
    Q_passive = data.qfrc_passive.copy()
    Q_constraint = data.qfrc_constraint.copy()
    a = data.qacc.copy()

    foot_pose, J = foot_pose_and_jacobian(data)
    J_fd = numerical_jacobian(data.qpos.copy())
    residual = M @ a + h - Q_act - Q_passive - Q_constraint

    print("\n" + "=" * 74)
    print(label, f"| t={data.time:.3f} s | contacts={data.ncon}")
    print("qpos =", data.qpos)
    print("qvel =", data.qvel)
    print("foot [x,z,phi] =", foot_pose)
    print("\nM(q) =\n", M)
    print("symmetry error ||M-M.T||_inf =", np.linalg.norm(M - M.T, ord=np.inf))
    print("eigenvalues(M) =", np.linalg.eigvalsh(M))
    print("\nh = qfrc_bias =", h)
    print("Q_passive =", Q_passive)
    print("Q_actuator =", Q_act)
    print("Q_constraint =", Q_constraint)
    print("\nJ_foot =\n", J)
    print("Jacobian finite-difference error (max abs) =", np.max(np.abs(J - J_fd)))
    print("\nqacc =", a)
    print("dynamics residual =", residual)
    print("dynamics residual norm =", np.linalg.norm(residual))

    assert np.linalg.norm(M - M.T, ord=np.inf) < 1e-9, "M not symmetric"
    assert np.min(np.linalg.eigvalsh(M)) > 0, "M not positive definite"
    assert np.max(np.abs(J - J_fd)) < 1e-5, "Foot Jacobian mismatch"
    assert np.linalg.norm(residual) < 1e-5, "Dynamics identity mismatch"


print("MuJoCo version:", mujoco.__version__)
print("nq, nv, nu:", model.nq, model.nv, model.nu)
print("joint order:", [mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_JOINT, i)
                       for i in range(model.njnt)])

# Case A: XML initial keyframe (foot bottom nominally touches the ground).
reset_initial()
diagnose("A. Initial pose")

# Case B: lift the whole robot off the ground and give it nonzero velocity.
reset_initial()
data.qpos[1] += 0.30  # base_z; leave all joint angles unchanged
data.qvel[:] = [0.2, -0.1, 0.3, -0.4, 0.2]
diagnose("B. Airborne pose with nonzero velocity")

# Case C: let it fall, then inspect the first state with active contacts.
reset_initial()
contact_found = False
for _ in range(500):
    mujoco.mj_step(model, data)
    if data.ncon > 0:
        diagnose("C. Contact during unactuated simulation")
        contact_found = True
        break
if not contact_found:
    print("\nNo contacts found within 1 s. Check the XML/initial pose.")

print("\nAll checked matrix/Jacobian/dynamics identities passed.")
