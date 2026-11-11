"""Day 10 / 10.4: compare planned QP wrench with MuJoCo contact wrench.

Put next to planar_leg.xml. Requires mujoco, cvxpy, osqp, matplotlib, numpy.
The actual contact wrench is sampled at the END of each 10 ms control interval,
so it is not simultaneous with the QP wrench computed at its beginning.
"""
from pathlib import Path
import numpy as np
import mujoco
import cvxpy as cp
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
model = mujoco.MjModel.from_xml_path(str(ROOT / 'planar_leg.xml'))
data = mujoco.MjData(model)
assert (model.nq, model.nv, model.nu) == (5, 5, 2)
site_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, 'foot_center')
foot_geom_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM, 'foot_geom')
ground_geom_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM, 'ground')
key_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_KEY, 'initial')
assert min(site_id, foot_geom_id, ground_geom_id, key_id) >= 0

MU, HALF_LENGTH, TAU_LIMIT = 0.8, 0.12, 200.0
DT_CONTROL, T_FINAL = 0.01, 3.0
KP_BASE = np.array([35.0, 55.0, 55.0])
KD_BASE = 2 * np.sqrt(KP_BASE)
KP_FOOT = np.array([80.0, 100.0, 80.0])
KD_FOOT = 2 * np.sqrt(KP_FOOT)


def mass_matrix():
    M = np.zeros((model.nv, model.nv))
    for j in range(model.nv):
        unit = np.zeros(model.nv)
        unit[j] = 1.0
        col = np.zeros(model.nv)
        mujoco.mj_mulM(model, data, col, unit)
        M[:, j] = col
    return M


def foot_pose_jac(d):
    jp = np.zeros((3, model.nv))
    jr = np.zeros((3, model.nv))
    mujoco.mj_jacSite(model, d, jp, jr, site_id)
    pose = np.array([d.site_xpos[site_id, 0], d.site_xpos[site_id, 2],
                     np.sum(d.qpos[2:5])])
    return pose, np.vstack((jp[0], jp[2], -jr[1]))


def jdot_qdot(q, v, eps=1e-5):
    scratch = mujoco.MjData(model)
    Js = []
    for sign in (1, -1):
        scratch.qpos[:] = q + sign * eps * v
        scratch.qvel[:] = v
        mujoco.mj_kinematics(model, scratch)
        mujoco.mj_comPos(model, scratch)
        Js.append(foot_pose_jac(scratch)[1])
    return ((Js[0] - Js[1]) / (2 * eps)) @ v


def measured_contact_wrench(d):
    """Sum foot-ground contact forces/moments about foot_center.

    MuJoCo contact.frame is a 3x3 matrix whose ROWS are contact-frame axes
    expressed in world coordinates. mj_contactForce returns [normal, tangent1,
    tangent2, torsion, rolling1, rolling2] in this contact frame, with the force
    acting on geom2. Sign is reversed when the foot is geom1.

    Returns w=[Fx,Fz,mc], with mc=-world_y_torque so J_foot.T @ w
    uses the same convention as the WBC-QP.
    """
    force_total = np.zeros(3)
    torque_total = np.zeros(3)
    other_contacts = 0
    used = 0
    origin = d.site_xpos[site_id].copy()
    for i in range(d.ncon):
        contact = d.contact[i]
        g1, g2 = int(contact.geom1), int(contact.geom2)
        if {g1, g2} != {ground_geom_id, foot_geom_id}:
            other_contacts += 1
            continue
        cf = np.zeros(6)
        mujoco.mj_contactForce(model, d, i, cf)
        R = np.asarray(contact.frame).reshape(3, 3)
        sign = 1.0 if g2 == foot_geom_id else -1.0
        f_world = sign * (R.T @ cf[:3])
        intrinsic_torque = sign * (R.T @ cf[3:6])
        r = np.asarray(contact.pos) - origin
        force_total += f_world
        torque_total += np.cross(r, f_world) + intrinsic_torque
        used += 1
    return np.array([force_total[0], force_total[2], -torque_total[1]]), used, other_contacts


# Same WBC-QP as experiment 10.3.
a = cp.Variable(5)
w = cp.Variable(3)
u = cp.Variable(2)
pM = cp.Parameter((5, 5))
ph = cp.Parameter(5)
pJ = cp.Parameter((3, 5))
pbias = cp.Parameter(3)
pfoot_acc = cp.Parameter(3)
pbase_acc = cp.Parameter(3)
pQpass = cp.Parameter(5)
S_T = np.zeros((5, 2))
S_T[3, 0] = S_T[4, 1] = 1.0
constraints = [
    pM @ a + ph == S_T @ u + pQpass + pJ.T @ w,
    pJ @ a + pbias == pfoot_acc,
    w[1] >= 0,
    w[0] <= MU * w[1], -w[0] <= MU * w[1],
    w[2] <= HALF_LENGTH * w[1], -w[2] <= HALF_LENGTH * w[1],
    u <= TAU_LIMIT, u >= -TAU_LIMIT,
]
objective = cp.Minimize(
    cp.sum_squares(cp.multiply([100.0, 100.0, 40.0], a[:3] - pbase_acc))
    + 1e-5 * cp.sum_squares(u)
    + 1e-7 * cp.sum_squares(w)
)
prob = cp.Problem(objective, constraints)

mujoco.mj_resetDataKeyframe(model, data, key_id)
data.ctrl[:] = 0.0
mujoco.mj_forward(model, data)
qref = data.qpos[:3].copy()
foot_ref, _ = foot_pose_jac(data)
steps_per_control = round(DT_CONTROL / model.opt.timestep)
assert abs(steps_per_control * model.opt.timestep - DT_CONTROL) < 1e-10
n_steps = round(T_FINAL / model.opt.timestep)

log = {k: [] for k in ('t', 'qp_fx', 'qp_fz', 'qp_mc', 'mj_fx', 'mj_fz',
                         'mj_mc', 'qp_cop', 'mj_cop', 'n_foot_contacts',
                         'n_other_contacts', 'generalized_force_error')}
fail = None

for step in range(n_steps):
    if step % steps_per_control == 0:
        # At the START of the control interval, solve for planned wrench/torques.
        mujoco.mj_forward(model, data)
        M = mass_matrix()
        p, J = foot_pose_jac(data)
        v = data.qvel.copy()
        base_acc_des = KP_BASE * (qref - data.qpos[:3]) - KD_BASE * v[:3]
        foot_acc_des = KP_FOOT * (foot_ref - p) - KD_FOOT * (J @ v)
        pM.value, ph.value = M, data.qfrc_bias.copy()
        pJ.value, pbias.value = J, jdot_qdot(data.qpos.copy(), v)
        pfoot_acc.value, pbase_acc.value = foot_acc_des, base_acc_des
        pQpass.value = data.qfrc_passive.copy()
        try:
            prob.solve(solver=cp.OSQP, warm_start=True, verbose=False,
                       eps_abs=1e-5, eps_rel=1e-5, max_iter=10000)
        except cp.error.SolverError as exc:
            fail = f'QP solver error at t={data.time:.3f}: {exc}'
            break
        if prob.status not in ('optimal', 'optimal_inaccurate') or w.value is None:
            fail = f'QP status={prob.status} at t={data.time:.3f}'
            break
        planned = np.asarray(w.value, dtype=float).reshape(3).copy()
        data.ctrl[:] = np.asarray(u.value, dtype=float).reshape(2)

    mujoco.mj_step(model, data)
    if not np.all(np.isfinite(data.qpos)):
        fail = f'Nonfinite state at t={data.time:.3f}'
        break

    if (step + 1) % steps_per_control == 0:
        # At the END of the interval, read actual contact wrench under the held
        # torque. mj_forward synchronizes contact force arrays with current q,v.
        mujoco.mj_forward(model, data)
        actual, n_foot, n_other = measured_contact_wrench(data)
        # Compare contact force mapping to generalized constraint force. This
        # check is exact only if all constraint forces come from foot-ground.
        _, J_now = foot_pose_jac(data)
        qforce_err = np.linalg.norm(J_now.T @ actual - data.qfrc_constraint)
        log['t'].append(float(data.time))
        for k, val in zip(('qp_fx', 'qp_fz', 'qp_mc'), planned):
            log[k].append(float(val))
        for k, val in zip(('mj_fx', 'mj_fz', 'mj_mc'), actual):
            log[k].append(float(val))
        log['qp_cop'].append(float(-planned[2] / planned[1]) if planned[1] > 1e-6 else np.nan)
        log['mj_cop'].append(float(-actual[2] / actual[1]) if actual[1] > 1e-6 else np.nan)
        log['n_foot_contacts'].append(n_foot)
        log['n_other_contacts'].append(n_other)
        log['generalized_force_error'].append(float(qforce_err))

print('MuJoCo', mujoco.__version__, '| controller', 1 / DT_CONTROL, 'Hz')
print('QP at interval START; MuJoCo wrench at interval END (10 ms later).')
if fail:
    print('STOP:', fail)
else:
    print('Completed:', round(data.time, 3), 's')
print('Final qpos:', np.round(data.qpos, 5))
print('Final qvel:', np.round(data.qvel, 5))

if log['t']:
    for component in ('fx', 'fz', 'mc'):
        qp = np.asarray(log['qp_' + component])
        mj = np.asarray(log['mj_' + component])
        print(f'{component}: final QP={qp[-1]:+.5f}, MJ={mj[-1]:+.5f}, '
              f'final diff={mj[-1] - qp[-1]:+.5f}, RMSE={np.sqrt(np.mean((mj-qp)**2)):.5f}')
    print('Max |J.T @ wrench_MJ - qfrc_constraint|:',
          np.max(log['generalized_force_error']))
    print('Other-geometry contact intervals:',
          np.count_nonzero(np.asarray(log['n_other_contacts']) > 0))
    print('No foot-ground contact intervals:',
          np.count_nonzero(np.asarray(log['n_foot_contacts']) == 0))

    t = np.asarray(log['t'])
    fig, axs = plt.subplots(4, 1, figsize=(11, 11), sharex=True)
    for ax, name, ylabel in zip(axs[:3], ('fx', 'fz', 'mc'),
                                ('Fx (N)', 'Fz (N)', 'mc (N m)')):
        ax.plot(t, log['qp_' + name], label='QP planned', lw=1.6)
        ax.plot(t, log['mj_' + name], label='MuJoCo actual', lw=1.2)
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.4)
        ax.legend()
    axs[3].plot(t, log['qp_cop'], label='QP CoP')
    axs[3].plot(t, log['mj_cop'], label='MuJoCo CoP')
    axs[3].axhline(HALF_LENGTH, color='gray', ls='--', lw=1)
    axs[3].axhline(-HALF_LENGTH, color='gray', ls='--', lw=1)
    axs[3].set_ylabel('CoP (m)')
    axs[3].set_xlabel('Time (s)')
    axs[3].grid(True, alpha=0.4)
    axs[3].legend()
    fig.suptitle('Day 10.4: QP wrench vs MuJoCo contact wrench (10 ms sample offset)')
    fig.tight_layout()
    out = ROOT / 'day10_4_contact_comparison.png'
    fig.savefig(out, dpi=160)
    print('Saved:', out)
    plt.show()
