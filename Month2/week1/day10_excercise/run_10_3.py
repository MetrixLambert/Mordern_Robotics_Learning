"""Day 10 / 10.3: closed-loop planar stance WBC-QP (MuJoCo 3.17 compatible).

Put next to planar_leg.xml. Install: pip install mujoco cvxpy osqp matplotlib
Run: python run_10_3.py

This is a *stance-mode* controller, not a contact-implicit controller. It assumes
that the foot remains in contact; stop interpreting results as stance WBC if it lifts.
"""
from pathlib import Path
import numpy as np
import mujoco
import cvxpy as cp
import matplotlib.pyplot as plt

MODEL_PATH = Path(__file__).with_name('planar_leg.xml')
model = mujoco.MjModel.from_xml_path(str(MODEL_PATH))
data = mujoco.MjData(model)
assert (model.nq, model.nv, model.nu) == (5, 5, 2)
site_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, 'foot_center')
key_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_KEY, 'initial')
assert site_id >= 0 and key_id >= 0

# Contact and control parameters. The physical half-length of the foot is 0.12 m.
MU, FOOT_HALF_LENGTH = 0.8, 0.12
TAU_LIMIT = 200.0
DT_CONTROL = 0.01     # 100 Hz controller; MuJoCo integrates at 500 Hz
T_FINAL = 3.0
KP_BASE = np.array([35.0, 55.0, 55.0])
KD_BASE = 2.0 * np.sqrt(KP_BASE)
KP_FOOT = np.array([80.0, 100.0, 80.0])
KD_FOOT = 2.0 * np.sqrt(KP_FOOT)


def mass_matrix():
    """Use mj_mulM to avoid depending on data.qM (not exposed in some bindings)."""
    M = np.zeros((model.nv, model.nv))
    for j in range(model.nv):
        e = np.zeros(model.nv)
        e[j] = 1.0
        col = np.zeros(model.nv)
        mujoco.mj_mulM(model, data, col, e)
        M[:, j] = col
    return M


def foot_pose_jac(d):
    jp = np.zeros((3, model.nv))
    jr = np.zeros((3, model.nv))
    mujoco.mj_jacSite(model, d, jp, jr, site_id)
    # In XML the hinges rotate around -Y; phi = theta + hip + knee.
    p = np.array([d.site_xpos[site_id, 0], d.site_xpos[site_id, 2],
                  np.sum(d.qpos[2:5])])
    J = np.vstack((jp[0], jp[2], -jr[1]))
    return p, J


def jdot_qdot(q, v, eps=1e-5):
    """Central directional difference: Jdot(q,v) @ v.

    Here all generalized coordinates are slide/hinge so q +/- eps*v is valid.
    Scratch MjData keeps the live simulation untouched.
    """
    scratch = mujoco.MjData(model)
    Jpm = []
    for sign in (1, -1):
        scratch.qpos[:] = q + sign * eps * v
        scratch.qvel[:] = v
        mujoco.mj_kinematics(model, scratch)
        mujoco.mj_comPos(model, scratch)
        Jpm.append(foot_pose_jac(scratch)[1])
    return ((Jpm[0] - Jpm[1]) / (2.0 * eps)) @ v


# QP compiled once, numeric parameters updated every control step.
a = cp.Variable(5)
w = cp.Variable(3)       # [fx, fz, mc] about foot_center
u = cp.Variable(2)       # joint motor torques
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
    w[2] <= FOOT_HALF_LENGTH * w[1],
    -w[2] <= FOOT_HALF_LENGTH * w[1],
    u <= TAU_LIMIT, u >= -TAU_LIMIT,
]
objective = cp.Minimize(
    cp.sum_squares(cp.multiply(np.array([100.0, 100.0, 40.0]), a[:3] - pbase_acc))
    + 1e-5 * cp.sum_squares(u)
    + 1e-7 * cp.sum_squares(w)
)
prob = cp.Problem(objective, constraints)

mujoco.mj_resetDataKeyframe(model, data, key_id)
data.ctrl[:] = 0
mujoco.mj_forward(model, data)
qref = data.qpos[:3].copy()
foot_ref, _ = foot_pose_jac(data)
print('MuJoCo', mujoco.__version__, '| initial qbase', qref, '| foot', foot_ref)

log = {k: [] for k in ('t', 'x', 'z', 'theta', 'fx_qp', 'fz_qp', 'mc_qp',
                         'tau_hip', 'tau_knee', 'ncon', 'cop_qp', 'foot_error')}
steps_per_control = max(1, round(DT_CONTROL / model.opt.timestep))
n_steps = int(round(T_FINAL / model.opt.timestep))
fail = None

for step in range(n_steps):
    if step % steps_per_control == 0:
        # Recompute model quantities at the CURRENT state (before stepping).
        mujoco.mj_forward(model, data)
        M = mass_matrix()
        p, J = foot_pose_jac(data)
        v = data.qvel.copy()
        bias = jdot_qdot(data.qpos.copy(), v)
        base_error = qref - data.qpos[:3]
        base_acc_des = KP_BASE * base_error - KD_BASE * v[:3]
        foot_vel = J @ v
        foot_acc_des = KP_FOOT * (foot_ref - p) - KD_FOOT * foot_vel

        pM.value, ph.value = M, data.qfrc_bias.copy()
        pJ.value, pbias.value = J, bias
        pfoot_acc.value, pbase_acc.value = foot_acc_des, base_acc_des
        pQpass.value = data.qfrc_passive.copy()
        try:
            prob.solve(solver=cp.OSQP, warm_start=True, verbose=False,
                       eps_abs=1e-5, eps_rel=1e-5, max_iter=10000)
        except cp.error.SolverError as exc:
            fail = f'QP solver error at t={data.time:.3f}: {exc}'
            break
        if prob.status not in ('optimal', 'optimal_inaccurate') or u.value is None:
            fail = f'QP {prob.status} at t={data.time:.3f}'
            break
        data.ctrl[:] = np.asarray(u.value).ravel()
        wc = np.asarray(w.value).ravel()
        log['t'].append(float(data.time))
        log['x'].append(float(data.qpos[0]))
        log['z'].append(float(data.qpos[1]))
        log['theta'].append(float(data.qpos[2]))
        log['fx_qp'].append(float(wc[0]))
        log['fz_qp'].append(float(wc[1]))
        log['mc_qp'].append(float(wc[2]))
        log['tau_hip'].append(float(data.ctrl[0]))
        log['tau_knee'].append(float(data.ctrl[1]))
        log['ncon'].append(int(data.ncon))
        log['cop_qp'].append(float(-wc[2] / wc[1]) if wc[1] > 1e-6 else np.nan)
        log['foot_error'].append(float(np.linalg.norm(foot_ref - p)))
    mujoco.mj_step(model, data)
    if not np.all(np.isfinite(data.qpos)):
        fail = f'Nonfinite state at t={data.time:.3f}'
        break

if fail:
    print('STOP:', fail)
else:
    print('Completed', round(data.time, 3), 'seconds')
print('final qpos:', np.round(data.qpos, 5))
print('final qvel:', np.round(data.qvel, 5))
print('final contacts:', data.ncon)
if log['t']:
    print('max |theta| [rad]:', np.max(np.abs(log['theta'])))
    print('max foot pose error:', np.max(log['foot_error']))
    print('min QP normal force [N]:', np.min(log['fz_qp']))
    print('fraction of controller samples with no contact:',
          np.mean(np.asarray(log['ncon']) == 0))

fig, axs = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
t = np.asarray(log['t'])
if t.size:
    axs[0].plot(t, log['x'], label='base x')
    axs[0].plot(t, log['z'], label='base z displacement')
    axs[0].plot(t, log['theta'], label='base theta (rad)')
    axs[0].legend(); axs[0].grid(True)
    axs[1].plot(t, log['tau_hip'], label='hip torque')
    axs[1].plot(t, log['tau_knee'], label='knee torque')
    axs[1].legend(); axs[1].grid(True)
    axs[2].plot(t, log['fz_qp'], label='QP fz (N)')
    axs[2].plot(t, log['cop_qp'], label='QP CoP (m)')
    axs[2].plot(t, log['ncon'], label='contact count')
    axs[2].legend(); axs[2].grid(True)
axs[2].set_xlabel('time (s)')
fig.tight_layout()
output = Path(__file__).with_name('day10_3_results.png')
fig.savefig(output, dpi=150)
print('Saved plot:', output)
plt.show()
