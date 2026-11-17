"""Day 11 / 11.3: planar single-leg touchdown and impact impulse.

Place this file next to the Day 10 planar_leg.xml, or set XML_PATH below.
Tested for syntax; intended for MuJoCo 3.x (including 3.17.0).
No WBC controller: motors are set to zero.
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import mujoco

HERE = Path(__file__).resolve().parent
XML_PATH = HERE / 'planar_leg.xml'
if not XML_PATH.exists():
    raise FileNotFoundError('Place planar_leg.xml beside this script, or edit XML_PATH.')

model = mujoco.MjModel.from_xml_path(str(XML_PATH))
data = mujoco.MjData(model)
assert (model.nq, model.nv) == (5, 5), 'Expected the Day 10 planar robot'
site_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_SITE, 'foot_center')
ground_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_GEOM, 'ground')
key_id = mujoco.mj_name2id(model, mujoco.mjtObj.mjOBJ_KEY, 'initial')
assert min(site_id, ground_id, key_id) >= 0

# Start 0.20 m above the original standing pose, with zero motor torques.
HEIGHT = 0.20
T_END = 0.55
POST_WINDOW = 0.030   # 30 ms: inspect finite-time collision, not an instantaneous reset

mujoco.mj_resetDataKeyframe(model, data, key_id)
data.qpos[1] += HEIGHT
data.qvel[:] = 0
if model.nu:
    data.ctrl[:] = 0
mujoco.mj_forward(model, data)


def full_mass_matrix():
    """MuJoCo 3.17.0-compatible extraction without data.qM."""
    M = np.zeros((model.nv, model.nv))
    for k in range(model.nv):
        e = np.zeros(model.nv)
        e[k] = 1
        out = np.zeros(model.nv)
        mujoco.mj_mulM(model, data, out, e)
        M[:, k] = out
    return M


def normal_jacobian():
    jp = np.zeros((3, model.nv))
    jr = np.zeros((3, model.nv))
    mujoco.mj_jacSite(model, data, jp, jr, site_id)
    return jp[2, :].copy()


def ground_forces():
    """World force exerted by ground on ALL robot geometries.

    Contact.frame rows are contact axes in world coordinates; mj_contactForce
    gives force acting on geom2. Summing all contacts avoids the Day 10
    shank-vs-foot omission. In the intended model, ground is static.
    """
    total = np.zeros(3)
    n_ground = 0
    names = set()
    for k in range(data.ncon):
        c = data.contact[k]
        g1, g2 = int(c.geom1), int(c.geom2)
        if ground_id not in (g1, g2):
            continue
        cf = np.zeros(6)
        mujoco.mj_contactForce(model, data, k, cf)
        R = np.asarray(c.frame).reshape(3, 3)
        sign = 1.0 if g2 != ground_id else -1.0
        total += sign * (R.T @ cf[:3])
        other = g2 if g1 == ground_id else g1
        names.add(mujoco.mj_id2name(model, mujoco.mjtObj.mjOBJ_GEOM, other))
        n_ground += 1
    return total, n_ground, sorted(names)


# Log state at step boundaries, then use forces evaluated at those boundaries.
# Force integration is approximate because contact forces change within a step.
records = []
first_hit = None
pre = None
while data.time < T_END - 1e-12:
    # This is the state immediately before the next integration step.
    q_pre = data.qpos.copy()
    v_pre = data.qvel.copy()
    t_pre = float(data.time)
    mujoco.mj_step(model, data)
    F, n, names = ground_forces()
    records.append((float(data.time), data.qpos[1], data.qvel[1],
                    data.site_xpos[site_id, 2], F[2], n, *data.qvel.copy()))
    if first_hit is None and n > 0:
        first_hit = len(records) - 1
        pre = (t_pre, q_pre, v_pre)
        print(f'First contact: t={data.time:.6f} s; geoms={names}')

arr = np.asarray(records)
if first_hit is None:
    raise RuntimeError('No contact detected. Increase T_END or check the XML geometry.')

t_pre, q_minus, v_minus = pre
# Evaluate mass matrix and J at the PRE-impact configuration, in free flight.
# Use a separate MjData so the post-impact simulation state is unchanged.
scratch = mujoco.MjData(model)
scratch.qpos[:] = q_minus
scratch.qvel[:] = v_minus
mujoco.mj_forward(model, scratch)
M = np.zeros((model.nv, model.nv))
for k in range(model.nv):
    e = np.zeros(model.nv)
    e[k] = 1
    out = np.zeros(model.nv)
    mujoco.mj_mulM(model, scratch, out, e)
    M[:, k] = out
jp = np.zeros((3, model.nv))
jr = np.zeros((3, model.nv))
mujoco.mj_jacSite(model, scratch, jp, jr, site_id)
Jn = jp[2, :].copy()

# Ideal perfectly inelastic, frictionless, ONE-normal-constraint prediction.
Minv_Jt = np.linalg.solve(M, Jn)
normal_effective_inverse_mass = float(Jn @ Minv_Jt)
v_n_minus = float(Jn @ v_minus)
Lambda_ideal = max(0.0, -v_n_minus / normal_effective_inverse_mass)
v_plus_ideal = v_minus + Minv_Jt * Lambda_ideal

# Measured force impulse in a short, fixed window after first contact.
# This is not expected to match Lambda_ideal: finite contact duration, multiple
# contact points, gravity, passive forces, friction, and possible rebound.
t_hit = arr[first_hit, 0]
end_idx = int(np.searchsorted(arr[:, 0], t_hit + POST_WINDOW, side='right'))
window = arr[first_hit:end_idx]
Lambda_measured = float(np.sum(window[:, 4]) * model.opt.timestep)

print('\n--- Pre-impact / ideal single-normal collision ---')
print('Pre-impact qvel:', np.round(v_minus, 5))
print(f'Pre-impact foot normal speed Jn @ v-: {v_n_minus:+.5f} m/s')
print(f'Ideal normal impulse: {Lambda_ideal:.5f} N*s')
print('Ideal post-impact qvel:', np.round(v_plus_ideal, 5))
print(f'Ideal post-impact normal speed: {Jn @ v_plus_ideal:+.3e} m/s')
print('\n--- MuJoCo finite-time contact ---')
print(f'Window: {window[0, 0]:.4f} to {window[-1, 0]:.4f} s ({len(window)} samples)')
print(f'Integrated world vertical ground force: {Lambda_measured:.5f} N*s')
print(f'Max recorded vertical contact force: {np.max(window[:, 4]):.3f} N')
print(f'Foot-center height at first detected contact: {arr[first_hit, 3]:.5f} m')
print('qvel at window end:', np.round(arr[end_idx - 1, 6:11], 5))
print('NOTE: the measured 30 ms impulse includes gravity compensation and')
print('      potentially multiple contacts; it is NOT a direct equality test')
print('      against the ideal instantaneous single-normal impulse.')

fig, axs = plt.subplots(3, 1, figsize=(9, 9), sharex=True)
t = arr[:, 0]
axs[0].plot(t, arr[:, 1], label='base z')
axs[0].plot(t, arr[:, 3], label='foot-center z')
axs[0].set_ylabel('Height [m]')
axs[0].legend()
axs[1].plot(t, arr[:, 2], label='base vertical velocity')
axs[1].axhline(0, color='gray', lw=0.8)
axs[1].set_ylabel('Velocity [m/s]')
axs[1].legend()
axs[2].plot(t, arr[:, 4], label='total ground Fz')
axs[2].set_ylabel('Force [N]')
axs[2].set_xlabel('Time [s]')
axs[2].legend()
for ax in axs:
    ax.axvline(t_hit, color='tab:red', ls='--', lw=1, label='first contact')
    ax.grid(alpha=0.25)
fig.suptitle('Day 11: Single-leg touchdown (no control)')
fig.tight_layout()
out = HERE / 'day11_3_touchdown.png'
fig.savefig(out, dpi=160)
print(f'Plot saved: {out}')
plt.show()
