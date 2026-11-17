"""Day 14.2: Multi-step Compass Gait (point masses, passive walking).
Place this file next to run_14_1.py and run: python run_14_2.py
This is a deliberately minimal hybrid simulation, not a contact-complementarity solver.

Results: 
Failed simulation, invalid initial condition (not satisfying the absorbing boundary condition for the periodic motion)
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from run_14_1 import (ode, foot_gap, swing_jacobian, impact, energy, gamma)

MAX_STEPS = 12
MAX_STEP_DURATION = 3.0
MIN_STEP_DURATION = 1e-5
GAP_TOL = 1e-7
SPEED_TOL = 1e-7
MU = 0.8  # illustrative Coulomb friction coefficient at impact


def gap_rate(q, v):
    """Time derivative of vertical clearance y + x*tan(gamma)."""
    return np.array([np.tan(gamma), 1.0]) @ swing_jacobian(q) @ v


def integrate_one_step(t0, x0):
    """Find the next downward crossing of the swing-foot ground clearance."""
    # At the start of a step the old stance foot is exactly on the ground.
    # We only allow the event after it has lifted off; avoid a spurious t=t0 root.
    def touchdown(t, x):
        h = foot_gap(x[:2])
        if t - t0 < 1e-7 and abs(h) < 1e-9:
            return 1e-9
        return h
    touchdown.terminal = True
    touchdown.direction = -1
    return solve_ivp(ode, (t0, t0 + MAX_STEP_DURATION), x0,
                     events=touchdown, max_step=0.002,
                     rtol=1e-10, atol=1e-12)


def simulate_steps():
    x = np.array([np.deg2rad(15.), np.deg2rad(-15.), 0.8, -0.5])
    t = 0.0
    segments, impacts = [], []
    reason = f'Reached MAX_STEPS={MAX_STEPS}.'
    for k in range(MAX_STEPS):
        h0 = foot_gap(x[:2]); dh0 = gap_rate(x[:2], x[2:])
        if h0 < -GAP_TOL:
            reason = f'Step {k+1}: swing foot starts below the ground.'; break
        if abs(h0) < GAP_TOL and dh0 <= SPEED_TOL:
            reason = f'Step {k+1}: former stance foot does not lift off (gap rate {dh0:+.3g}).'; break
        sol = integrate_one_step(t, x)
        segments.append(sol)
        if not sol.t_events[0].size:
            reason = f'Step {k+1}: no touchdown within {MAX_STEP_DURATION:g} s.'; break
        tc = float(sol.t_events[0][0]); xm = sol.y_events[0][0]
        if tc-t < MIN_STEP_DURATION:
            reason = f'Step {k+1}: near-immediate touchdown; possible invalid contact mode.'; break
        if gap_rate(xm[:2], xm[2:]) >= -SPEED_TOL:
            reason = f'Step {k+1}: touchdown is not approaching the ground.'; break
        q_new, v_new, lam, d = impact(xm[:2], xm[2:])
        new_gap_rate = gap_rate(q_new, v_new)
        loss = d['kinetic_before'] - d['kinetic_after']
        friction_ratio = abs(lam[0]) / lam[1] if lam[1] > 0 else np.inf
        impacts.append(dict(step=k+1, time=tc, duration=tc-t,
                            pre=xm.copy(), post=np.r_[q_new,v_new],
                            impulse=lam.copy(), loss=loss,
                            liftoff_rate=new_gap_rate,
                            friction_ratio=friction_ratio))
        print(f'Step {k+1:2d}: t={tc:.6f}s, dt={tc-t:.4f}s, '
              f'Lambda=({lam[0]:+.3f},{lam[1]:+.3f}) Ns, '
              f'energy loss={loss:.3f}J, old-foot gap rate={new_gap_rate:+.4f}m/s')
        if lam[1] <= 0:
            reason = f'Step {k+1}: negative normal impulse; impact invalid.'; break
        if friction_ratio > MU:
            reason = f'Step {k+1}: friction ratio {friction_ratio:.3f} exceeds mu={MU}.'; break
        if new_gap_rate <= SPEED_TOL:
            reason = f'Step {k+1}: old stance foot fails to lift off.'; break
        if loss < -1e-8:
            reason = f'Step {k+1}: impact increased kinetic energy.'; break
        t, x = tc, np.r_[q_new, v_new]
    print('\nStop reason:', reason)
    return segments, impacts


def plot_results(segments, impacts):
    if not segments: return
    fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    for seg in segments:
        ax[0].plot(seg.t, np.rad2deg(seg.y[0]), color='C0')
        ax[0].plot(seg.t, np.rad2deg(seg.y[1]), color='C1')
        ax[1].plot(seg.t, seg.y[2], color='C0')
        ax[1].plot(seg.t, seg.y[3], color='C1')
        ax[2].plot(seg.t, [foot_gap(seg.y[:2,i]) for i in range(seg.t.size)], color='C2')
    for d in impacts:
        for a in ax: a.axvline(d['time'], ls=':', color='0.45', lw=1)
        # Plot the velocity jump (newly relabeled angles) at impact.
        ax[0].plot([d['time'], d['time']],
                   np.rad2deg([d['pre'][0],d['post'][0]]), color='0.6', ls='--', alpha=.5)
    ax[0].set_ylabel('Angles (deg)')
    ax[1].set_ylabel('Angular speeds (rad/s)')
    ax[2].set_ylabel('Swing-foot gap (m)')
    ax[2].set_xlabel('Time (s)')
    ax[2].axhline(0, color='k', ls='--', lw=.8)
    for a in ax: a.grid(alpha=.25)
    ax[0].set_title('Blue: stance q1, Orange: swing q2 (relabel at each impact)')
    fig.suptitle('Day 14.2: Compass Gait hybrid simulation')
    fig.tight_layout()
    fig.savefig('day14_2_multistep.png', dpi=160)
    print('Saved day14_2_multistep.png')


if __name__ == '__main__':
    segs, imps = simulate_steps()
    plot_results(segs, imps)
