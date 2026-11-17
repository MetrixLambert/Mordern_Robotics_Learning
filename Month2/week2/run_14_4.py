"""Day 14.4: Direct shooting for passive Compass Gait periodic orbits.

Place beside run_14_1.py. Requires numpy, scipy, matplotlib.
Unknowns: z=(post-impact q1, qdot1, qdot2) and step duration T.
Residual: post-impact periodicity (3) + swing-foot touchdown gap (1).
Unlike event-based root finding, this residual is smooth near valid trajectories.
"""
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import run_14_1 as model

OUT = Path(__file__).resolve().parent
MU = 0.8
# slope angles, degrees. The model uses absolute angles, so gamma enters only in contact geometry.
SLOPES_DEG = [3.0, 5.0, 7.0]
MAX_NFEV = 65


def unpack(z, gamma):
    a, da, db = z
    return np.array([a, 2 * gamma - a, da, db], dtype=float)


def gap(q, gamma):
    px, py = model.swing_position(q)
    return py + px * np.tan(gamma)


def gap_speed(x, gamma):
    return np.array([np.tan(gamma), 1.0]) @ model.swing_jacobian(x[:2]) @ x[2:]


def integrate(z, T, gamma, dense=False):
    x0 = unpack(z, gamma)
    return solve_ivp(model.ode, (0.0, float(T)), x0,
                     rtol=2e-8, atol=2e-10, max_step=0.035,
                     dense_output=dense)


def shoot(y, gamma, details=False):
    z, T = y[:3], y[3]
    sol = integrate(z, T, gamma, dense=details)
    xm = sol.y[:, -1]
    qn, vn, lam, d = model.impact(xm[:2], xm[2:])
    zn = np.r_[qn[0], vn]
    # scale each residual to balance angles, velocities, and foot gap
    r = np.r_[(zn - z) / np.array([0.3, 1.5, 2.0]), gap(xm[:2], gamma) / 0.3]
    if details:
        return r, sol, xm, np.r_[qn, vn], lam, d
    return r


def feasible(y, gamma, tol=1e-5):
    r, sol, xm, xn, lam, d = shoot(y, gamma, details=True)
    z, T = y[:3], y[3]
    normal = lam @ np.array([np.sin(gamma), np.cos(gamma)])
    tangent = lam @ np.array([np.cos(gamma), -np.sin(gamma)])
    times = np.linspace(0, T, 51)
    gaps = np.array([gap(sol.sol(t)[:2], gamma) for t in times])
    start_rate = gap_speed(unpack(z, gamma), gamma)
    end_rate = gap_speed(xm, gamma)
    after_rate = gap_speed(xn, gamma)
    valid = (np.linalg.norm(r) < tol and T > 0.12 and
             np.all(gaps[1:-1] > 1e-7) and
             start_rate > 1e-5 and end_rate < -1e-5 and
             after_rate > 1e-5 and normal > 1e-6 and
             abs(tangent) <= MU * normal and
             d['kinetic_before'] >= d['kinetic_after'] - 1e-7 and
             xn[0] < gamma - 1e-3 and
             model.swing_position(xm[:2])[0] > 0.02)
    diag = dict(residual=np.linalg.norm(r), T=T, normal=normal, tangent=tangent,
                start_gap_rate=start_rate, touchdown_gap_rate=end_rate,
                liftoff_rate=after_rate, min_interior_gap=np.min(gaps[1:-1]),
                energy_loss=d['kinetic_before']-d['kinetic_after'])
    return valid, diag


def search():
    rng = np.random.default_rng(1404)
    bounds = (np.array([-0.85, -0.5, -4.0, 0.15]),
              np.array([-0.015, 3.5, 4.0, 2.5]))
    # Deterministic seeds plus random guesses; can expand N_RANDOM for a broader search.
    seeds = [
        [-0.22, 0.8, -0.5, 0.6],
        [-0.30, 0.7, -0.4, 0.8],
        [-0.16, 0.7, -1.0, 0.6],
        [-0.12, 0.5, -0.5, 0.5],
        [-0.4, 1.2, 0.0, 0.8],
        [-0.2, 1.1, 0.5, 0.55],
    ]
    N_RANDOM = 2
    seeds.extend([[rng.uniform(-0.6, -0.05), rng.uniform(0.15, 2.0),
                   rng.uniform(-2, 2), rng.uniform(0.25, 1.4)]
                  for _ in range(N_RANDOM)])
    for deg in SLOPES_DEG:
        gamma = np.deg2rad(deg)
        model.gamma = gamma
        print(f'\n--- slope {deg:g} deg; {len(seeds)} starts ---', flush=True)
        best = (np.inf, None)
        for k, seed in enumerate(seeds, 1):
            fit = least_squares(shoot, seed, args=(gamma,), bounds=bounds,
                                max_nfev=MAX_NFEV, ftol=1e-9,
                                xtol=1e-9, gtol=1e-9)
            err = np.linalg.norm(fit.fun)
            if err < best[0]:
                best = (err, fit.x.copy())
            if err < 1e-5:
                ok, d = feasible(fit.x, gamma)
                if ok:
                    print(f'Found admissible orbit on start {k}, residual={err:.3e}')
                    return gamma, fit.x
            if k % 8 == 0:
                print(f'  {k:2d}/{len(seeds)}, best scaled residual={best[0]:.3e}', flush=True)
        print(f'No admissible orbit at {deg:g} deg; best residual={best[0]:.3e}', flush=True)
    return None


def jacobian_map(z, gamma, Tguess):
    """Event-based Poincare map only used near an already validated periodic orbit."""
    def P(z0):
        x0 = unpack(z0, gamma)
        def event(t, x):
            # no initial t=0 event: physically lift-off has positive gap rate
            if t < 1e-8:
                return max(1e-10, gap(x[:2], gamma))
            return gap(x[:2], gamma)
        event.terminal, event.direction = True, -1
        sol = solve_ivp(model.ode, (0, max(2*Tguess, 1.0)), x0,
                        events=event, rtol=1e-10, atol=1e-12, max_step=0.01)
        if not sol.t_events[0].size:
            raise RuntimeError('No touchdown in perturbed Poincare trajectory')
        xm = sol.y_events[0][0]
        qn, vn, _, _ = model.impact(xm[:2], xm[2:])
        return np.r_[qn[0], vn]
    h = 2e-5
    eye = np.eye(3)
    A = np.column_stack([(P(z+h*eye[j])-P(z-h*eye[j]))/(2*h) for j in range(3)])
    return A, np.linalg.eigvals(A)


def plot_orbit(gamma, y):
    r, sol, xm, xn, lam, d = shoot(y, gamma, details=True)
    ts = np.linspace(0, y[3], 250)
    xs = sol.sol(ts)
    fig, axs = plt.subplots(3, 1, figsize=(9, 8), sharex=True)
    axs[0].plot(ts, np.rad2deg(xs[0]), label='stance q1')
    axs[0].plot(ts, np.rad2deg(xs[1]), label='swing q2')
    axs[0].set_ylabel('Angle (deg)'); axs[0].legend()
    axs[1].plot(ts, xs[2], label='stance qdot1')
    axs[1].plot(ts, xs[3], label='swing qdot2')
    axs[1].set_ylabel('Angular speed (rad/s)'); axs[1].legend()
    axs[2].plot(ts, [gap(xs[:2, i], gamma) for i in range(xs.shape[1])])
    axs[2].axhline(0, color='k', lw=0.8, ls='--')
    axs[2].set_ylabel('Swing-foot gap (m)'); axs[2].set_xlabel('Time (s)')
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.axvline(y[3], color='r', ls=':', alpha=0.6)
    fig.suptitle(f'Day 14.4: periodic passive gait, slope {np.rad2deg(gamma):g} deg')
    fig.tight_layout()
    path = OUT / 'day14_4_periodic_orbit.png'
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print('Saved:', path)


def main():
    result = search()
    if result is None:
        print('\nNo admissible periodic orbit found in this bounded search.')
        print('This does not prove nonexistence. Expand seeds, slopes, or change mass ratios.')
        return
    gamma, y = result
    ok, d = feasible(y, gamma)
    print('\nPHYSICALLY ADMISSIBLE PERIODIC ORBIT (within tested impact assumptions)')
    print('slope (deg):', np.rad2deg(gamma))
    print('postimpact z=[q1, qdot1, qdot2]:', y[:3])
    print('postimpact q1, q2 (deg):', np.rad2deg(unpack(y[:3], gamma)[:2]))
    print('step time (s):', y[3])
    for key, value in d.items():
        print(f'{key}: {value:+.9g}')
    try:
        A, eig = jacobian_map(y[:3], gamma, y[3])
        print('Poincare eigenvalues:', eig)
        print('spectral radius:', np.max(np.abs(eig)))
    except RuntimeError as e:
        print('Stability analysis skipped:', e)
    np.savez(OUT / 'day14_4_solution.npz', z=y[:3], T=y[3], gamma=gamma)
    plot_orbit(gamma, y)


if __name__ == '__main__':
    main()
