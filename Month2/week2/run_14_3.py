"""Day 14.3: Search for a passive Compass Gait periodic orbit.
Run beside run_14_1.py. The numerical search may find no admissible orbit;
that is a valid result, not evidence that passive walking is impossible.\

Failed attempt: not able to find a physically admissible periodic orbit in the search, maybe due to the choice of angle of the platform. 
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import root
from run_14_1 import ode, foot_gap, swing_jacobian, impact, gamma

MU = 0.8
MAX_TIME = 2.0

def unpack(z):
    a, va, vb = z
    return np.array([a, 2*gamma-a, va, vb], dtype=float)

def pack(x):
    return x[[0, 2, 3]]

def gap_rate(x):
    return np.array([np.tan(gamma), 1.0]) @ swing_jacobian(x[:2]) @ x[2:]

def one_step(z, diagnostics=False):
    x0 = unpack(z)
    if not (-1.2 < x0[0] < gamma-0.005):
        raise ValueError('Postimpact stance angle outside forward-step region')
    if gap_rate(x0) <= 1e-7:
        raise ValueError('Old foot does not lift off')

    def event(t, x):
        # At t=0, both feet touch the ground; find the NEXT descending crossing.
        if t < 1e-7:
            return max(1e-9, foot_gap(x[:2]))
        return foot_gap(x[:2])
    event.terminal = True
    event.direction = -1
    sol = solve_ivp(ode, (0, MAX_TIME), x0, events=event,
                    rtol=2e-9, atol=1e-11, max_step=0.01)
    if not sol.t_events[0].size:
        raise ValueError('No next touchdown')
    t = float(sol.t_events[0][0])
    if t < 1e-4:
        raise ValueError('Near-immediate touchdown')
    xm = sol.y_events[0][0]
    if gap_rate(xm) >= -1e-7:
        raise ValueError('Not an approaching touchdown')
    qnew, vnew, impulse, d = impact(xm[:2], xm[2:])
    xn = np.r_[qnew, vnew]
    if diagnostics:
        return pack(xn), dict(time=t, impulse=impulse, liftoff=gap_rate(xn),
                              loss=d['kinetic_before']-d['kinetic_after'],
                              pre=xm, post=xn)
    return pack(xn)

def residual(z):
    try:
        return one_step(z)-z
    except (ValueError, np.linalg.LinAlgError, FloatingPointError):
        # Avoid treating failed trajectory integration as a valid fixed point.
        return np.array([10., 10., 10.])

def admissible(d):
    lam = d['impulse']
    normal = np.dot(lam, [np.sin(gamma), np.cos(gamma)])
    tangent = np.dot(lam, [np.cos(gamma), -np.sin(gamma)])
    return (normal > 0 and abs(tangent) <= MU*normal
            and d['liftoff'] > 1e-7 and d['loss'] >= -1e-7)

def main():
    rng = np.random.default_rng(14)
    guesses = [np.array([np.deg2rad(-15), 0.8, -0.5]),
               np.array([np.deg2rad(-20), 1.0, -1.0])]
    for _ in range(100):
        guesses.append(np.array([rng.uniform(-0.6, -0.04),
                                 rng.uniform(-0.2, 2.5),
                                 rng.uniform(-3, 2)]))
    roots = []
    for k, guess in enumerate(guesses, 1):
        result = root(residual, guess, method='hybr', options={'maxfev': 100})
        if np.linalg.norm(residual(result.x)) > 1e-6:
            continue
        if any(np.linalg.norm(result.x-r) < 1e-4 for r in roots):
            continue
        roots.append(result.x.copy())
        _, d = one_step(result.x, diagnostics=True)
        print(f'Candidate {len(roots)}: z={result.x}, time={d["time"]:.6f}, '
              f'impulse={d["impulse"]}, liftoff={d["liftoff"]:.6f}, '
              f'admissible={admissible(d)}')
        if admissible(d):
            zstar = result.x
            print('\nPHYSICALLY ADMISSIBLE PERIODIC ORBIT FOUND')
            print('q1, q2 [deg]:', np.rad2deg(unpack(zstar)[:2]))
            print('qdot1, qdot2 [rad/s]:', unpack(zstar)[2:])
            print('fixed-point residual:', np.linalg.norm(residual(zstar)))
            # Central finite differences for the Poincare Jacobian.
            eps = 1e-5
            jac = np.column_stack([(one_step(zstar+eps*np.eye(3)[j])-
                                    one_step(zstar-eps*np.eye(3)[j]))/(2*eps)
                                   for j in range(3)])
            eig = np.linalg.eigvals(jac)
            print('Poincare eigenvalues:', eig)
            print('spectral radius:', max(abs(eig)))
            np.save('day14_3_fixed_point.npy', zstar)
            return
    print(f'\nSearch finished: {len(roots)} mathematical candidate(s), '
          'no admissible periodic orbit found in this search.' )
    print('This does NOT prove nonexistence. Try other parameter values or search methods.')

if __name__ == '__main__':
    main()
