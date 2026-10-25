"""ROBO_Learn Day 15 | Dynamics validation, energy and integration.

Run: python day15_dynamics_validation_exercises.py
Requires: numpy; matplotlib optional (for figures).

Complete TODO 1-5. A self-contained 2R reference model is included; no imports
from Day 14 are required. Coordinates: XY vertical plane, Y up, q1 from +X.
"""
import numpy as np


class Planar2R:
    def __init__(self):
        self.L1, self.L2 = 2.0, 1.0
        self.m1, self.m2 = 2.0, 1.0
        self.I1, self.I2 = 0.25, 0.10
        self.gravity = 9.81


def mass_matrix(robot, q):
    """Working Day 14 analytic reference."""
    L1, r1, r2 = robot.L1, robot.L1 / 2, robot.L2 / 2
    m1, m2, I1, I2 = robot.m1, robot.m2, robot.I1, robot.I2
    a = I1 + I2 + m1*r1**2 + m2*(L1**2 + r2**2 + 2*L1*r2*np.cos(q[1]))
    b = I2 + m2*(r2**2 + L1*r2*np.cos(q[1]))
    d = I2 + m2*r2**2
    return np.array([[a, b], [b, d]])


def coriolis_matrix(robot, q, qd):
    h = robot.m2 * robot.L1 * robot.L2/2 * np.sin(q[1])
    return np.array([[-h*qd[1], -h*(qd[0]+qd[1])], [h*qd[0], 0.]])


def potential_energy(robot, q):
    g = robot.gravity
    return g*(robot.m1*robot.L1/2*np.sin(q[0]) +
              robot.m2*(robot.L1*np.sin(q[0]) + robot.L2/2*np.sin(q[0]+q[1])))


def gravity_vector(robot, q):
    g = robot.gravity
    return g*np.array([
        (robot.m1*robot.L1/2 + robot.m2*robot.L1)*np.cos(q[0]) +
        robot.m2*robot.L2/2*np.cos(q[0]+q[1]),
        robot.m2*robot.L2/2*np.cos(q[0]+q[1]),
    ])


def forward_dynamics(robot, q, qd, tau):
    return np.linalg.solve(mass_matrix(robot, q),
                           tau - coriolis_matrix(robot, q, qd) @ qd - gravity_vector(robot, q))


# ========== Fill these functions (leave the tests unchanged) ==========

def total_energy(robot, q, qd):
    """TODO 1: return K+U, where K=0.5*qd.T@M(q)@qd."""
    M = mass_matrix(robot, q) 
    K = 0.5 * qd.T @ M @ qd 
    U = potential_energy(robot, q) 
    return K + U 


def semi_implicit_euler_step(robot, q, qd, tau, dt):
    """TODO 2: qdd=f(q,qd,tau); qd_next=qd+dt*qdd; q_next=q+dt*qd_next.

    Return (q_next, qd_next). Never mutate input q or qd.
    """
    qdd = forward_dynamics(robot, q, qd, tau) 
    qd_next = qd + dt * qdd 
    q_next = q + dt * qd_next 
    return q_next, qd_next 


def simulate(robot, q0, qd0, torque_fn, dt, duration):
    """TODO 3: integrate with semi_implicit_euler_step.

    torque_fn(t,q,qd) -> (2,) torque. Evaluate at the BEGINNING of each step.
    n_steps = round(duration/dt); assert it is consistent with duration.
    Return t:(N+1,), q:(N+1,2), qd:(N+1,2), tau:(N,2).
    Record initial state and every subsequent state; do not mutate q0/qd0.
    """
    t = np.linspace(0, duration, int(duration/dt) + 1) 
    n_steps = len(t) - 1 
    qs = np.zeros((n_steps + 1, 2)) 
    qs[0] = q0 
    qds = np.zeros((n_steps + 1, 2))  
    qds[0] = qd0 
    taus = np.zeros((n_steps, 2)) 

    for k in range(n_steps): 
        tau = torque_fn(t[k], qs[k], qds[k])
        taus[k] = tau 
        qs[k + 1], qds[k + 1] = semi_implicit_euler_step(robot, qs[k], qds[k], tau, dt) 

    return t, qs, qds, taus 


def energy_balance_residual(robot, t, q, qd, tau):
    """TODO 4: energy change minus discrete input work.

    E[k]=total_energy(robot,q[k],qd[k]).
    Approximate work on each step by left-endpoint power:
        W[k+1]=W[k]+dt * tau[k]@qd[k].
    Return residual array R[k]=E[k]-E[0]-W[k], length len(t).

    For finite dt, R need not be exactly zero: integration AND quadrature
    errors contribute. This is a diagnostic, not an exact invariant.
    """
    dt = t[1] - t[0] 
    R = np.zeros(len(t)) 
    W = np.zeros(len(t)) 
    E = np.zeros(len(t)) 
    W[0] = 0 
    E[0] = total_energy(robot, q[0], qd[0]) 

    for k in range(len(t) -1): 
        W[k + 1] = W[k] + dt * tau[k] @ qd[k] 
        E[k + 1] = total_energy(robot, q[k + 1], qd[k + 1]) 
        R[k + 1] = E[k + 1] - E[0] - W[k + 1] 

    return R 


def compare_time_steps(robot, q0, qd0, dt_values, duration):
    """TODO 5: simulate FREE motion (tau=0) for each dt.

    Return a dict mapping float(dt) to max(abs(E-E[0])) across trajectory.
    Use total_energy and simulate. No plotting inside this function.
    """
    
    results = {} 
    for dt in dt_values: 
        t, q, qd, tau = simulate(robot, q0, qd0, lambda t,q,qd: np.zeros(2), dt, duration)
        E = np.array([total_energy(robot, qi, vi) for qi,vi in zip(q, qd)])
        results[dt] = np.max(np.abs(E - E[0]))
    return results 


# ========== Tests ==========

def run_tests():
    robot = Planar2R()
    rng = np.random.default_rng(1500)

    print('[1] Energy expression: 100 random states')
    for _ in range(100):
        q = rng.uniform(-2, 2, 2)
        qd = rng.uniform(-3, 3, 2)
        expected = 0.5*qd @ mass_matrix(robot, q) @ qd + potential_energy(robot, q)
        np.testing.assert_allclose(total_energy(robot, q, qd), expected, atol=1e-11)
    print('  PASS')

    print('[2] Semi-implicit Euler step: exact one-step identity')
    for _ in range(100):
        q = rng.uniform(-2, 2, 2)
        qd = rng.uniform(-2, 2, 2)
        tau = rng.uniform(-10, 10, 2)
        dt = 0.001
        q_old, qd_old = q.copy(), qd.copy()
        a = forward_dynamics(robot, q, qd, tau)
        v_ref = qd + dt*a
        q_ref = q + dt*v_ref
        q_new, v_new = semi_implicit_euler_step(robot, q, qd, tau, dt)
        np.testing.assert_allclose(q_new, q_ref, atol=1e-12)
        np.testing.assert_allclose(v_new, v_ref, atol=1e-12)
        np.testing.assert_array_equal(q, q_old)
        np.testing.assert_array_equal(qd, qd_old)
    print('  PASS')

    print('[3] Static gravity compensation and disturbed initial state')
    qstar = np.zeros(2)
    t, q, qd, tau = simulate(robot, qstar, np.zeros(2),
                             lambda t,q,qd: gravity_vector(robot, qstar),
                             dt=0.002, duration=1.0)
    np.testing.assert_allclose(q, 0, atol=1e-10)
    np.testing.assert_allclose(qd, 0, atol=1e-10)
    assert q.shape == (501, 2) and tau.shape == (500, 2)
    td, qd_hist, vd_hist, _ = simulate(robot, np.array([0.1, 0.0]), np.zeros(2),
                                       lambda t,q,qd: gravity_vector(robot, qstar),
                                       dt=0.002, duration=0.1)
    # At a disturbed configuration, constant gravity compensation is not equilibrium.
    assert np.linalg.norm(vd_hist[-1]) > 1e-4
    print('  PASS')

    print('[4] Free-motion energy residual and step-size comparison')
    initial_q = np.array([0.4, -0.6])
    initial_v = np.array([0.2, 0.1])
    t, q, qd, tau = simulate(robot, initial_q, initial_v,
                             lambda t,q,qd: np.zeros(2), dt=0.001, duration=0.5)
    residual = energy_balance_residual(robot, t, q, qd, tau)
    energies = np.array([total_energy(robot, qi, vi) for qi,vi in zip(q, qd)])
    np.testing.assert_allclose(residual, energies-energies[0], atol=1e-10)
    assert np.isfinite(residual).all()
    errors = compare_time_steps(robot, initial_q, initial_v, [0.004, 0.002, 0.001], 0.5)
    assert set(errors) == {0.004, 0.002, 0.001}
    assert all(np.isfinite(x) and x >= 0 for x in errors.values())
    # This specific short experiment should show decreasing max drift.
    assert errors[0.001] < errors[0.002] < errors[0.004], errors
    print('  max |E(t)-E(0)|:', {k: f'{v:.3e}' for k,v in errors.items()})
    print('  PASS')

    print('[5] Work-energy balance with applied torque')
    t, q, qd, tau = simulate(robot, np.array([0.3, -0.2]), np.array([0.1, 0.2]),
                             lambda t,q,qd: np.array([1.5, -0.7]),
                             dt=0.001, duration=0.5)
    r = energy_balance_residual(robot, t, q, qd, tau)
    assert r.shape == t.shape and np.isfinite(r).all()
    assert abs(r[-1]) < 0.1  # nonzero discretization error is expected
    print(f'  terminal energy-work residual = {r[-1]:.4e} J')
    print('  PASS')
    print('\nDay 15: all tests passed.')


def optional_plots():
    """Optional figures; install matplotlib if needed: pip install matplotlib."""
    import matplotlib.pyplot as plt
    robot = Planar2R()
    q0, v0 = np.array([0.4, -0.6]), np.array([0.2, 0.1])
    fig, ax = plt.subplots()
    for dt in [0.004, 0.002, 0.001]:
        t, q, qd, tau = simulate(robot, q0, v0, lambda t,q,qd: np.zeros(2), dt, 0.5)
        E = np.array([total_energy(robot, qi, vi) for qi,vi in zip(q, qd)])
        ax.plot(t, E-E[0], label=f'dt={dt:g}')
    ax.set(xlabel='time (s)', ylabel='E(t)-E(0) (J)', title='Day 15: free-motion energy drift')
    ax.legend()
    ax.grid(True)
    plt.show()


if __name__ == '__main__':
    run_tests()
    # After all TODOs pass, uncomment for an optional energy-drift plot:
    # optional_plots()
