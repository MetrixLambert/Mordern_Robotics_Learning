"""Day 17 | Computed-torque regulation with gravity-model mismatch (NumPy only).

Fill TODO 1-4 and run:
    python day17_gravity_mismatch_exercises.py

Convention: e = q - q_target (opposite sign to Day 16's controller helper).
True plant uses `robot`; controller uses `gravity_scale` * true gravity.
No friction, torque saturation, or sensor noise. Semi-implicit Euler integration.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Planar2R:
    L1: float = 2.0
    m1: float = 2.0
    m2: float = 1.0
    r1: float = 1.0
    r2: float = 0.5
    I1: float = 0.25
    I2: float = 0.10
    gravity: float = 9.81


def mass_matrix(robot: Planar2R, q: np.ndarray) -> np.ndarray:
    q2 = float(q[1])
    a = robot.I1 + robot.I2 + robot.m1 * robot.r1**2 + robot.m2 * (robot.L1**2 + robot.r2**2)
    b = robot.m2 * robot.L1 * robot.r2
    d = robot.I2 + robot.m2 * robot.r2**2
    return np.array([[a + 2*b*np.cos(q2), d + b*np.cos(q2)],
                     [d + b*np.cos(q2), d]], dtype=float)


def coriolis_vector(robot: Planar2R, q: np.ndarray, qd: np.ndarray) -> np.ndarray:
    h = robot.m2 * robot.L1 * robot.r2 * np.sin(float(q[1]))
    v1, v2 = np.asarray(qd, dtype=float)
    return np.array([-h*(2*v1*v2 + v2**2), h*v1**2], dtype=float)


def gravity_vector(robot: Planar2R, q: np.ndarray) -> np.ndarray:
    q1, q2 = np.asarray(q, dtype=float)
    a = (robot.m1*robot.r1 + robot.m2*robot.L1)*np.cos(q1)
    b = robot.m2*robot.r2*np.cos(q1+q2)
    return robot.gravity*np.array([a+b, b], dtype=float)


def forward_dynamics(robot: Planar2R, q: np.ndarray, qd: np.ndarray, tau: np.ndarray) -> np.ndarray:
    rhs = tau - coriolis_vector(robot, q, qd) - gravity_vector(robot, q)
    return np.linalg.solve(mass_matrix(robot, q), rhs)


# TODO 1 --------------------------------------------------------------
def computed_torque_regulation(robot, q, qd, q_target, Kp, Kd, gravity_scale=1.0):
    """Implement computed-torque position regulation:

    e = q - q_target
    v = -Kd @ qd - Kp @ e
    tau = M(q) @ v + c(q,qd) + gravity_scale * g(q)

    Only gravity compensation is mismatched. Return shape (2,).
    """
    e = q - q_target 
    v = -Kd @ qd - Kp @ e 
    M = mass_matrix(robot, q) 
    c = coriolis_vector(robot, q, qd)
    g = gravity_vector(robot, q) 
    tau = M @ v + c + gravity_scale * g 
    return tau 


# Simulation (provided) -----------------------------------------------
def simulate(robot, controller, q0, q_target, Kp, Kd, dt=0.002, duration=8.0, **kwargs):
    """Return time, q, qd, tau; true robot dynamics used for propagation."""
    n = int(round(duration/dt))
    assert dt > 0 and np.isclose(n*dt, duration)
    t = np.arange(n+1)*dt
    q = np.zeros((n+1, 2)); qd = np.zeros((n+1, 2))
    tau = np.zeros((n, 2))
    q[0] = np.asarray(q0, dtype=float)
    for k in range(n):
        tau[k] = controller(robot, q[k], qd[k], q_target, Kp, Kd, **kwargs)
        qdd = forward_dynamics(robot, q[k], qd[k], tau[k])
        qd[k+1] = qd[k] + dt*qdd
        q[k+1] = q[k] + dt*qd[k+1]
    return dict(t=t, q=q, qd=qd, tau=tau)


# TODO 2 --------------------------------------------------------------
def compare_gravity_models(robot, q0, q_target, Kp, Kd):
    """Simulate two cases, both with dt=0.002 and duration=8.0:

    - 'exact': gravity_scale=1.0
    - 'mismatch_10pct': gravity_scale=1.1

    Return dict with exactly those two keys and values from simulate().
    """
    dt = 0.002 
    duration = 8.0 
    results = {} 
    results['exact'] = simulate(robot, computed_torque_regulation, q0, q_target, Kp, Kd, dt, duration, gravity_scale=1.0)
    results['mismatch_10pct'] = simulate(robot, computed_torque_regulation, q0, q_target, Kp, Kd, dt, duration, gravity_scale=1.1)
    return results 


# TODO 3 --------------------------------------------------------------
def equilibrium_residual(robot, q_ss, q_target, Kp, gravity_scale):
    """Compute the static equilibrium residual (2-vector):

    M(q_ss) @ Kp @ (q_ss - q_target) - (gravity_scale-1)*g(q_ss)

    Note: d(q) = (gravity_scale-1)*g(q) is state dependent, not constant.
    """
    M = mass_matrix(robot, q_ss) 
    qd_ss = np.zeros(2)
    c = coriolis_vector(robot, q_ss, qd_ss)
    g = gravity_vector(robot, q_ss) 
    equilibrium_residual = M @ Kp @ (q_ss - q_target) - (gravity_scale-1) * g 
    return equilibrium_residual 


# TODO 4 --------------------------------------------------------------
def summarize_results(robot, results, q_target, Kp, gravity_scale=1.1):
    """Return dict with keys:
      'exact_final_error': norm(q_exact[-1]-q_target)
      'mismatch_final_error': norm(q_mismatch[-1]-q_target)
      'mismatch_final_speed': norm(qd_mismatch[-1])
      'mismatch_equilibrium_residual': norm(equilibrium_residual(...))

    Use final states from results. Do not solve for equilibrium numerically.
    """
    exact_final_error = np.linalg.norm(results['exact']['q'][-1] - q_target)
    mismatch_final_error = np.linalg.norm(results['mismatch_10pct']['q'][-1] - q_target)
    mismatch_final_speed = np.linalg.norm(results['mismatch_10pct']['qd'][-1])
    mismatch_equilibrium_residual = np.linalg.norm(equilibrium_residual(robot, results['mismatch_10pct']['q'][-1], q_target, Kp, 1.1))
    return {
        'exact_final_error': exact_final_error,
        'mismatch_final_error': mismatch_final_error,
        'mismatch_final_speed': mismatch_final_speed,
        'mismatch_equilibrium_residual': mismatch_equilibrium_residual
    }


def run_tests():
    robot = Planar2R()
    q0 = np.array([-0.25, -0.3])
    q_target = np.array([0.55, 0.25])
    Kp = np.diag([36.0, 36.0])
    Kd = np.diag([14.0, 14.0])

    print('[1] Computed-torque control torque')
    q = np.array([0.2, -0.4]); qd = np.array([0.1, -0.2])
    v = -Kd@qd - Kp@(q-q_target)
    expected = mass_matrix(robot,q)@v + coriolis_vector(robot,q,qd) + 1.1*gravity_vector(robot,q)
    np.testing.assert_allclose(computed_torque_regulation(robot,q,qd,q_target,Kp,Kd,1.1),expected,atol=1e-11)
    print('  PASS')

    print('[2] Simulation comparison')
    results = compare_gravity_models(robot,q0,q_target,Kp,Kd)
    assert set(results)=={'exact','mismatch_10pct'}
    for name, r in results.items():
        assert r['q'].shape == (4001,2) and r['qd'].shape == (4001,2)
        assert np.all(np.isfinite(r['q'])) and np.all(np.isfinite(r['tau']))
        print(f"  {name:15s} final error = {np.linalg.norm(r['q'][-1]-q_target):.6f} rad")
    print('  PASS')

    print('[3] Static equilibrium residual')
    q_ss = results['mismatch_10pct']['q'][-1]
    expected_residual = mass_matrix(robot,q_ss)@Kp@(q_ss-q_target)-0.1*gravity_vector(robot,q_ss)
    np.testing.assert_allclose(equilibrium_residual(robot,q_ss,q_target,Kp,1.1),expected_residual,atol=1e-10)
    print('  PASS')

    print('[4] End-state checks')
    metrics = summarize_results(robot,results,q_target,Kp)
    assert set(metrics)=={'exact_final_error','mismatch_final_error','mismatch_final_speed','mismatch_equilibrium_residual'}
    assert metrics['exact_final_error'] < 1e-3
    assert metrics['mismatch_final_error'] > 1e-3
    assert metrics['mismatch_final_speed'] < 1e-3
    assert metrics['mismatch_equilibrium_residual'] < 1e-2
    for name,value in metrics.items():
        print(f'  {name}: {value:.7g}')
    print('  PASS\n\nDay 17: all tests passed.')
    return results


if __name__ == '__main__':
    run_tests()
