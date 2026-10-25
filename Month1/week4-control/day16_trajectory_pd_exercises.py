"""Day 16 | Quintic joint trajectory and PD tracking (NumPy only).

Complete TODO 1-5, then run:
    python day16_trajectory_pd_exercises.py

Optional: python day16_trajectory_pd_exercises.py --plot

Sign convention: e = q_des - q, edot = qd_des - qd.
Robot: vertical XY plane, gravity along -Y, joint angles in radians.
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class Planar2R:
    L1: float = 2.0
    L2: float = 1.0
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
    a = (robot.m1*robot.r1 + robot.m2*robot.L1) * np.cos(q1)
    b = robot.m2*robot.r2*np.cos(q1+q2)
    return robot.gravity * np.array([a+b, b], dtype=float)


def inverse_dynamics(robot: Planar2R, q: np.ndarray, qd: np.ndarray, qdd: np.ndarray) -> np.ndarray:
    return mass_matrix(robot, q) @ qdd + coriolis_vector(robot, q, qd) + gravity_vector(robot, q)


def forward_dynamics(robot: Planar2R, q: np.ndarray, qd: np.ndarray, tau: np.ndarray) -> np.ndarray:
    rhs = np.asarray(tau) - coriolis_vector(robot, q, qd) - gravity_vector(robot, q)
    return np.linalg.solve(mass_matrix(robot, q), rhs)


# -----------------------------------------------------------------------------
# TODO 1: quintic trajectory
# -----------------------------------------------------------------------------
def quintic_trajectory(t: float, q0: np.ndarray, qf: np.ndarray, T: float):
    """Return (q_des, qd_des, qdd_des), each shape (2,).

    Rest-to-rest with zero endpoint acceleration:
        q(0)=q0, q(T)=qf,
        qdot(0)=qdot(T)=0,
        qddot(0)=qddot(T)=0.
    For t<0, return (q0,0,0); for t>T, return (qf,0,0).
    Require T>0. Use s=t/T and vectorized NumPy expressions.
    """
    if t < 0: 
        return q0, np.zeros(2), np.zeros(2)
    if t > T:
        return qf, np.zeros(2), np.zeros(2)
    s = t / T
    q_des = q0 + (qf - q0) * (10*s**3 - 15*s**4 + 6*s**5)
    qd_des = (qf - q0) / T * (30*s**2 - 60*s**3 + 30*s**4)
    qdd_des = (qf - q0) /T**2 * (60*s - 180*s**2 + 120*s**3)
    return q_des, qd_des, qdd_des
    


# -----------------------------------------------------------------------------
# TODO 2: pure PD
# -----------------------------------------------------------------------------
def control_pd(robot, q, qd, q_des, qd_des, qdd_des, Kp, Kd):
    """Return tau from pure PD tracking. qdd_des is intentionally unused.

    e=q_des-q, edot=qd_des-qd; Kp and Kd are (2,2) gain matrices.
    """
    e = q_des - q 
    edot = qd_des - qd 
    tau = Kp @ e + Kd @ edot 
    return tau 
    


# -----------------------------------------------------------------------------
# TODO 3: inverse dynamics feedforward + PD
# -----------------------------------------------------------------------------
def control_ff_pd(robot, q, qd, q_des, qd_des, qdd_des, Kp, Kd):
    """tau = inverse_dynamics(q_des, qd_des, qdd_des) + PD feedback.

    IMPORTANT: inverse dynamics feedforward is evaluated at DESIRED state.
    """
    M = mass_matrix(robot, q_des) 
    c = coriolis_vector(robot, q_des, qd_des)
    g = gravity_vector(robot, q_des) 
    tau = M @ qdd_des + c + g + Kp @ (q_des - q) + Kd @ (qd_des - qd) 
    return tau 


# -----------------------------------------------------------------------------
# TODO 4: gravity compensation at actual position + PD
# -----------------------------------------------------------------------------
def control_gravity_pd(robot, q, qd, q_des, qd_des, qdd_des, Kp, Kd):
    """tau = gravity_vector(robot,q) + PD feedback.

    This is not full inverse dynamics feedforward: no inertia/Coriolis feedforward.
    """
    g = gravity_vector(robot, q)  
    tau = g + Kp @ (q_des - q) + Kd @ (qd_des - qd) 
    return tau 


# -----------------------------------------------------------------------------
# Simulation (provided; semi-implicit Euler as in Day 15)
# -----------------------------------------------------------------------------
def simulate_tracking(robot, controller, q0, qd0, traj_q0, traj_qf, T, Kp, Kd, dt=0.001, settle=1.0):
    """Track for T seconds, then hold target for `settle` seconds.

    Returns dict with t, q, qd, q_des, qd_des, tau, tracking_rmse, final_error.
    Input arrays are not mutated. No actuator torque limits in this exercise.
    """
    duration = T + settle
    n = round(duration / dt)
    assert np.isclose(n*dt, duration), "duration must be divisible by dt"
    times = np.arange(n+1)*dt
    qs = np.zeros((n+1, 2)); qds = np.zeros((n+1, 2))
    qrefs = np.zeros((n+1, 2)); vrefs = np.zeros((n+1, 2))
    taus = np.zeros((n, 2))
    qs[0] = q0; qds[0] = qd0
    for k, t in enumerate(times):
        qref, vref, aref = quintic_trajectory(t, traj_q0, traj_qf, T)
        qrefs[k], vrefs[k] = qref, vref
        if k == n:
            break
        tau = controller(robot, qs[k], qds[k], qref, vref, aref, Kp, Kd)
        taus[k] = tau
        qdd = forward_dynamics(robot, qs[k], qds[k], tau)
        qds[k+1] = qds[k] + dt*qdd
        qs[k+1] = qs[k] + dt*qds[k+1]
    errors = qs-qrefs
    return dict(t=times, q=qs, qd=qds, q_des=qrefs, qd_des=vrefs,
                tau=taus, tracking_rmse=float(np.sqrt(np.mean(errors**2))),
                final_error=float(np.linalg.norm(errors[-1])))


# -----------------------------------------------------------------------------
# TODO 5: compare controllers
# -----------------------------------------------------------------------------
def compare_controllers(robot, q_start, q_goal, T, Kp, Kd, dt=0.001):
    """Return a dictionary mapping controller names to simulation results.

    Use these exact keys: 'PD', 'FF+PD', 'Gravity+PD'.
    Use initial q0=q_start, qd0=zeros(2), desired start=q_start,
    desired finish=q_goal, and settle=1.0 for all three runs.
    Each value is the dict returned by simulate_tracking.
    """
    results = {} 
    results['PD'] = simulate_tracking(robot, control_pd, q_start, np.zeros(2), q_start, q_goal, T, Kp, Kd, dt)
    results['FF+PD'] = simulate_tracking(robot, control_ff_pd, q_start, np.zeros(2), q_start, q_goal, T, Kp, Kd, dt)
    results['Gravity+PD'] = simulate_tracking(robot, control_gravity_pd, q_start, np.zeros(2), q_start, q_goal, T, Kp, Kd, dt)
    return results 


# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------
def run_tests():
    robot = Planar2R()
    q0 = np.array([0.0, -0.2]); qf = np.array([0.8, 0.45]); T = 1.5
    Kp = np.diag([120.0, 65.0]); Kd = np.diag([24.0, 14.0])
    print('[1] Quintic trajectory boundary conditions and derivatives')
    for t, q_expected in [(0., q0), (T, qf)]:
        q, v, a = quintic_trajectory(t, q0, qf, T)
        np.testing.assert_allclose(q, q_expected, atol=1e-12)
        np.testing.assert_allclose(v, np.zeros(2), atol=1e-12)
        np.testing.assert_allclose(a, np.zeros(2), atol=1e-12)
    for t in [0.15, 0.5, 0.9, 1.3]:
        h = 1e-5
        q, v, a = quintic_trajectory(t, q0, qf, T)
        qplus, vplus, _ = quintic_trajectory(t+h, q0, qf, T)
        qminus, vminus, _ = quintic_trajectory(t-h, q0, qf, T)
        np.testing.assert_allclose((qplus-qminus)/(2*h), v, atol=2e-8)
        np.testing.assert_allclose((vplus-vminus)/(2*h), a, atol=2e-8)
    for t, q_expected in [(-1., q0), (T+1., qf)]:
        q, v, a = quintic_trajectory(t, q0, qf, T)
        np.testing.assert_allclose(q, q_expected, atol=1e-12)
        np.testing.assert_allclose(v, np.zeros(2), atol=1e-12)
        np.testing.assert_allclose(a, np.zeros(2), atol=1e-12)
    print('  PASS')

    print('[2] Pure PD tracking torque')
    q=np.array([0.2,-0.1]); qd=np.array([0.3,-0.2]); qr=np.array([0.4,0.15]); vr=np.array([-0.1,0.2]); ar=np.array([1.,-1.])
    pd_ref=Kp@(qr-q)+Kd@(vr-qd)
    np.testing.assert_allclose(control_pd(robot,q,qd,qr,vr,ar,Kp,Kd),pd_ref,atol=1e-12)
    print('  PASS')

    print('[3] Inverse dynamics feedforward + PD torque')
    ff_ref=inverse_dynamics(robot,qr,vr,ar)+pd_ref
    np.testing.assert_allclose(control_ff_pd(robot,q,qd,qr,vr,ar,Kp,Kd),ff_ref,atol=1e-12)
    np.testing.assert_allclose(control_ff_pd(robot,qr,vr,qr,vr,ar,Kp,Kd),inverse_dynamics(robot,qr,vr,ar),atol=1e-12)
    print('  PASS')

    print('[4] Actual-state gravity compensation + PD torque')
    g_ref=gravity_vector(robot,q)+pd_ref
    np.testing.assert_allclose(control_gravity_pd(robot,q,qd,qr,vr,ar,Kp,Kd),g_ref,atol=1e-12)
    print('  PASS')

    print('[5] Closed-loop controller comparison')
    results=compare_controllers(robot,q0,qf,T,Kp,Kd,dt=0.001)
    assert set(results)=={'PD','FF+PD','Gravity+PD'}
    for name, r in results.items():
        assert np.all(np.isfinite(r['q'])) and np.all(np.isfinite(r['tau']))
        assert r['q'].shape==r['q_des'].shape and r['tau'].shape[0]==len(r['t'])-1
        print(f"  {name:12s} RMSE={r['tracking_rmse']:.5f} rad  "
              f"final_error={r['final_error']:.5f} rad  "
              f"max|tau|={np.max(np.abs(r['tau'])):.2f} N m")
    assert results['FF+PD']['tracking_rmse'] < results['PD']['tracking_rmse']
    assert results['Gravity+PD']['final_error'] < results['PD']['final_error']
    print('  PASS')
    print('\nDay 16: all tests passed.')
    return results


def optional_plots(results):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print('matplotlib not installed; skipping optional plots')
        return
    fig, axes=plt.subplots(2,1,figsize=(9,7),sharex=True)
    for name,r in results.items():
        err=r['q']-r['q_des']
        for j in range(2):
            axes[j].plot(r['t'],err[:,j],label=name)
            axes[j].set_ylabel(f'joint {j+1} error (rad)')
            axes[j].grid(alpha=.25)
            axes[j].legend()
    axes[-1].set_xlabel('time (s)')
    fig.tight_layout()
    plt.show()


if __name__ == '__main__':
    results=run_tests()
    optional_plots(results)
