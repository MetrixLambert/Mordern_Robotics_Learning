"""Day 21 — 2R robot: computed-torque control under actuator saturation.

NumPy only. Fill TODO 1–4 and run:
    python day21_torque_saturation_exercises.py

Convention: e=q-q_des, ed=qd-qd_des. Gravity is included.
The true dynamics are M(q)qdd + c(q,qd) + g(q) = tau_actual.
"""
import numpy as np


def mass_matrix(q):
    q2 = float(q[1])
    return np.array([[6.60 + 2.0*np.cos(q2), 0.35 + np.cos(q2)],
                     [0.35 + np.cos(q2), 0.35]], dtype=float)


def coriolis_vector(q, qd):
    v1, v2 = qd
    s = np.sin(q[1])
    return s * np.array([-2*v1*v2-v2*v2, v1*v1])


def gravity_vector(q):
    q1, q2 = q
    return 9.81*np.array([4*np.cos(q1)+0.5*np.cos(q1+q2),
                          0.5*np.cos(q1+q2)])


def quintic_reference(t, duration=2.0, start=None, goal=None):
    """Provided rest-to-rest joint-space reference; returns (q,qd,qdd)."""
    if start is None: start = np.array([0.0, 0.0])
    if goal is None: goal = np.array([0.8, -0.6])
    s = np.clip(t/duration, 0.0, 1.0)
    h = 10*s**3 - 15*s**4 + 6*s**5
    hd = (30*s**2 - 60*s**3 + 30*s**4)/duration if t < duration else 0.0
    hdd = (60*s - 180*s**2 + 120*s**3)/duration**2 if t < duration else 0.0
    dq = goal-start
    return start+dq*h, dq*hd, dq*hdd


# TODO 1 --------------------------------------------------------------------
def computed_torque(q, qd, q_des, qd_des, qdd_des, kp=25.0, kd=10.0):
    """Compute unconstrained CTC command (shape (2,)).

    v=qdd_des-kd*(qd-qd_des)-kp*(q-q_des)
    tau_cmd=M(q)@v+c(q,qd)+g(q)
    """
    v = qdd_des - kd * (qd - qd_des) - kp * (q - q_des)
    tau_cmd = mass_matrix(q) @ v + coriolis_vector(q, qd) + gravity_vector(q)
    return tau_cmd


# TODO 2 --------------------------------------------------------------------
def saturate_torque(tau_cmd, tau_limit):
    """Componentwise clipping to [-tau_limit_i,+tau_limit_i].

    tau_limit can be scalar or a length-2 vector; np.asarray is useful.
    Return (tau_actual, delta_tau), where delta_tau=tau_actual-tau_cmd.
    Do not modify the input array.
    """
    tau_limit = np.asarray(tau_limit)
    tau_actual = np.clip(tau_cmd, -tau_limit, tau_limit) 
    delta_tau = tau_actual - tau_cmd
    return tau_actual, delta_tau


# TODO 3 --------------------------------------------------------------------
def actual_error_acceleration(q, qd, q_des, qd_des, qdd_des,
                              tau_cmd, tau_actual, kp=25.0, kd=10.0):
    """Return the *actual* tracking error acceleration using the CTC identity:

    e_dd = -kd*ed-kp*e + solve(M(q), tau_actual-tau_cmd).
    (Use np.linalg.solve, not np.linalg.inv.)
    """
    e_dd = -kd * (qd - qd_des) - kp * (q - q_des) + np.linalg.solve(mass_matrix(q), tau_actual - tau_cmd)
    return e_dd


def simulate(tau_limit, kp=25.0, kd=10.0, dt=0.002, duration=4.0,
             initial_error=(0.15, -0.10)):
    """Provided semi-implicit Euler simulation of TRUE 2R dynamics.

    No clipping of gravity/feedback separately: saturate total commanded torque.
    """
    assert dt > 0 and duration > 0 and kp > 0 and kd > 0
    n = int(round(duration/dt))
    t = np.arange(n+1)*dt
    q = np.zeros((n+1, 2)); qd = np.zeros_like(q)
    q_des = np.zeros_like(q)
    tau_cmd = np.zeros((n,2)); tau_actual = np.zeros((n,2))
    saturated = np.zeros((n,2), dtype=bool)
    q0, v0, _ = quintic_reference(0.0)
    q[0] = q0 + np.asarray(initial_error)
    qd[0] = v0
    for k in range(n+1):
        qr, vr, ar = quintic_reference(t[k])
        q_des[k] = qr
        if k == n: break
        cmd = computed_torque(q[k], qd[k], qr, vr, ar, kp, kd)
        actual, delta = saturate_torque(cmd, tau_limit)
        tau_cmd[k] = cmd
        tau_actual[k] = actual
        saturated[k] = np.abs(delta) > 1e-10
        qdd = np.linalg.solve(mass_matrix(q[k]),
                              actual-coriolis_vector(q[k],qd[k])-gravity_vector(q[k]))
        qd[k+1] = qd[k]+dt*qdd
        q[k+1] = q[k]+dt*qd[k+1]
    e = q-q_des
    return dict(t=t, q=q, qd=qd, q_des=q_des, error=e,
                tau_cmd=tau_cmd, tau_actual=tau_actual, saturated=saturated,
                rmse=float(np.sqrt(np.mean(np.sum(e*e, axis=1)))),
                final_error=float(np.linalg.norm(e[-1])),
                saturation_fraction=float(np.mean(np.any(saturated,axis=1))))


# TODO 4 --------------------------------------------------------------------
def compare_saturation_cases():
    """Return exactly these four runs (use simulate() defaults otherwise):

    'unlimited' : tau_limit=np.inf, kp=25, kd=10
    'moderate'  : tau_limit=[80, 15], kp=25, kd=10
    'tight'     : tau_limit=[45,  8], kp=25, kd=10
    'high_gain' : tau_limit=[45,  8], kp=100, kd=20

    Note: gravity at q=(0,0) is about [44.145,4.905] Nm;
    limits must constrain TOTAL motor torque, including gravity.
    """
    results = {} 
    results['unlimited'] = simulate(np.inf, 25, 10)
    results['moderate'] = simulate(np.array([80, 15]), 25, 10)
    results['tight'] = simulate(np.array([45, 8]), 25, 10)
    results['high_gain'] = simulate(np.array([45, 8]), 100, 20)
    return results


def run_tests():
    print('[1] CTC command and zero-error feedforward')
    q=np.array([0.3,-0.4]); z=np.zeros(2)
    got=computed_torque(q,z,q,z,z)
    np.testing.assert_allclose(got,gravity_vector(q),atol=1e-10)
    got2=computed_torque(q+np.array([0.1,-0.2]),z,q,z,z)
    assert not np.allclose(got2,got)
    print('  PASS')

    print('[2] Componentwise saturation and sign')
    actual, delta=saturate_torque(np.array([100.,-20.]),np.array([45.,8.]))
    np.testing.assert_allclose(actual,[45.,-8.])
    np.testing.assert_allclose(delta,[-55.,12.])
    actual2,delta2=saturate_torque(np.array([1.,-2.]),np.inf)
    np.testing.assert_allclose(actual2,[1.,-2.])
    np.testing.assert_allclose(delta2,[0.,0.])
    print('  PASS')

    print('[3] Actual error acceleration identity')
    q=np.array([0.3,-0.5]); qd=np.array([0.2,-0.1])
    qr=np.array([0.1,-0.3]); vr=np.array([0.05,0.04]); ar=np.array([0.3,-0.2])
    cmd=computed_torque(q,qd,qr,vr,ar)
    actual,_=saturate_torque(cmd,[40.,6.])
    err_dd=actual_error_acceleration(q,qd,qr,vr,ar,cmd,actual)
    true_qdd=np.linalg.solve(mass_matrix(q),actual-coriolis_vector(q,qd)-gravity_vector(q))
    np.testing.assert_allclose(err_dd,true_qdd-ar,atol=1e-10)
    print('  PASS')

    print('[4] Four saturation experiments')
    runs=compare_saturation_cases()
    assert set(runs)=={'unlimited','moderate','tight','high_gain'}
    for name, r in runs.items():
        assert np.isfinite(r['rmse']) and np.isfinite(r['final_error'])
        assert np.isfinite(r['q']).all()
        print(f"  {name:10s} RMSE={r['rmse']:.5f} rad  "
              f"final={r['final_error']:.5f} rad  "
              f"sat={100*r['saturation_fraction']:.1f}%")
    assert runs['unlimited']['saturation_fraction']==0.0
    assert runs['tight']['saturation_fraction']>0.0
    assert runs['high_gain']['saturation_fraction']>0.0
    assert np.max(np.abs(runs['tight']['tau_actual'][:,0]))<=45+1e-10
    assert np.max(np.abs(runs['tight']['tau_actual'][:,1]))<=8+1e-10
    print('  PASS')
    print('\nDay 21: all tests passed.')


if __name__=='__main__':
    run_tests()
