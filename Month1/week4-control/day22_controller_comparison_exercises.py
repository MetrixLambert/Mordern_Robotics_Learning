"""Day 22 — Four 2R robot controllers: unified comparison (NumPy only).

Run: python day22_controller_comparison_exercises.py
Complete TODO 1–4. All other functions and tests are provided.

Convention: e = q - q_des, ed = qd - qd_des.
Real plant: M(q) qdd + c(q,qd) + g(q) = tau_actual.
Controller model may be mismatched; saturation is applied to TOTAL torque.

Controllers:
  pd     : -Kp e - Kd ed
  gpd    : g_hat(q) - Kp e - Kd ed
  ffpd   : M_hat(q_des) qdd_des + c_hat(q_des,qd_des) + g_hat(q_des)
           - Kp e - Kd ed
  ctc    : M_hat(q) (qdd_des - Kd ed - Kp e) + c_hat(q,qd) + g_hat(q)

IMPORTANT: Kp/Kd have torque units in PD/GPD/FFPD but acceleration
units in CTC. Numerical gains are NOT intrinsically comparable.
"""
import numpy as np


def mass_matrix(q):
    q2 = float(q[1])
    return np.array([[6.60 + 2.0*np.cos(q2), 0.35 + np.cos(q2)],
                     [0.35 + np.cos(q2), 0.35]], dtype=float)


def coriolis_vector(q, qd):
    v1, v2 = qd
    return np.sin(q[1]) * np.array([-2*v1*v2-v2*v2, v1*v1])


def gravity_vector(q):
    q1, q2 = q
    return 9.81*np.array([4*np.cos(q1)+0.5*np.cos(q1+q2),
                          0.5*np.cos(q1+q2)])


def quintic_reference(t, duration=2.0, start=None, goal=None):
    if start is None: start = np.array([0.0, 0.0])
    if goal is None: goal = np.array([0.8, -0.6])
    s = np.clip(t/duration, 0.0, 1.0)
    h = 10*s**3 - 15*s**4 + 6*s**5
    hd = (30*s**2 - 60*s**3 + 30*s**4)/duration if t < duration else 0.0
    hdd = (60*s - 180*s**2 + 120*s**3)/duration**2 if t < duration else 0.0
    dq = goal-start
    return start+dq*h, dq*hd, dq*hdd


def desired_state(t, task):
    if task == 'tracking':
        return quintic_reference(t)
    if task == 'regulation':
        return np.array([0.8, -0.6]), np.zeros(2), np.zeros(2)
    raise ValueError('task must be tracking or regulation')


def model_terms(q, qd, model_scale=(1., 1., 1.)):
    """Return nominal (M_hat, c_hat, g_hat).

    Each scale is a deliberately simplified mismatch knob; not necessarily
    a physically consistent alternative robot parameterization.
    """
    a, b, d = model_scale
    return a*mass_matrix(q), b*coriolis_vector(q, qd), d*gravity_vector(q)


# TODO 1 --------------------------------------------------------------------
def controller_torque(kind, q, qd, qr, vr, ar, kp=25., kd=10.,
                      model_scale=(1., 1., 1.)):
    """Return unconstrained total torque command, shape (2,).

    kind in {'pd','gpd','ffpd','ctc'}; use definitions in module docstring.
    Use model_terms(...) for estimated dynamics, not true terms directly.
    Scalar kp/kd are sufficient; NumPy broadcasting supports vectors too.
    Raise ValueError for unknown controller kind.
    """
    e = q - qr 
    ed = qd - vr 
    
    if kind == 'pd':
        return -kp * e - kd * ed
    elif kind == 'gpd': 
        M_hat, c_hat, g_hat = model_terms(q, qd, model_scale)
        return g_hat - kp * e - kd * ed 
    elif kind == 'ffpd': 
        M_hat, c_hat, g_hat = model_terms(qr, vr, model_scale)
        return M_hat @ ar + c_hat + g_hat - kp * e - kd * ed 
    elif kind == 'ctc': 
        M_hat, c_hat, g_hat = model_terms(q, qd, model_scale)
        v = ar - kd * ed - kp * e 
        return M_hat @ v + c_hat + g_hat 
    else: 
        raise ValueError('unknown controller kind')



# TODO 2 --------------------------------------------------------------------
def saturate_torque(tau_cmd, tau_limit):
    """Componentwise saturation; return (tau_actual, saturation_mask).

    tau_limit is scalar or length-2 vector, can be np.inf.
    saturation_mask[i] is True iff actual torque differs from command.
    Never mutate tau_cmd.
    """
    tau_limit = np.asarray(tau_limit)
    tau_actual = np.clip(tau_cmd, -tau_limit, tau_limit) 
    saturation_mask = tau_actual != tau_cmd
    return tau_actual, saturation_mask


def simulate(kind, task='tracking', model_scale=(1., 1., 1.),
             tau_limit=np.inf, kp=25., kd=10., duration=4., dt=0.002,
             initial_error=(0.15, -0.10)):
    """Provided semi-implicit Euler integration of TRUE robot dynamics."""
    if kind not in {'pd','gpd','ffpd','ctc'}:
        raise ValueError(kind)
    if task not in {'tracking','regulation'}:
        raise ValueError(task)
    if duration <= 0 or dt <= 0:
        raise ValueError('duration and dt must be positive')
    n = int(round(duration/dt))
    t = np.arange(n+1)*dt
    q = np.zeros((n+1, 2)); qd = np.zeros_like(q)
    q_des = np.zeros_like(q)
    tau_cmd = np.zeros((n, 2)); tau_actual = np.zeros_like(tau_cmd)
    saturated = np.zeros((n, 2), dtype=bool)
    q0, v0, _ = desired_state(0., task)
    q[0] = q0 + np.asarray(initial_error)
    qd[0] = v0
    for k in range(n+1):
        qr, vr, ar = desired_state(t[k], task)
        q_des[k] = qr
        if k == n:
            break
        cmd = controller_torque(kind, q[k], qd[k], qr, vr, ar,
                                kp=kp, kd=kd, model_scale=model_scale)
        actual, mask = saturate_torque(cmd, tau_limit)
        tau_cmd[k], tau_actual[k], saturated[k] = cmd, actual, mask
        qdd = np.linalg.solve(mass_matrix(q[k]),
                              actual-coriolis_vector(q[k], qd[k])-gravity_vector(q[k]))
        qd[k+1] = qd[k] + dt*qdd
        q[k+1] = q[k] + dt*qd[k+1]
    e = q-q_des
    return dict(t=t, q=q, qd=qd, q_des=q_des, error=e,
                tau_cmd=tau_cmd, tau_actual=tau_actual, saturated=saturated,
                rmse=float(np.sqrt(np.mean(np.sum(e**2, axis=1)))),
                final_error=float(np.linalg.norm(e[-1])),
                peak_actual_torque=float(np.max(np.abs(tau_actual))),
                saturation_fraction=float(np.mean(np.any(saturated, axis=1))))


# TODO 3 --------------------------------------------------------------------
def run_comparison():
    """Return nested dict: results[scenario][controller] = simulate result.

    Run all four controllers for each scenario, in this order:
      'ideal'      : tracking, exact model, unlimited torque
      'mismatch'   : tracking, model_scale=(1.15, 1.15, 0.85), unlimited
      'saturation' : tracking, exact model, tau_limit=[45., 8.]
      'regulation' : regulation, exact model, unlimited

    Use kp=25., kd=10. and all other simulate defaults in every run.
    This is a SAME-NUMERICAL-GAIN baseline, not a fair gain-tuned ranking.
    """
    results = {} 
    for scenario in ['ideal', 'mismatch', 'saturation', 'regulation']:
        results[scenario] = {}
        for controller in ['pd', 'gpd', 'ffpd', 'ctc']:
            if scenario == 'ideal':
                results[scenario][controller] = simulate(controller, task='tracking', model_scale=(1., 1., 1.), tau_limit=np.inf, kp=25., kd=10.)
            elif scenario == 'mismatch':
                results[scenario][controller] = simulate(controller, task='tracking', model_scale=(1.15, 1.15, 0.85), tau_limit=np.inf, kp=25., kd=10.)
            elif scenario == 'saturation':
                results[scenario][controller] = simulate(controller, task='tracking', model_scale=(1., 1., 1.), tau_limit=np.array([45., 8.]), kp=25., kd=10.)
            elif scenario == 'regulation':
                results[scenario][controller] = simulate(controller, task='regulation', model_scale=(1., 1., 1.), tau_limit=np.inf, kp=25., kd=10.)
    return results


# TODO 4 --------------------------------------------------------------------
def summarize(results):
    """Return summary[scenario][controller] as a length-4 tuple:

    (rmse, final_error, peak_actual_torque, saturation_fraction).
    Do not recompute simulations; extract fields from the supplied results.
    """
    summary = {} 
    for scenario in results:
        summary[scenario] = {}
        for controller in results[scenario]:
            summary[scenario][controller] = (results[scenario][controller]['rmse'], results[scenario][controller]['final_error'], results[scenario][controller]['peak_actual_torque'], results[scenario][controller]['saturation_fraction'])
    return summary


# ---------------------------- TESTS ---------------------------------------
def _close(a, b, tol=1e-9):
    return np.allclose(a, b, rtol=tol, atol=tol)


def test_controllers():
    q=np.array([0.2,-0.3]); qd=np.array([0.1,-0.2])
    qr=np.array([0.4,-0.5]); vr=np.array([0.2,0.1]); ar=np.array([0.3,-0.1])
    e=q-qr; ed=qd-vr
    fb=-25*e-10*ed
    assert _close(controller_torque('pd',q,qd,qr,vr,ar),fb)
    assert _close(controller_torque('gpd',q,qd,qr,vr,ar),fb+gravity_vector(q))
    ff=mass_matrix(qr)@ar+coriolis_vector(qr,vr)+gravity_vector(qr)
    assert _close(controller_torque('ffpd',q,qd,qr,vr,ar),ff+fb)
    ctc=mass_matrix(q)@(ar-10*ed-25*e)+coriolis_vector(q,qd)+gravity_vector(q)
    assert _close(controller_torque('ctc',q,qd,qr,vr,ar),ctc)
    sc=(1.15,1.15,0.85)
    assert _close(controller_torque('gpd',q,qd,qr,vr,ar,model_scale=sc),
                  fb+0.85*gravity_vector(q))
    assert _close(controller_torque('ffpd',q,qd,qr,vr,ar,model_scale=sc),
                  1.15*mass_matrix(qr)@ar+1.15*coriolis_vector(qr,vr)
                  +0.85*gravity_vector(qr)+fb)
    assert _close(controller_torque('ctc',q,qd,qr,vr,ar,model_scale=sc),
                  1.15*mass_matrix(q)@(ar-10*ed-25*e)
                  +1.15*coriolis_vector(q,qd)+0.85*gravity_vector(q))
    try:
        controller_torque('unknown',q,qd,qr,vr,ar)
        raise AssertionError('unknown kind must raise ValueError')
    except ValueError:
        pass


def test_saturation():
    cmd=np.array([60.,-12.]); before=cmd.copy()
    act,mask=saturate_torque(cmd,[45.,8.])
    assert _close(act,[45.,-8.]) and np.array_equal(mask,[True,True])
    assert _close(cmd,before), 'do not mutate input'
    act,mask=saturate_torque(cmd,np.inf)
    assert _close(act,cmd) and not np.any(mask)
    act,mask=saturate_torque(cmd,20.)
    assert _close(act,[20.,-12.]) and np.array_equal(mask,[True,False])


def test_comparison():
    res=run_comparison()
    expected={'ideal','mismatch','saturation','regulation'}
    assert set(res)==expected
    for scenario in expected:
        assert set(res[scenario])=={'pd','gpd','ffpd','ctc'}
        for kind,r in res[scenario].items():
            assert np.isfinite(r['rmse']) and np.isfinite(r['final_error'])
            assert np.isfinite(r['peak_actual_torque'])
            assert 0 <= r['saturation_fraction'] <= 1
            if scenario!='saturation':
                assert r['saturation_fraction']==0
            else:
                assert np.all(np.abs(r['tau_actual']) <= np.array([45.,8.])+1e-10)
    # Exact model + no saturation: CTC obeys decoupled error dynamics.
    r=res['ideal']['ctc']
    assert r['final_error'] < 0.005, 'ideal CTC should track well'
    assert res['regulation']['pd']['final_error'] > 0.01, 'PD has gravity bias'
    return res


def test_summary(res):
    summary=summarize(res)
    assert set(summary)==set(res)
    for scenario in res:
        assert set(summary[scenario])==set(res[scenario])
        for kind in res[scenario]:
            r=res[scenario][kind]
            expected=(r['rmse'],r['final_error'],r['peak_actual_torque'],
                      r['saturation_fraction'])
            assert _close(summary[scenario][kind],expected)
    return summary


if __name__=='__main__':
    print('[1] Four controller laws and model mismatch')
    test_controllers(); print('  PASS')
    print('[2] Total-torque saturation')
    test_saturation(); print('  PASS')
    print('[3] Four scenarios x four controllers')
    results=test_comparison(); print('  PASS')
    print('[4] Extract and print performance metrics')
    summary=test_summary(results)
    for scenario, cases in summary.items():
        print('  '+scenario)
        for kind, (rmse,final,peak,sat) in cases.items():
            print(f'    {kind:4s}  RMSE={rmse:.5f} rad  '
                  f'final={final:.5f} rad  peak_tau={peak:.2f} Nm  sat={sat:.1%}')
    print('  PASS\nDay 22: all tests passed.')
