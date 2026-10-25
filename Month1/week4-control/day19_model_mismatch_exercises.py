"""Day 19 | Model mismatch: regulation and sinusoidal tracking (NumPy only).

Fill TODO 1-4 and run:
    python day19_model_mismatch_exercises.py

Single-DOF inertial joint: m*qdd = tau, with m>0.
Controller: tau = m_hat * (qdd_des - kd*ed - kp*e).
Error: e=q-q_des.  No friction, gravity, saturation or noise.
This isolates inertia mismatch before returning to the 2R robot.
"""
from __future__ import annotations
import numpy as np


# TODO 1 ---------------------------------------------------------------
def computed_torque_1d(q, qd, q_des, qd_des, qdd_des, m_hat, kp, kd):
    """Return scalar CTC torque using the estimated inertia m_hat.

    e=q-q_des; ed=qd-qd_des
    tau=m_hat*(qdd_des-kd*ed-kp*e)
    """
    e = q - q_des
    ed = qd - qd_des 
    tau = m_hat * (qdd_des - kd * ed - kp * e) 
    return tau 


def reference(t, kind, amplitude=0.5, omega=3.0, target=0.7):
    """Provided reference: either fixed regulation or persistent sinusoid."""
    if kind == 'regulation':
        return float(target), 0.0, 0.0
    if kind == 'sinusoid':
        return (amplitude*np.sin(omega*t),
                amplitude*omega*np.cos(omega*t),
                -amplitude*omega**2*np.sin(omega*t))
    raise ValueError('kind must be regulation or sinusoid')


def simulate(m, m_hat, kp, kd, kind, *, dt=0.001, duration=12.0,
             initial_error=0.2, initial_velocity_error=0.0,
             amplitude=0.5, omega=3.0, target=0.7):
    """Provided semi-implicit Euler integration of true dynamics m*qdd=tau."""
    assert m>0 and m_hat>0 and kp>0 and kd>0 and dt>0
    n=int(round(duration/dt))
    ts=np.arange(n+1)*dt
    q=np.zeros(n+1); qd=np.zeros(n+1); tau=np.zeros(n)
    qref=np.zeros(n+1); qdref=np.zeros(n+1)
    qref[0],qdref[0],_=reference(0,kind,amplitude,omega,target)
    q[0]=qref[0]+initial_error
    qd[0]=qdref[0]+initial_velocity_error
    for k,t in enumerate(ts):
        qr,vr,ar=reference(t,kind,amplitude,omega,target)
        qref[k],qdref[k]=qr,vr
        if k==n: break
        tau[k]=computed_torque_1d(q[k],qd[k],qr,vr,ar,m_hat,kp,kd)
        qdd=tau[k]/m
        qd[k+1]=qd[k]+dt*qdd
        q[k+1]=q[k]+dt*qd[k+1]
    return dict(t=ts,q=q,qd=qd,tau=tau,q_des=qref,
                error=q-qref,velocity_error=qd-qdref)


# TODO 2 ---------------------------------------------------------------
def effective_second_order(m, m_hat, kp, kd):
    """Return (omega_eff, zeta_eff) for regulation error dynamics:

    e_dd + alpha*kd*e_d + alpha*kp*e = 0, alpha=m_hat/m.
    omega_eff=sqrt(alpha*kp)
    zeta_eff=alpha*kd/(2*omega_eff)
    """
    alpha = m_hat / m 
    omega_eff = np.sqrt(alpha * kp) 
    zeta_eff = alpha * kd / (2 * omega_eff) 
    return omega_eff, zeta_eff 


# TODO 3 ---------------------------------------------------------------
def predicted_sinusoidal_error_amplitude(m, m_hat, kp, kd, amplitude, omega):
    """Return nonnegative theoretical steady-state error amplitude.

    e_dd+alpha*kd*e_d+alpha*kp*e=(alpha-1)*qdd_des
    q_des=A*sin(omega*t)
    E_amp=abs(alpha-1)*A*omega**2 /
          sqrt((alpha*kp-omega**2)**2+(alpha*kd*omega)**2)
    """
    alpha = m_hat / m 
    E_amp = abs(alpha-1) * amplitude * omega**2 / np.sqrt((alpha*kp-omega**2)**2+(alpha*kd*omega)**2) 
    return E_amp 


# TODO 4 ---------------------------------------------------------------
def compare_mismatch_cases(m=1.0, kp=25.0, kd=10.0):
    """Return a dictionary with exactly these keys:

      'reg_under':  regulation, m_hat=0.7*m
      'reg_exact':  regulation, m_hat=1.0*m
      'reg_over':   regulation, m_hat=1.5*m
      'sin_under':  sinusoid,   m_hat=0.7*m
      'sin_exact':  sinusoid,   m_hat=1.0*m
      'sin_over':   sinusoid,   m_hat=1.5*m

    All runs: dt=0.001, duration=12.0, initial_error=0.2,
    initial_velocity_error=0.0; sinusoid A=0.5, omega=3.0.
    Use provided simulate(); return its output dictionaries.
    """
    results = {} 
    results['reg_under'] = simulate(m, 0.7*m, kp, kd, 'regulation')
    results['reg_exact'] = simulate(m, 1.0*m, kp, kd, 'regulation')
    results['reg_over'] = simulate(m, 1.5*m, kp, kd, 'regulation')
    results['sin_under'] = simulate(m, 0.7*m, kp, kd, 'sinusoid')
    results['sin_exact'] = simulate(m, 1.0*m, kp, kd, 'sinusoid')
    results['sin_over'] = simulate(m, 1.5*m, kp, kd, 'sinusoid')
    return results 


def steady_state_amplitude(result, *, omega=3.0, window_start=8.0):
    """Provided: sinusoidal least-squares amplitude over final time window.

    Fits e(t) = a*sin(omega*t)+b*cos(omega*t)+c to remove DC/numeric bias.
    """
    t=result['t']; e=result['error']; mask=t>=window_start
    x=np.column_stack((np.sin(omega*t[mask]),np.cos(omega*t[mask]),
                       np.ones(mask.sum())))
    coeff=np.linalg.lstsq(x,e[mask],rcond=None)[0]
    return float(np.hypot(coeff[0],coeff[1]))


def run_tests():
    print('[1] CTC torque and signs')
    tau=computed_torque_1d(0.4,-0.2,0.1,0.3,1.5,2.0,25.0,10.0)
    expected=2.0*(1.5-10.0*(-0.5)-25.0*0.3)
    assert np.isclose(tau,expected,atol=1e-12)
    print('  PASS')

    print('[2] Effective frequency and damping')
    for alpha in (0.7,1.0,1.5):
        wn,zeta=effective_second_order(1.0,alpha,25.0,10.0)
        assert np.isclose(wn,5*np.sqrt(alpha),rtol=1e-10)
        assert np.isclose(zeta,np.sqrt(alpha),rtol=1e-10)
        print(f'  alpha={alpha:.1f}  omega_eff={wn:.4f}  zeta_eff={zeta:.4f}')
    print('  PASS')

    print('[3] Sinusoidal steady-state amplitude formula')
    assert np.isclose(predicted_sinusoidal_error_amplitude(1,1,25,10,.5,3),0)
    for alpha in (.7,1.5):
        amp=predicted_sinusoidal_error_amplitude(1,alpha,25,10,.5,3)
        expected=abs(alpha-1)*.5*9/np.hypot(alpha*25-9,alpha*10*3)
        assert np.isclose(amp,expected,rtol=1e-10)
    print('  PASS')

    print('[4] Six closed-loop experiments')
    results=compare_mismatch_cases()
    assert set(results)=={'reg_under','reg_exact','reg_over','sin_under','sin_exact','sin_over'}
    for key in ('reg_under','reg_exact','reg_over'):
        r=results[key]
        assert np.isclose(r['t'][-1],12.0)
        assert abs(r['error'][-1])<2e-4, (key,r['error'][-1])
        print(f'  {key:10s} final_error={r["error"][-1]:+.3e} rad')
    for key,alpha in [('sin_under',.7),('sin_exact',1.),('sin_over',1.5)]:
        r=results[key]
        measured=steady_state_amplitude(r)
        predicted=predicted_sinusoidal_error_amplitude(1,alpha,25,10,.5,3)
        print(f'  {key:10s} amplitude measured={measured:.6f}, predicted={predicted:.6f} rad')
        assert abs(measured-predicted)<0.0015, (key,measured,predicted)
    print('  PASS')
    print('\nDay 19: all tests passed.')


if __name__=='__main__':
    run_tests()
