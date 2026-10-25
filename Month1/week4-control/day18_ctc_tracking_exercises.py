"""Day 18 | FF+PD versus Computed Torque Tracking (NumPy only).

Fill TODO 1-4, then run:
    python day18_ctc_tracking_exercises.py

2R vertical-plane arm, same parameters as Days 16-17.
Error convention: e = q - q_desired.
IMPORTANT: FF+PD gains have torque/rad units; CTC gains have 1/s^2 and 1/s.
The comparison uses separately tuned gains, not numerically identical ones.
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


def mass_matrix(robot, q):
    q2 = float(q[1])
    a = robot.I1 + robot.I2 + robot.m1*robot.r1**2 + robot.m2*(robot.L1**2+robot.r2**2)
    b = robot.m2*robot.L1*robot.r2
    d = robot.I2+robot.m2*robot.r2**2
    return np.array([[a+2*b*np.cos(q2), d+b*np.cos(q2)],
                     [d+b*np.cos(q2), d]], dtype=float)


def coriolis_vector(robot, q, qd):
    h = robot.m2*robot.L1*robot.r2*np.sin(float(q[1]))
    v1,v2 = np.asarray(qd,dtype=float)
    return np.array([-h*(2*v1*v2+v2*v2), h*v1*v1])


def gravity_vector(robot,q):
    q1,q2=np.asarray(q,dtype=float)
    a=(robot.m1*robot.r1+robot.m2*robot.L1)*np.cos(q1)
    b=robot.m2*robot.r2*np.cos(q1+q2)
    return robot.gravity*np.array([a+b,b])


def forward_dynamics(robot,q,qd,tau):
    return np.linalg.solve(mass_matrix(robot,q),
                           tau-coriolis_vector(robot,q,qd)-gravity_vector(robot,q))


# TODO 1 ---------------------------------------------------------------
def quintic_trajectory(t, q0, qf, T):
    """Return (q_des, qd_des, qdd_des), each shape (2,).

    s=clip(t/T,0,1)
    h=10*s**3-15*s**4+6*s**5
    dh/dt=(30*s**2-60*s**3+30*s**4)/T
    d2h/dt2=(60*s-180*s**2+120*s**3)/T**2
    Outside [0,T], return zero desired velocity and acceleration.
    """
    if t < 0: 
        return q0, np.zeros(2), np.zeros(2)
    if t > T:
        return qf, np.zeros(2), np.zeros(2)
    s = np.clip(t/T, 0, 1)
    h = 10*s**3 - 15*s**4 + 6*s**5
    dh_dt = (30*s**2 - 60*s**3 + 30*s**4) / T
    d2h_dt2 = (60*s - 180*s**2 + 120*s**3) / T**2
    return q0 + (qf - q0) * h, (qf - q0) * dh_dt, (qf - q0) * d2h_dt2


# TODO 2 ---------------------------------------------------------------
def control_ff_pd(robot, q, qd, q_des, qd_des, qdd_des, Kp_tau, Kd_tau):
    """Desired-state inverse-dynamics feedforward + torque-level PD:

    e=q-q_des; ed=qd-qd_des
    tau=M(q_des)@qdd_des+c(q_des,qd_des)+g(q_des)-Kp_tau@e-Kd_tau@ed
    """
    e = q - q_des 
    ed = qd - qd_des 
    tau = mass_matrix(robot, q_des) @ qdd_des + coriolis_vector(robot, q_des, qd_des) + gravity_vector(robot, q_des) - Kp_tau @ e - Kd_tau @ ed 
    return tau 


# TODO 3 ---------------------------------------------------------------
def control_ctc(robot, q, qd, q_des, qd_des, qdd_des, Kp_acc, Kd_acc):
    """Actual-state computed torque control (acceleration-level PD):

    e=q-q_des; ed=qd-qd_des
    v=qdd_des-Kd_acc@ed-Kp_acc@e
    tau=M(q)@v+c(q,qd)+g(q)
    """
    e = q - q_des 
    ed = qd - qd_des 
    v = qdd_des - Kd_acc @ ed - Kp_acc @ e 
    tau = mass_matrix(robot, q) @ v + coriolis_vector(robot, q, qd) + gravity_vector(robot, q) 
    return tau 


def simulate_tracking(robot, controller, q_start, qd_start, q_traj_start,
                      q_traj_final, T, Kp, Kd, dt=0.002, duration=4.0):
    """Provided: semi-implicit Euler simulation of the true robot."""
    n=int(round(duration/dt))
    assert dt>0 and T>0 and duration>=T and np.isclose(n*dt,duration)
    ts=np.arange(n+1)*dt
    q=np.zeros((n+1,2)); qd=np.zeros((n+1,2)); tau=np.zeros((n,2))
    q_des=np.zeros_like(q); qd_des=np.zeros_like(q); qdd_des=np.zeros_like(q)
    q[0]=q_start; qd[0]=qd_start
    for k,t in enumerate(ts):
        q_des[k],qd_des[k],qdd_des[k]=quintic_trajectory(t,q_traj_start,q_traj_final,T)
        if k==n: break
        tau[k]=controller(robot,q[k],qd[k],q_des[k],qd_des[k],qdd_des[k],Kp,Kd)
        qdd=forward_dynamics(robot,q[k],qd[k],tau[k])
        qd[k+1]=qd[k]+dt*qdd
        q[k+1]=q[k]+dt*qd[k+1]
    return dict(t=ts,q=q,qd=qd,tau=tau,q_des=q_des,qd_des=qd_des,
                error=q-q_des, velocity_error=qd-qd_des)


# TODO 4 ---------------------------------------------------------------
def compare_controllers(robot, q_traj_start, q_traj_final, T,
                        initial_position_offset, initial_velocity_offset,
                        Kp_tau, Kd_tau, Kp_acc, Kd_acc):
    """Return dict with keys 'ff_pd' and 'ctc'.

    Both use same actual initial state:
      q_start = q_traj_start + initial_position_offset
      qd_start = initial_velocity_offset  # desired initial velocity is zero

    Simulate with dt=0.002, duration=4.0.
    'ff_pd' uses Kp_tau/Kd_tau; 'ctc' uses Kp_acc/Kd_acc.
    """
    q_start = q_traj_start + initial_position_offset 
    qd_start = initial_velocity_offset 
    dt = 0.002 
    duration = 4.0 
    results = {} 
    results['ff_pd'] = simulate_tracking(robot, control_ff_pd, q_start, qd_start, q_traj_start, q_traj_final, T, Kp_tau, Kd_tau, dt, duration)
    results['ctc'] = simulate_tracking(robot, control_ctc, q_start, qd_start, q_traj_start, q_traj_final, T, Kp_acc, Kd_acc, dt, duration)
    return results 


def metrics(result):
    """Provided: metrics computed over the entire trajectory."""
    err=result['error']
    return dict(rmse=float(np.sqrt(np.mean(np.sum(err**2,axis=1)))),
                final_error=float(np.linalg.norm(err[-1])),
                max_torque=float(np.max(np.abs(result['tau']))))


def run_tests():
    robot=Planar2R()
    q0=np.array([-0.4,0.2]); qf=np.array([0.8,-0.5]); T=2.0
    # CTC: wn=5 rad/s, zeta=1 -> kp=25 s^-2, kd=10 s^-1.
    Kp_acc=np.diag([25.,25.]); Kd_acc=np.diag([10.,10.])
    # FF+PD: independent torque-level gains (N m/rad, N m s/rad).
    Kp_tau=np.diag([100.,35.]); Kd_tau=np.diag([30.,12.])

    print('[1] Quintic trajectory boundary conditions')
    for t,q_expected in [(0.,q0),(T,qf),(-0.2,q0),(T+0.2,qf)]:
        q,v,a=quintic_trajectory(t,q0,qf,T)
        np.testing.assert_allclose(q,q_expected,atol=1e-12)
        if t<=0 or t>=T:
            np.testing.assert_allclose(v,np.zeros(2),atol=1e-12)
            np.testing.assert_allclose(a,np.zeros(2),atol=1e-12)
    q,v,a=quintic_trajectory(T/2,q0,qf,T)
    np.testing.assert_allclose(q,(q0+qf)/2,atol=1e-12)
    np.testing.assert_allclose(v,1.875*(qf-q0)/T,atol=1e-12)
    np.testing.assert_allclose(a,np.zeros(2),atol=1e-12)
    print('  PASS')

    print('[2] FF+PD control torque')
    q=np.array([-0.2,0.1]); qd=np.array([0.15,-0.05]);
    q_des,v_des,a_des=quintic_trajectory(0.7,q0,qf,T)
    expected=(mass_matrix(robot,q_des)@a_des+coriolis_vector(robot,q_des,v_des)
              +gravity_vector(robot,q_des)-Kp_tau@(q-q_des)-Kd_tau@(qd-v_des))
    np.testing.assert_allclose(control_ff_pd(robot,q,qd,q_des,v_des,a_des,Kp_tau,Kd_tau),expected,atol=1e-11)
    print('  PASS')

    print('[3] CTC control torque and exact acceleration linearization')
    v=a_des-Kd_acc@(qd-v_des)-Kp_acc@(q-q_des)
    expected=mass_matrix(robot,q)@v+coriolis_vector(robot,q,qd)+gravity_vector(robot,q)
    tau=control_ctc(robot,q,qd,q_des,v_des,a_des,Kp_acc,Kd_acc)
    np.testing.assert_allclose(tau,expected,atol=1e-11)
    np.testing.assert_allclose(forward_dynamics(robot,q,qd,tau),v,atol=1e-11)
    print('  PASS')

    print('[4] Tracking with nonzero initial position error')
    initial_offset=np.array([0.2,-0.15]); initial_velocity=np.zeros(2)
    results=compare_controllers(robot,q0,qf,T,initial_offset,initial_velocity,
                                Kp_tau,Kd_tau,Kp_acc,Kd_acc)
    assert set(results)=={'ff_pd','ctc'}
    for name,res in results.items():
        assert res['q'].shape==(2001,2)
        np.testing.assert_allclose(res['error'][0],initial_offset,atol=1e-12)
        assert np.all(np.isfinite(res['q'])) and np.all(np.isfinite(res['tau']))
        print(f"  {name:6s}: {metrics(res)}")
    assert metrics(results['ctc'])['final_error']<1e-3
    print('  PASS')

    print('[5] CTC trajectory error against analytical critical-damping solution')
    r=results['ctc']; t=r['t'][:,None]
    # With e_dot(0)=0: e(t)=e(0)*(1+5t)*exp(-5t).
    e_theory=initial_offset[None,:]*(1+5*t)*np.exp(-5*t)
    max_deviation=float(np.max(np.abs(r['error']-e_theory)))
    print(f'  max absolute deviation = {max_deviation:.6g} rad')
    assert max_deviation<0.015  # integration and desired-trajectory discretization error
    print('  PASS\n\nDay 18: all tests passed.')
    return results


if __name__=='__main__':
    run_tests()
