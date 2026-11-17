"""Day 14.1: Compass gait, single support -> plastic heel strike -> leg relabel.
Point masses: hip Mh, each leg midpoint m. Absolute angles from upward vertical,
positive toward +x. Ground is y=-x*tan(gamma).
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

l, Mh, m, g = 1.0, 10.0, 1.0, 9.81
gamma = np.deg2rad(3.0)


def mass_reduced(q):
    a, b = q
    return l*l*np.array([[Mh+1.25*m, -.5*m*np.cos(a-b)],
                         [-.5*m*np.cos(a-b), .25*m]])


def velocity_terms(q, v):
    a, b = q
    da, db = v
    s = np.sin(a-b)
    return np.array([-.5*m*l*l*s*db*db, .5*m*l*l*s*da*da])


def gravity_terms(q):
    a, b = q
    return np.array([-(Mh+1.5*m)*g*l*np.sin(a), .5*m*g*l*np.sin(b)])


def ode(t, state):
    q, v = state[:2], state[2:]
    acc = np.linalg.solve(mass_reduced(q),
                          -velocity_terms(q, v)-gravity_terms(q))
    return np.r_[v, acc]


def swing_position(q):
    a, b = q
    return l*np.array([np.sin(a)-np.sin(b), np.cos(a)-np.cos(b)])


def swing_jacobian(q):
    a, b = q
    return l*np.array([[np.cos(a), -np.cos(b)],
                       [-np.sin(a), np.sin(b)]])


def foot_gap(q):
    px, py = swing_position(q)
    return py + px*np.tan(gamma)


def touchdown(t, state):
    return foot_gap(state[:2])

touchdown.terminal = True
touchdown.direction = -1


def floating_jacobians(q):
    a, b = q
    J1 = np.array([[1., 0., -l*np.cos(a), 0.],
                   [0., 1.,  l*np.sin(a), 0.]])
    J2 = np.array([[1., 0., 0., -l*np.cos(b)],
                   [0., 1., 0.,  l*np.sin(b)]])
    return J1, J2


def floating_mass(q):
    """Sum m_i * J_COM_i.T @ J_COM_i; no rod rotational inertia."""
    a, b = q
    Jh = np.array([[1.,0.,0.,0.],[0.,1.,0.,0.]])
    # stance COM: hip - (l/2)*[sin(a),cos(a)]
    Jc1 = np.array([[1.,0.,-.5*l*np.cos(a),0.],
                    [0.,1., .5*l*np.sin(a),0.]])
    # swing COM: hip - (l/2)*[sin(b),cos(b)]
    Jc2 = np.array([[1.,0.,0.,-.5*l*np.cos(b)],
                    [0.,1.,0., .5*l*np.sin(b)]])
    return Mh*(Jh.T@Jh) + m*(Jc1.T@Jc1) + m*(Jc2.T@Jc2)


def impact(q, v):
    """Old stance released; new foot fully sticks. Returns q_new,v_new,impulse, diagnostics."""
    J1, J2 = floating_jacobians(q)
    # Preimpact hip velocity from fixed old stance foot
    a, b = q
    vhip = l*np.array([np.cos(a)*v[0], -np.sin(a)*v[0]])
    vbar_minus = np.r_[vhip, v]
    Mbar = floating_mass(q)
    KKT = np.block([[Mbar, -J2.T], [J2, np.zeros((2,2))]])
    rhs = np.r_[Mbar@vbar_minus, np.zeros(2)]
    sol = np.linalg.solve(KKT, rhs)
    vbar_plus, impulse = sol[:4], sol[4:]
    # New leg1 is old leg2; new leg2 is old leg1
    q_new, v_new = q[::-1].copy(), vbar_plus[2:][::-1].copy()
    Tminus = .5*vbar_minus@Mbar@vbar_minus
    Tplus = .5*vbar_plus@Mbar@vbar_plus
    return q_new, v_new, impulse, {
        'old_foot_speed_before':np.linalg.norm(J1@vbar_minus),
        'new_foot_speed_before':np.linalg.norm(J2@vbar_minus),
        'new_foot_speed_after':np.linalg.norm(J2@vbar_plus),
        'kinetic_before':Tminus, 'kinetic_after':Tplus,
        'impulse_balance_error':np.linalg.norm(Mbar@(vbar_plus-vbar_minus)-J2.T@impulse),
        'vbar_before':vbar_minus, 'vbar_after':vbar_plus,
        'old_foot_speed_after':np.linalg.norm(J1@vbar_plus),
    }


def energy(q, v):
    a,b=q
    T=.5*v@mass_reduced(q)@v
    V=(Mh+1.5*m)*g*l*np.cos(a)-.5*m*g*l*np.cos(b)
    return T+V


def main():
    # Day 13 single-support initial condition, with front swing foot above slope.
    x0=np.array([np.deg2rad(15.),np.deg2rad(-15.),.8,-.5])
    print(f'Initial foot gap: {foot_gap(x0[:2]):.8f} m')
    sol=solve_ivp(ode,(0.,2.),x0,events=touchdown,max_step=.002,
                  rtol=1e-10,atol=1e-12,dense_output=True)
    if not len(sol.t_events[0]):
        raise RuntimeError('No touchdown found; adjust initial condition.')
    ti=float(sol.t_events[0][0]); xm=sol.y_events[0][0]
    qm,vm=xm[:2],xm[2:]
    qnew,vnew,lam,d=impact(qm,vm)
    print(f'Touchdown: t={ti:.8f} s, gap={foot_gap(qm):+.3e} m')
    print('Preimpact q [deg]:',np.round(np.rad2deg(qm),6))
    print('Preimpact qdot:',np.round(vm,6))
    print('Impulse [x,y] (N*s):',np.round(lam,8))
    print('Relabeled postimpact q [deg]:',np.round(np.rad2deg(qnew),6))
    print('Relabeled postimpact qdot:',np.round(vnew,8))
    for key in ('old_foot_speed_before','new_foot_speed_before',
                'new_foot_speed_after','old_foot_speed_after',
                'impulse_balance_error','kinetic_before','kinetic_after'):
        print(f'{key}: {d[key]:.10g}')
    print('Kinetic energy loss:',d['kinetic_before']-d['kinetic_after'],'J')
    print('Single-support energy drift:',np.max(np.abs(
        np.array([energy(sol.y[:2,i],sol.y[2:,i]) for i in range(sol.t.size)])
        -energy(x0[:2],x0[2:]))),'J')
    if lam[1]<0:
        print('WARNING: negative vertical impulse; assumed impact is not unilateral-feasible.')
    if lam[1]>0:
        print('Required impact friction ratio |Lambda_x|/Lambda_y:',abs(lam[0])/lam[1])
    print('NOTE: check old-foot liftoff and collision feasibility before multi-step simulation.')
    fig,axs=plt.subplots(3,1,figsize=(9,8),sharex=True)
    axs[0].plot(sol.t,np.rad2deg(sol.y[0]),label='stance q1')
    axs[0].plot(sol.t,np.rad2deg(sol.y[1]),label='swing q2')
    axs[0].set_ylabel('Angle (deg)'); axs[0].legend()
    axs[1].plot(sol.t,sol.y[2],label='stance qdot1')
    axs[1].plot(sol.t,sol.y[3],label='swing qdot2')
    axs[1].set_ylabel('Angular speed (rad/s)'); axs[1].legend()
    axs[2].plot(sol.t,[foot_gap(sol.y[:2,i]) for i in range(sol.t.size)])
    axs[2].axhline(0,color='k',ls='--',lw=.8)
    axs[2].set_ylabel('Swing-foot gap (m)'); axs[2].set_xlabel('Time (s)')
    for ax in axs: ax.grid(alpha=.3); ax.axvline(ti,color='r',ls=':',alpha=.7)
    fig.suptitle('Day 14.1: Single support to first touchdown')
    fig.tight_layout(); fig.savefig('day14_1_impact.png',dpi=160)
    print('Saved: day14_1_impact.png')

if __name__=='__main__': main()
