"""Day 13.2: Compass gait single-support dynamics (point-mass legs)."""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

l, mh, m, g = 1.0, 10.0, 1.0, 9.81
gamma = np.deg2rad(3.0)  # downhill toward +x: y_ground=-tan(gamma)*x

def dynamics(t, x):
    q1, q2, w1, w2 = x
    d = q1-q2
    mass = l*l*np.array([[mh+1.25*m, -0.5*m*np.cos(d)],
                          [-0.5*m*np.cos(d), 0.25*m]])
    coriolis = np.array([-0.5*m*l*l*np.sin(d)*w2*w2,
                          0.5*m*l*l*np.sin(d)*w1*w1])
    gravity = np.array([-(mh+1.5*m)*g*l*np.sin(q1),
                         0.5*m*g*l*np.sin(q2)])
    acc = np.linalg.solve(mass, -(coriolis+gravity))
    return [w1, w2, *acc]

def swing_foot(x):
    q1,q2=x[:2]
    return np.array([l*(np.sin(q1)-np.sin(q2)), l*(np.cos(q1)-np.cos(q2))])

def gap(x):
    px,py=swing_foot(x)
    return py+np.tan(gamma)*px

def touchdown(t,x):
    return gap(x)
touchdown.direction=-1
touchdown.terminal=True

def energy(x):
    q1,q2,w1,w2=x
    d=q1-q2
    kinetic=0.5*l*l*((mh+1.25*m)*w1*w1+0.25*m*w2*w2-m*np.cos(d)*w1*w2)
    potential=(mh+1.5*m)*g*l*np.cos(q1)-0.5*m*g*l*np.cos(q2)
    return kinetic+potential

# Front swing foot is initially above a slope descending toward +x.
x0=np.array([np.deg2rad(15),np.deg2rad(-15),0.8,-0.5])
if gap(x0)<=0:
    raise ValueError('Initial swing foot must be above the slope.')
sol=solve_ivp(dynamics,(0,1.5),x0,events=touchdown,max_step=0.002,rtol=1e-10,atol=1e-12,dense_output=True)
E=np.array([energy(x) for x in sol.y.T]); err=np.max(np.abs(E-E[0]))
print(f'Initial swing-foot clearance: {gap(x0):.6f} m')
print(f'Integration ended at t={sol.t[-1]:.6f} s')
print(f'Touchdown detected: {len(sol.t_events[0])>0}')
if len(sol.t_events[0]):
    print(f'Pre-impact state [q1,q2,w1,w2]: {sol.y_events[0][0]}')
print(f'Max mechanical-energy drift: {err:.3e} J')

fig,ax=plt.subplots(3,1,figsize=(10,9),sharex=True)
ax[0].plot(sol.t,np.rad2deg(sol.y[0]),label='stance q1');ax[0].plot(sol.t,np.rad2deg(sol.y[1]),label='swing q2');ax[0].set_ylabel('Angle (deg)');ax[0].legend()
ax[1].plot(sol.t,sol.y[2],label='stance w1');ax[1].plot(sol.t,sol.y[3],label='swing w2');ax[1].set_ylabel('Angular speed (rad/s)');ax[1].legend()
ax[2].plot(sol.t,E-E[0],label='E(t)-E(0)');ax[2].set_ylabel('Energy drift (J)');ax[2].set_xlabel('Time (s)');ax[2].legend()
for a in ax:a.grid(True,alpha=.3)
fig.suptitle('Compass gait: single support until touchdown')
fig.tight_layout();fig.savefig('day13_2_single_support.png',dpi=160)
print('Plot saved: day13_2_single_support.png')
