"""Day 9 extension: compare point contact and finite flat-foot contact WBC-QP.

Run: python day09_flat_foot_wbc.py
Requires: numpy, cvxpy, clarabel.
The flat foot is rigidly attached to the second link; it adds a contact moment
and a no-angular-acceleration contact constraint. Model is instantaneous and
uses massless legs; results are illustrative, not a dynamic walking simulation.
"""
import numpy as np
import cvxpy as cp
import matplotlib.pyplot as plt


G, MASS, INERTIA = 9.81, 10.0, 0.6
L1, L2 = 0.5, 0.5
MU, TAU_MAX = 0.6, 100.0
FOOT_HALF_LENGTH = 0.12
WX, WZ, WTHETA = 100.0, 100.0, 10.0
RF, RTAU = 1e-5, 1e-5


def model(q, v):
    x, z, th, q1, q2 = q
    alpha1, alpha2 = th + q1, th + q1 + q2
    foot = np.array([x + L1*np.sin(alpha1) + L2*np.sin(alpha2),
                     z - L1*np.cos(alpha1) - L2*np.cos(alpha2)])
    J = np.zeros((2, 5))
    J[:, 0] = [1, 0]
    J[:, 1] = [0, 1]
    col1 = [L1*np.cos(alpha1) + L2*np.cos(alpha2),
            L1*np.sin(alpha1) + L2*np.sin(alpha2)]
    col2 = [L2*np.cos(alpha2), L2*np.sin(alpha2)]
    J[:, 2] = col1
    J[:, 3] = col1
    J[:, 4] = col2
    da1, da2 = v[2]+v[3], v[2]+v[3]+v[4]
    jdotv = np.array([-L1*np.sin(alpha1)*da1**2-L2*np.sin(alpha2)*da2**2,
                      L1*np.cos(alpha1)*da1**2+L2*np.cos(alpha2)*da2**2])
    M = np.diag([MASS, MASS, INERTIA, 0.0, 0.0])
    h = np.array([0.0, MASS*G, 0.0, 0.0, 0.0])
    return foot, M, h, J, jdotv


def solve_wbc(ax_cmd, flat_foot=False, q=None, v=None):
    if q is None:
        # alpha2=0 so a rigid horizontal foot is compatible with the ground.
        # q = np.array([0.0, L1*np.cos(0.4)+L2, 0.0, 0.4, -0.4])
        q = np.array([0.0, L1*np.cos(0.4)+L2, 0.0, 0.1, -0.1])
    if v is None:
        v = np.zeros(5)
    foot, M, h, J, jdotv = model(q, v)
    if flat_foot:
        Jw = np.vstack([J, [0.0, 0.0, 1.0, 1.0, 1.0]])
        bias = np.r_[jdotv, 0.0]
        w = cp.Variable(3)  # fx, fz, contact moment about foot center
    else:
        Jw = J
        bias = jdotv
        w = cp.Variable(2)  # fx, fz
    a = cp.Variable(5)
    tau = cp.Variable(2)
    S_T = np.vstack([np.zeros((3, 2)), np.eye(2)])

    constraints = [
        M @ a + h == S_T @ tau + Jw.T @ w,
        Jw @ a + bias == 0,
        w[1] >= 0,
        w[0] <= MU*w[1],
        w[0] >= -MU*w[1],
        tau <= TAU_MAX,
        tau >= -TAU_MAX,
    ]
    if flat_foot:
        constraints += [w[2] <= FOOT_HALF_LENGTH*w[1],
                        w[2] >= -FOOT_HALF_LENGTH*w[1]]

    objective = (
        WX*cp.square(a[0]-ax_cmd)
        + WZ*cp.square(a[1])
        + WTHETA*cp.square(a[2])
        + RF*cp.sum_squares(w[:2])
        + RTAU*cp.sum_squares(tau)
    )
    prob = cp.Problem(cp.Minimize(objective), constraints)
    prob.solve(solver=cp.CLARABEL)
    if prob.status not in ('optimal', 'optimal_inaccurate'):
        return {'status': prob.status, 'ax_cmd': ax_cmd, 'flat_foot': flat_foot}
    av, wv, tv = a.value, w.value, tau.value
    result = {
        'status': prob.status,
        'ax_cmd': ax_cmd,
        'flat_foot': flat_foot,
        'a': av,
        'wrench': wv,
        'tau': tv,
        'objective': prob.value,
        'dynamics_residual': np.linalg.norm(M@av+h-S_T@tv-Jw.T@wv),
        'contact_residual': np.linalg.norm(Jw@av+bias),
        'friction_margin': MU*wv[1]-abs(wv[0]),
        'foot_position': foot,
    }
    if flat_foot:
        result['cop_relative_x'] = -wv[2]/wv[1] if wv[1] > 1e-8 else np.nan
        result['cop_margin'] = FOOT_HALF_LENGTH*wv[1]-abs(wv[2])
        result['foot_orientation'] = q[2]+q[3]+q[4]
    return result


def show(result):
    name = 'FLAT FOOT' if result['flat_foot'] else 'POINT CONTACT'
    print(f"\n{name}, commanded ax={result['ax_cmd']:.2f}, status={result['status']}")
    if 'a' not in result:
        return
    print('  a     =', np.round(result['a'], 5))
    print('  wrench=', np.round(result['wrench'], 5))
    print('  tau   =', np.round(result['tau'], 5))
    print('  objective =', round(result['objective'], 6))
    print('  dynamics residual =', result['dynamics_residual'])
    print('  contact residual  =', result['contact_residual'])
    print('  friction margin   =', result['friction_margin'])
    if result['flat_foot']:
        print('  CoP relative x (m)=', result['cop_relative_x'])
        print('  CoP margin (N m)  =', result['cop_margin'])


if __name__ == '__main__':
    # for ax in [0.0, 1.0, 3.0, 8.0]:
    #     show(solve_wbc(ax, flat_foot=False))
    #     show(solve_wbc(ax, flat_foot=True))

    # plot the results
    ax_actual_list = [] 
    atheta_actual_list = []
    cop_actual_list = []
    ax_cmd_list = np.linspace(0.0, 1.2, 25) 

    for ax_cmd in ax_cmd_list:
        result = solve_wbc(ax_cmd, flat_foot=True)
        ax_actual_list.append(result['a'][0])
        atheta_actual_list.append(result['a'][2])
        cop_actual_list.append(result['cop_relative_x'])

    # plot the results
    plt.figure(figsize=(10, 5))
    plt.subplot(1, 3, 1)
    plt.plot(ax_cmd_list, ax_actual_list, label='ax_actual')
    plt.subplot(1, 3, 2)
    plt.plot(ax_cmd_list, atheta_actual_list, label='atheta_actual')
    plt.subplot(1, 3, 3)
    plt.plot(ax_cmd_list, cop_actual_list, label='cop_actual')
    plt.legend()
    plt.tight_layout()
    plt.show()

