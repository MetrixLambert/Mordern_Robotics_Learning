"""Day 9: planar floating-base WBC-QP exercise.

Dependencies: numpy, cvxpy. Run: python day09_wbc_exercise.py
Fill the TODOs in solve_wbc(). The provided model is a physically consistent
planar floating-base point-mass torso plus two massless leg links, with a
massless foot in stationary point contact. Joint torque limits apply.

补充：9个优化变量，7等式约束，7不等式约束
虽然优化变量中包含f_c，但是实际上最终用作控制的只有tau
只要满足约束，最终可以存在f_c满足该解，那么在控制输入为tau的情况下，其余变量/非变量之间的关系还是会
符合mpc的描述
"""
import numpy as np
import cvxpy as cp

G = 9.81
MASS = 10.0
INERTIA = 0.6
L1, L2 = 0.5, 0.5
MU = 0.6
TAU_MAX = 100.0

W_X = 100.0  # weights for wbc objective
W_Z = 100.0
W_THETA = 0.01


def model(q, v):
    """q=[x,z,theta,q1,q2], v=qdot. Mass concentrated in torso.

    Foot position = base position + links pointing downward; theta and
    joint angles rotate successive links. Coordinates are planar (x,z).
    """
    x, z, th, q1, q2 = q
    a1, a2 = th + q1, th + q1 + q2
    foot = np.array([x + L1*np.sin(a1) + L2*np.sin(a2),
                     z - L1*np.cos(a1) - L2*np.cos(a2)])
    J = np.zeros((2, 5))
    J[:, 0] = [1, 0]
    J[:, 1] = [0, 1]
    col1 = np.array([L1*np.cos(a1)+L2*np.cos(a2), L1*np.sin(a1)+L2*np.sin(a2)])
    col2 = np.array([L2*np.cos(a2), L2*np.sin(a2)])
    J[:, 2] = col1
    J[:, 3] = col1
    J[:, 4] = col2
    da1, da2 = v[2]+v[3], v[2]+v[3]+v[4]
    jdotv = np.array([-L1*np.sin(a1)*da1**2-L2*np.sin(a2)*da2**2,
                      L1*np.cos(a1)*da1**2+L2*np.cos(a2)*da2**2])
    M = np.diag([MASS, MASS, INERTIA, 0.0, 0.0])
    h = np.array([0.0, MASS*G, 0.0, 0.0, 0.0])
    return foot, M, h, J, jdotv


def solve_wbc(ax_cmd, q=None, v=None):
    if q is None:
        q = np.array([0.0, np.cos(0.4), 0.0, 0.4, -0.8])
    if v is None:
        v = np.zeros(5)
    foot, M, h, J, jdotv = model(q, v)
    a = cp.Variable(5)
    f = cp.Variable(2)  # [fx, fz]
    tau = cp.Variable(2)
    S_T = np.vstack([np.zeros((3, 2)), np.eye(2)])

    # TODO 1: complete the floating-base dynamics equality
    constraints = [
        M @ a + h == S_T @ tau + J.T @ f, # dynamics equality
        J @ a + jdotv == 0, # fixed-foot acceleration equality
        f[1] >= 0, # unilateral normal force
        f[0] <= MU * f[1], # planar friction cone
        f[0] >= -MU * f[1], # planar friction cone
        tau[0] <= TAU_MAX, # torque bounds
        tau[0] >= -TAU_MAX, # torque bounds
        tau[1] <= TAU_MAX, # torque bounds
        tau[1] >= -TAU_MAX, # torque bounds
    ]
    # COM coincides with torso origin in this simplified model.
    # TODO 5: construct weighted acceleration-tracking objective:
    # horizontal a[0] tracks ax_cmd, vertical a[1] tracks 0,
    # torso angular a[2] tracks 0. Add small force and torque penalties.
    objective = W_X * cp.sum_squares(a[0] - ax_cmd) + W_Z * cp.sum_squares(a[1]) + W_THETA * cp.sum_squares(a[2]) + 1e-5 * cp.sum_squares(f) + 1e-5 * cp.sum_squares(tau) # weighted acceleration-tracking objective
    if objective is None:
        raise NotImplementedError('Fill TODOs in solve_wbc() first')
    problem = cp.Problem(cp.Minimize(objective), constraints)
    problem.solve(solver=cp.CLARABEL)
    if problem.status not in ('optimal', 'optimal_inaccurate'):
        return {'status': problem.status, 'ax_cmd': ax_cmd}
    av, fv, tv = a.value, f.value, tau.value
    return {
        'status': problem.status, 'ax_cmd': ax_cmd,
        'a': av, 'f': fv, 'tau': tv, 'foot': foot,
        'dynamics_residual': np.linalg.norm(M @ av + h - S_T @ tv - J.T @ fv),
        'contact_residual': np.linalg.norm(J @ av + jdotv),
        'friction_margin': MU*fv[1] - abs(fv[0]),
    }


if __name__ == '__main__':
    for ax in [0.0, 1.0, 3.0, 8.0]:
        result = solve_wbc(ax)
        print('\ncommanded horizontal acceleration:', ax)
        for key, value in result.items():
            if key != 'ax_cmd':
                print(f'  {key}: {value}')
