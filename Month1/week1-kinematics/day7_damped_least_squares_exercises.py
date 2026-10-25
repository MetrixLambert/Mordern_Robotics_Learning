"""ROBO_Learn | Day 7 extension: Damped Least Squares (DLS) inverse kinematics.

Run: python day7_damped_least_squares_exercises.py
Dependency: numpy

Complete ONLY the TODOs in dls_step() and ik_body_dls().
All earlier Day 5-7 utilities are provided, so this file runs independently.
Twist order: [omega_x, omega_y, omega_z, v_x, v_y, v_z].
"""
import numpy as np


# ============== Provided: geometry and kinematics utilities ==============
def skew(w):
    x, y, z = np.asarray(w, dtype=float).reshape(3)
    return np.array([[0., -z, y], [z, 0., -x], [-y, x, 0.]])


def so3_exp(w, theta):
    W = skew(w)
    return np.eye(3) + np.sin(theta) * W + (1 - np.cos(theta)) * W @ W


def se3_exp(S, theta):
    S = np.asarray(S, dtype=float).reshape(6)
    w, v = S[:3], S[3:]
    T = np.eye(4)
    if np.linalg.norm(w) < 1e-12:
        T[:3, 3] = v * theta
        return T
    W = skew(w)
    T[:3, :3] = so3_exp(w, theta)
    G = theta * np.eye(3) + (1 - np.cos(theta)) * W + (theta - np.sin(theta)) * (W @ W)
    T[:3, 3] = G @ v
    return T


def trans_inv(T):
    R, p = T[:3, :3], T[:3, 3]
    result = np.eye(4)
    result[:3, :3] = R.T
    result[:3, 3] = -R.T @ p
    return result


def adjoint(T):
    R, p = T[:3, :3], T[:3, 3]
    A = np.zeros((6, 6))
    A[:3, :3] = R
    A[3:, 3:] = R
    A[3:, :3] = skew(p) @ R
    return A


def fk_body(M, B_list, q):
    T = np.array(M, dtype=float, copy=True)
    for i in range(len(q)):
        T = T @ se3_exp(B_list[:, i], q[i])
    return T


def jacobian_body(B_list, q):
    n = len(q)
    J = np.zeros((6, n))
    T = np.eye(4)
    for i in range(n - 1, -1, -1):
        J[:, i] = adjoint(T) @ B_list[:, i]
        T = T @ se3_exp(B_list[:, i], -q[i])
    return J


def se3_log(T):
    """Principal SE(3) log; sufficient for exercise poses (no exactly-pi rotation)."""
    R, p = T[:3, :3], T[:3, 3]
    theta = np.arccos(np.clip((np.trace(R) - 1) / 2, -1, 1))
    if theta < 1e-8:
        return np.r_[np.zeros(3), p]
    phi_hat = theta / (2 * np.sin(theta)) * (R - R.T)
    phi = np.array([phi_hat[2, 1], phi_hat[0, 2], phi_hat[1, 0]])
    c = 1 / theta**2 - 1 / (2 * theta * np.tan(theta / 2))
    J_inv = np.eye(3) - 0.5 * phi_hat + c * (phi_hat @ phi_hat)
    return np.r_[phi, J_inv @ p]


def example_robot():
    """Planar 2R robot, L1=2, L2=1; body screw axes are columns."""
    B_list = np.array([[0, 0], [0, 0], [1, 1], [0, 0], [3, 1], [0, 0]], dtype=float)
    M = np.eye(4)
    M[0, 3] = 3.0
    return B_list, M


# =================== Exercise 1: DLS update ===================
def dls_step(J, error, damping):
    """Solve min_dq ||J @ dq - error||^2 + damping^2 ||dq||^2.

    J: shape (m, n); error: shape (m,); damping: positive float.
    Returns dq: shape (n,).
    Hint: np.linalg.solve(A, b) is preferable to np.linalg.inv(A) @ b.
    """
    J = np.asarray(J, dtype=float)
    error = np.asarray(error, dtype=float)
    if damping <= 0:
        raise ValueError('damping must be positive for DLS')

    # TODO 1: form the n x n regularized normal-equation matrix
    A = J.T @ J + damping **2 * np.eye(J.shape[1]) 

    # TODO 2: form the n-dimensional right-hand side
    b = J.T @ error 

    # TODO 3: solve A @ dq = b (do not explicitly invert A)
    dq = np.linalg.solve(A, b)
    return dq


# =================== Exercise 2: Body IK with DLS ===================
def ik_body_dls(B_list, M, T_d, q0, damping=0.1,
                eomg=1e-4, ev=1e-4, max_iter=100):
    """Return (q, success, history).

    history entries: (angular_error_norm, linear_error_norm, joint_step_norm).
    Fixed damping is used for now; no line search or joint limits.
    """
    q = np.asarray(q0, dtype=float).copy()
    history = []
    for _ in range(max_iter):
        # TODO 4: compute current pose
        T = fk_body(M, B_list, q) 

        # TODO 5: compute body-frame relative transform
        T_err = np.linalg.inv(T) @ T_d 

        # TODO 6: compute six-dimensional pose error
        V_err = se3_log(T_err) 

        err_omega = np.linalg.norm(V_err[:3])
        err_v = np.linalg.norm(V_err[3:])
        if err_omega < eomg and err_v < ev:
            history.append((err_omega, err_v, 0.0))
            return q, True, history

        # TODO 7: compute Body Jacobian
        Jb = jacobian_body(B_list, q) 

        # TODO 8: compute damped joint update with dls_step()
        dq = dls_step(Jb, V_err, damping)

        history.append((err_omega, err_v, np.linalg.norm(dq)))

        # TODO 9: update joints coordinates
        q = q + dq 

    return q, False, history


# =================== Tests: do not modify ===================
def test_dls_step():
    print('\n[Exercise 1] DLS linear algebra tests')
    J = np.diag([1., 0.001])
    e = np.array([0., 1.])
    dq = dls_step(J, e, damping=0.1)
    expected = np.array([0., 0.001 / (0.001**2 + 0.1**2)])
    np.testing.assert_allclose(dq, expected, atol=1e-10)
    print('Near-singular example:')
    print('  pseudoinverse dq =', np.linalg.pinv(J) @ e)
    print('  DLS dq           =', dq)

    # Rectangular matrix test: 6D twist error and 2 joints.
    rng = np.random.default_rng(2026)
    J = rng.normal(size=(6, 2))
    e = rng.normal(size=6)
    lam = 0.2
    dq = dls_step(J, e, lam)
    # Equivalent dual form for positive damping.
    expected = J.T @ np.linalg.solve(J @ J.T + lam**2 * np.eye(6), e)
    np.testing.assert_allclose(dq, expected, atol=1e-10)
    print('PASS: square and rectangular Jacobians')


def test_ik_body_dls():
    print('\n[Exercise 2] Damped Body IK')
    B_list, M = example_robot()
    q_true = np.array([0.5, 0.7])
    T_d = fk_body(M, B_list, q_true)
    q0 = np.array([0.4, 0.6])
    q, success, history = ik_body_dls(B_list, M, T_d, q0, damping=0.1)
    print('  success:', success, '| iterations:', len(history))
    print('  target q:', q_true, '| recovered q:', q)
    print('  final errors:', history[-1][:2])
    assert success, 'DLS IK did not converge'
    np.testing.assert_allclose(fk_body(M, B_list, q), T_d, atol=2e-4)
    print('PASS: reachable target recovered')


def compare_near_singularity():
    """Optional exploration: near-singular POSITION Jacobian for planar 2R.

    This is a task-specific positional singularity, NOT a singularity of the
    full 6 x 2 Body Jacobian. The distinction is important.
    """
    print('\n[Optional] Near-singular planar 2R POSITION Jacobian')
    L1, L2 = 2., 1.
    q1, q2 = 0., 0.001
    Jp = np.array([
        [-L1*np.sin(q1)-L2*np.sin(q1+q2), -L2*np.sin(q1+q2)],
        [ L1*np.cos(q1)+L2*np.cos(q1+q2),  L2*np.cos(q1+q2)]
    ])
    desired_dp = np.array([0.01, 0.])
    print('  singular values:', np.linalg.svd(Jp, compute_uv=False))
    print('  pinv step:', np.linalg.pinv(Jp) @ desired_dp)
    for lam in (0.01, 0.1, 0.5):
        dq = dls_step(Jp, desired_dp, lam)
        print(f'  damping={lam:<4}: dq={dq}, achieved_linear_dp={Jp @ dq}')


if __name__ == '__main__':
    test_dls_step()
    test_ik_body_dls()
    compare_near_singularity()
