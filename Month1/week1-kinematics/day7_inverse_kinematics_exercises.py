"""ROBO_Learn — Day 7: SE(3) logarithm and Body inverse kinematics.

Run:
    python day7_inverse_kinematics_exercises.py

Dependencies: numpy only.
Complete the TODOs in se3_log() and ik_body(). All other functions are provided.
Convention: 6D twists are ordered [omega_x, omega_y, omega_z, v_x, v_y, v_z].
"""

import numpy as np


# ===================== Provided utilities =====================

def skew(w):
    w = np.asarray(w, dtype=float).reshape(3)
    return np.array([[0., -w[2], w[1]], [w[2], 0., -w[0]], [-w[1], w[0], 0.]])


def so3_exp(w, theta):
    """Rodrigues formula; w must be a unit axis (unless theta is zero)."""
    w = np.asarray(w, dtype=float).reshape(3)
    W = skew(w)
    return np.eye(3) + np.sin(theta) * W + (1 - np.cos(theta)) * (W @ W)


def se3_exp(S, theta):
    """SE(3) exponential for a unit revolute screw or unit prismatic screw."""
    S = np.asarray(S, dtype=float).reshape(6)
    w, v = S[:3], S[3:]
    T = np.eye(4)
    if np.linalg.norm(w) < 1e-12:
        T[:3, 3] = v * theta
        return T
    W = skew(w)
    R = so3_exp(w, theta)
    G = np.eye(3) * theta + (1 - np.cos(theta)) * W + (theta - np.sin(theta)) * (W @ W)
    T[:3, :3] = R
    T[:3, 3] = G @ v
    return T


def trans_inv(T):
    R, p = T[:3, :3], T[:3, 3]
    T_inv = np.eye(4)
    T_inv[:3, :3] = R.T
    T_inv[:3, 3] = -R.T @ p
    return T_inv


def adjoint(T):
    R, p = T[:3, :3], T[:3, 3]
    Ad = np.zeros((6, 6))
    Ad[:3, :3] = R
    Ad[3:, 3:] = R
    Ad[3:, :3] = skew(p) @ R
    return Ad


def fk_body(M, B_list, q):
    """T(q) = M exp([B1]q1) ... exp([Bn]qn)."""
    T = np.array(M, dtype=float, copy=True)
    for i in range(len(q)):
        T = T @ se3_exp(B_list[:, i], q[i])
    return T


def jacobian_body(B_list, q):
    n = len(q)
    Jb = np.zeros((6, n))
    T = np.eye(4)
    for i in range(n - 1, -1, -1):
        Jb[:, i] = adjoint(T) @ B_list[:, i]
        T = T @ se3_exp(B_list[:, i], -q[i])
    return Jb


def vec_to_se3(V):
    V = np.asarray(V, dtype=float).reshape(6)
    result = np.zeros((4, 4))
    result[:3, :3] = skew(V[:3])
    result[:3, 3] = V[3:]
    return result


# ===================== Exercise 1: implement this =====================

def se3_log(T):
    """Return the 6D logarithm [phi, rho] of a homogeneous transform T.

    This exercise handles angles in [0, pi), excluding the exactly-pi case.
    For small angles, the provided tests only require exact pure translation.
    """
    R = T[:3, :3]
    p = T[:3, 3]

    # TODO 1: theta = arccos(clip((trace(R)-1)/2, -1, 1))
    theta = np.arccos(np.clip((np.trace(R)-1)/2, -1, 1))

    # TODO 2: when theta is near zero, return [0,0,0, p_x,p_y,p_z].
    if theta < 1e-8:
        return np.array([0,0,0, p[0], p[1], p[2]])

    # TODO 3: [phi]_x = theta/(2*sin(theta)) * (R - R.T)
    # Extract phi = [phi_hat[2,1], phi_hat[0,2], phi_hat[1,0]].
    phi_hat = theta / (2 * np.sin(theta)) * (R - R.T)

    # TODO 4: construct [phi]_x
    phi = np.array([phi_hat[2,1], phi_hat[0,2], phi_hat[1,0]])

    # TODO 5: J_l(phi)^(-1) = I - 0.5*[phi]_x + c*[phi]_x^2
    # Use c = 1/theta^2 - 1/(2*theta*tan(theta/2)) for numerical stability.
    c = 1/theta**2 - 1/(2*theta*np.tan(theta/2))
    J_inv = np.eye(3) - 0.5*phi_hat + c*phi_hat @ phi_hat

    # TODO 6: rho = J_inv @ p; return np.r_[phi, rho]
    rho = J_inv @ p
    return np.r_[phi, rho]


# ===================== Exercise 2: implement this =====================

def ik_body(B_list, M, T_d, q0, eomg=1e-4, ev=1e-4, max_iter=100):
    """Numerical Body IK using a pseudoinverse Newton update.

    Returns (q, success, history), with history entries (rotation_error, translation_error).
    """
    q = np.asarray(q0, dtype=float).copy()
    history = []

    for _ in range(max_iter):
        # TODO 1: compute current end-effector pose using fk_body
        T = fk_body(M, B_list, q) 

        # TODO 2: body-frame relative pose T_err = T^(-1) @ T_d
        T_err = trans_inv(T) @ T_d 

        # TODO 3: body pose-error twist using se3_log
        V_err = se3_log(T_err) 

        err_omega = np.linalg.norm(V_err[:3])
        err_v = np.linalg.norm(V_err[3:])
        history.append((err_omega, err_v))

        # TODO 4: return success when BOTH errors are within tolerance
        if err_omega < eomg and err_v < ev:
            return q, True, history

        # TODO 5: compute body Jacobian
        Jb = jacobian_body(B_list, q) 

        # TODO 6: dq = pseudoinverse(Jb) @ V_err
        dq = np.linalg.pinv(Jb) @ V_err

        # TODO 7: update q
        q = q + dq 

    return q, False, history


# ===================== Tests (do not modify) =====================

def test_se3_log():
    print('\n[Exercise 1] Testing se3_log...')
    np.testing.assert_allclose(se3_log(np.eye(4)), np.zeros(6), atol=1e-10)

    T = np.eye(4)
    T[:3, 3] = [2., -1., 3.]
    np.testing.assert_allclose(se3_log(T), [0, 0, 0, 2, -1, 3], atol=1e-10)

    T = np.eye(4)
    T[:3, :3] = so3_exp([0, 0, 1], np.pi / 3)
    np.testing.assert_allclose(se3_log(T), [0, 0, np.pi / 3, 0, 0, 0], atol=1e-9)

    T = np.array([[0, -1, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=float)
    np.testing.assert_allclose(
        se3_log(T), [0, 0, np.pi / 2, np.pi / 4, -np.pi / 4, 0], atol=1e-9
    )

    rng = np.random.default_rng(42)
    for _ in range(100):
        axis = rng.normal(size=3)
        axis /= np.linalg.norm(axis)
        theta = rng.uniform(0.1, 2.5)  # positive angle avoids principal-log sign ambiguity
        rho = rng.uniform(-1., 1., size=3)
        V = np.r_[axis * theta, rho]
        T = se3_exp(V / theta, theta)
        np.testing.assert_allclose(se3_log(T), V, atol=1e-8)
    print('PASS: identity, translation, rotation, coupled pose, 100 random poses')


def example_robot():
    B_list = np.array([[0, 0], [0, 0], [1, 1], [0, 0], [3, 1], [0, 0]], dtype=float)
    M = np.eye(4)
    M[0, 3] = 3.
    return B_list, M


def test_ik_body():
    print('\n[Exercise 2] Testing ik_body...')
    B_list, M = example_robot()
    q_true = np.array([0.5, 0.7])
    T_d = fk_body(M, B_list, q_true)
    q_sol, success, history = ik_body(B_list, M, T_d, np.array([0.4, 0.6]))
    print('Success:', success, '| iterations:', len(history))
    print('True q:', q_true, '| solved q:', q_sol)
    assert success, 'IK did not converge'
    np.testing.assert_allclose(fk_body(M, B_list, q_sol), T_d, atol=1e-4)
    print('PASS: IK recovered the target pose')


def explore_initial_guesses():
    print('\n[Optional] Initial-guess sensitivity:')
    B_list, M = example_robot()
    T_d = fk_body(M, B_list, np.array([0.5, 0.7]))
    for q0 in ([0.4, 0.6], [0., 0.], [-0.5, 1.], [1., -1.], [2., 2.]):
        q_sol, success, history = ik_body(B_list, M, T_d, np.array(q0))
        print(f'q0={q0}, success={success}, steps={len(history)}, final_error={history[-1]}')


if __name__ == '__main__':
    # Work on Exercise 1 first. After it passes, complete Exercise 2.
    test_se3_log()
    test_ik_body()
    explore_initial_guesses()
