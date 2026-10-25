"""Day 10: Independent numerical checks for robot kinematics.

Run: python day10_jacobian_validation.py
Dependencies: numpy, matplotlib
Conventions: angular-first Twist [omega, v]; T_sb maps body to space.

Complete the three TODO functions; all helpers and tests are provided.
"""
import numpy as np
import matplotlib.pyplot as plt


def skew(v):
    x, y, z = np.asarray(v, dtype=float).reshape(3)
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=float)


def vee_se3(X):
    """Extract [omega, v] from an se(3) matrix."""
    return np.array([X[2, 1], X[0, 2], X[1, 0], *X[:3, 3]])


def hat_se3(V):
    X = np.zeros((4, 4))
    X[:3, :3] = skew(V[:3])
    X[:3, 3] = V[3:]
    return X


def trans_inv(T):
    R, p = T[:3, :3], T[:3, 3]
    X = np.eye(4)
    X[:3, :3] = R.T
    X[:3, 3] = -R.T @ p
    return X


def adjoint(T):
    R, p = T[:3, :3], T[:3, 3]
    A = np.zeros((6, 6))
    A[:3, :3] = R
    A[3:, :3] = skew(p) @ R
    A[3:, 3:] = R
    return A


def se3_exp(Vtheta):
    """Exp of se(3) vector; supports revolute and prismatic twists."""
    w, v = np.asarray(Vtheta[:3], dtype=float), np.asarray(Vtheta[3:], dtype=float)
    theta = np.linalg.norm(w)
    T = np.eye(4)
    if theta < 1e-12:
        T[:3, 3] = v
        return T
    W = skew(w)
    A = np.sin(theta) / theta
    B = (1 - np.cos(theta)) / theta**2
    C = (theta - np.sin(theta)) / theta**3
    T[:3, :3] = np.eye(3) + A * W + B * (W @ W)
    T[:3, 3] = (np.eye(3) + B * W + C * (W @ W)) @ v
    return T


def so3_log_vector(R):
    """Rotation logarithm vector; intended for rotations away from pi."""
    c = np.clip((np.trace(R) - 1) / 2, -1.0, 1.0)
    theta = np.arccos(c)
    if theta < 1e-7:
        return np.array([R[2, 1]-R[1, 2], R[0, 2]-R[2, 0], R[1, 0]-R[0, 1]]) / 2
    return theta / (2 * np.sin(theta)) * np.array([
        R[2, 1]-R[1, 2], R[0, 2]-R[2, 0], R[1, 0]-R[0, 1]
    ])


def se3_log_vector(T):
    """Log_6(T) as angular-first vector, using inverse SO(3) left Jacobian."""
    phi = so3_log_vector(T[:3, :3])
    theta = np.linalg.norm(phi)
    W = skew(phi)
    if theta < 1e-7:
        Jinv = np.eye(3) - 0.5 * W + (W @ W) / 12
    else:
        coeff = 1 / theta**2 - (1 + np.cos(theta)) / (2 * theta * np.sin(theta))
        Jinv = np.eye(3) - 0.5 * W + coeff * (W @ W)
    return np.r_[phi, Jinv @ T[:3, 3]]


def fk_space(Slist, M, q):
    T = np.eye(4)
    for i in range(len(q)):
        T = T @ se3_exp(Slist[:, i] * q[i])
    return T @ M


def jacobian_space(Slist, q):
    J = np.zeros_like(Slist)
    T = np.eye(4)
    for i in range(len(q)):
        J[:, i] = adjoint(T) @ Slist[:, i]
        T = T @ se3_exp(Slist[:, i] * q[i])
    return J


def jacobian_body(Blist, q):
    J = np.zeros_like(Blist)
    T = np.eye(4)
    for i in range(len(q) - 1, -1, -1):
        J[:, i] = adjoint(T) @ Blist[:, i]
        T = T @ se3_exp(-Blist[:, i] * q[i])
    return J


def planar_2r_position(q, L1=2., L2=1.):
    a, b = q
    return np.array([L1*np.cos(a)+L2*np.cos(a+b), L1*np.sin(a)+L2*np.sin(a+b)])


def planar_2r_jacobian(q, L1=2., L2=1.):
    a, b = q
    return np.array([
        [-L1*np.sin(a)-L2*np.sin(a+b), -L2*np.sin(a+b)],
        [ L1*np.cos(a)+L2*np.cos(a+b),  L2*np.cos(a+b)]
    ])


def robot_2r():
    """Space screws (columns), home pose, and corresponding body screws."""
    Slist = np.array([[0., 0.], [0., 0.], [1., 1.],
                      [0., 0.], [0., -2.], [0., 0.]])
    M = np.eye(4)
    M[0, 3] = 3.
    Blist = adjoint(trans_inv(M)) @ Slist
    return Slist, Blist, M


# ===================== EXERCISES =====================
def numerical_position_jacobian(position_fn, q, eps=1e-6):
    """Return m x n Jacobian via central differences of position_fn(q)."""
    q = np.asarray(q, dtype=float)
    # TODO 1: allocate J with shape (len(position_fn(q)), len(q))
    J = np.zeros((len(position_fn(q)), len(q)))

    # TODO 2: for each joint i, evaluate position_fn(q +/- eps*e_i)
    for i in range(len(q)): 
        q_plus = q.copy() 
        q_minus = q.copy() 
        q_plus[i] += eps 
        q_minus[i] -= eps 
        J[:, i] = (position_fn(q_plus) - position_fn(q_minus)) / (2 * eps)

    return J 


def numerical_space_jacobian(fk_fn, q, eps=1e-5):
    """Return 6 x n Space Jacobian using central differences and SE(3) Log.

    For each i: Log_6(T_plus @ inv(T_minus)) / (2*eps).
    """
    q = np.asarray(q, dtype=float)
    # TODO 1: allocate J with shape (6, len(q))
    J = np.zeros((6, len(q)))
    # TODO 2: form T_plus and T_minus for each joint perturbation
    for i in range(len(q)): 
        q_plus = q.copy() 
        q_minus = q.copy() 
        q_plus[i] += eps 
        q_minus[i] -= eps 
        T_plus = fk_fn(q_plus) 
        T_minus = fk_fn(q_minus) 
        J[:, i] = se3_log_vector(T_plus @ trans_inv(T_minus)) / (2 * eps) 
    # TODO 3: relative transform T_plus @ trans_inv(T_minus)
    # TODO 4: apply se3_log_vector and divide by 2*eps

    return J 


def verify_jacobian_consistency(Slist, Blist, M, q, eps=1e-5):
    """Return dict with 'space_fd', 'body_adjoint', 'position_fd' errors."""
    # TODO 1: compute T, analytic Js, analytic Jb
    T = fk_space(Slist, M, q)  
    Js = jacobian_space(Slist, q) 
    Jb = jacobian_body(Blist, q) 

    # TODO 2: compare Js with numerical_space_jacobian(...)
    Js_diff = numerical_space_jacobian(lambda x: fk_space(Slist, M, x), q, eps)
    Js_fd_error = np.linalg.norm(Js_diff - Js, ord='fro') 

    # TODO 3: compare Js with adjoint(T) @ Jb
    Js_adjoint = adjoint(T) @ Jb 
    Js_adjoint_error = np.linalg.norm(Js - Js_adjoint, ord='fro') 
    # Js_error = np.linalg.norm(Js_diff - Js_adjoint, ord='fro') 
    
    # TODO 4: compare planar_2r_jacobian(q) with numerical_position_jacobian(...)
    # Use np.linalg.norm(..., ord='fro') for each error.
    Jp_diff = numerical_position_jacobian(planar_2r_position, q, eps)
    Jp_error = np.linalg.norm(Jp_diff - planar_2r_jacobian(q), ord='fro')

    return {
        'space_fd': Js_fd_error, 
        'body_adjoint': Js_adjoint_error, 
        'position_fd': Jp_error, 
    }

# ===================== TESTS =====================
def test_position_fd():
    q = np.array([0.3, -0.8])
    J_fd = numerical_position_jacobian(planar_2r_position, q)
    J_true = planar_2r_jacobian(q)
    np.testing.assert_allclose(J_fd, J_true, atol=1e-8, rtol=1e-8)
    print('Exercise 1 passed! Position FD error:', np.linalg.norm(J_fd-J_true))


def test_space_fd():
    S, _, M = robot_2r()
    q = np.array([0.3, -0.8])
    J_fd = numerical_space_jacobian(lambda x: fk_space(S, M, x), q)
    J_true = jacobian_space(S, q)
    np.testing.assert_allclose(J_fd, J_true, atol=1e-7, rtol=1e-7)
    print('Exercise 2 passed! Space FD error:', np.linalg.norm(J_fd-J_true))


def test_consistency():
    S, B, M = robot_2r()
    q = np.array([0.3, -0.8])
    errors = verify_jacobian_consistency(S, B, M, q)
    for key in ('space_fd', 'body_adjoint', 'position_fd'):
        assert key in errors, f'Missing result: {key}'
        assert np.isfinite(errors[key]) and errors[key] < 1e-7, (key, errors[key])
        print(f'{key}: {errors[key]:.3e}')
    print('Exercise 3 passed!')


def plot_step_size_study():
    S, _, M = robot_2r()
    q = np.array([0.3, -0.8])
    eps_values = np.logspace(-2, -8, 13)
    position_errors, space_errors = [], []
    for eps in eps_values:
        Jp = numerical_position_jacobian(planar_2r_position, q, eps)
        Js = numerical_space_jacobian(lambda x: fk_space(S, M, x), q, eps)
        position_errors.append(np.linalg.norm(Jp-planar_2r_jacobian(q)))
        space_errors.append(np.linalg.norm(Js-jacobian_space(S, q)))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.loglog(eps_values, position_errors, 'o-', label='Position Jacobian FD')
    ax.loglog(eps_values, space_errors, 's-', label='Space Jacobian SE(3) Log')
    ax.set_xlabel('Finite-difference step epsilon')
    ax.set_ylabel('Frobenius error')
    ax.set_title('Day 10: Finite-difference step-size study')
    ax.grid(True, which='both', alpha=0.3)
    ax.legend()
    fig.tight_layout()
    fig.savefig('day10_fd_errors.png', dpi=180)
    plt.close(fig)
    print('Saved: day10_fd_errors.png')


if __name__ == '__main__':
    np.set_printoptions(precision=6, suppress=True)
    test_position_fd()
    test_space_fd()
    test_consistency()
    plot_step_size_study()
    print('\nAll Day 10 tests passed!')
