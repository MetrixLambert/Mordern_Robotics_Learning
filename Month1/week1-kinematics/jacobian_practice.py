"""ROBO_Learn — Day 6: Space and Body Jacobians.

Run: python day6_jacobian_practice.py

Convention:
  twist = [omega_x, omega_y, omega_z, v_x, v_y, v_z]
  screw arrays have shape (6, n), each COLUMN is a screw axis.

Only implement jacobian_space and jacobian_body. Other functions are provided.
"""
import numpy as np

EPS = 1e-10


def skew(w):
    x, y, z = np.asarray(w, dtype=float).reshape(3)
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=float)


def so3_exp(w, theta):
    w = np.asarray(w, dtype=float).reshape(3)
    a = np.linalg.norm(w)
    if a < EPS:
        return np.eye(3)
    u = w / a
    phi = a * theta
    W = skew(u)
    return np.eye(3) + np.sin(phi) * W + (1 - np.cos(phi)) * (W @ W)


def se3_exp(V, theta):
    """Return exp([V] * theta); works for revolute and prismatic twists."""
    V = np.asarray(V, dtype=float).reshape(6)
    w, v = V[:3], V[3:]
    a = np.linalg.norm(w)
    T = np.eye(4)
    if a < EPS:
        T[:3, 3] = v * theta
        return T
    u = w / a
    phi = a * theta
    W = skew(u)
    R = so3_exp(w, theta)
    G = np.eye(3) * phi + (1 - np.cos(phi)) * W + (phi - np.sin(phi)) * (W @ W)
    T[:3, :3] = R
    T[:3, 3] = G @ (v / a)
    return T


def adjoint(T):
    R, p = T[:3, :3], T[:3, 3]
    Ad = np.zeros((6, 6))
    Ad[:3, :3] = R
    Ad[3:, :3] = skew(p) @ R
    Ad[3:, 3:] = R
    return Ad


def trans_inv(T):
    T_inv = np.eye(4)
    R, p = T[:3, :3], T[:3, 3]
    T_inv[:3, :3] = R.T
    T_inv[:3, 3] = -R.T @ p
    return T_inv


def fk_space(M, S_list, q):
    T = np.eye(4)
    for i in range(len(q)):
        T = T @ se3_exp(S_list[:, i], q[i])
    return T @ M


def fk_body(M, B_list, q):
    T = M.copy()
    for i in range(len(q)):
        T = T @ se3_exp(B_list[:, i], q[i])
    return T


def jacobian_space(S_list, q):
    """Compute J_s(q), shape (6, n). S_list columns are HOME space screws."""
    S_list = np.asarray(S_list, dtype=float)
    q = np.asarray(q, dtype=float).reshape(-1)
    assert S_list.shape == (6, len(q))
    # TODO: Implement using the product of exponentials and Adjoint.
    # Hint: for column i, use the transform from joints 0,...,i-1 ONLY.
    J = np.zeros((6, len(q))) 
    for i in range(len(q)): 
        T = np.eye(4) 
        for j in range(i + 1): 
            T = T @ se3_exp(S_list[:, j], q[j]) 
        J[:, i] = adjoint(T) @ S_list[:, i]
    return J 



def jacobian_body(B_list, q):
    """Compute J_b(q), shape (6, n). B_list columns are HOME body screws."""
    B_list = np.asarray(B_list, dtype=float)
    q = np.asarray(q, dtype=float).reshape(-1)
    assert B_list.shape == (6, len(q))
    # TODO: Implement using reverse traversal of downstream joints.
    # Hint: column n-1 is B_list[:, n-1].
    J = np.zeros((6, len(q))) 
    for i in range(len(q)): 
        T = np.eye(4) 
        for j in range(i, len(q)): 
            T = T @ se3_exp(B_list[:, j], q[j]) 
        J[:, i] = adjoint(trans_inv(T)) @ B_list[:, i]
    return J 


def make_2r():
    M = np.eye(4)
    M[:3, 3] = [3, 0, 0]
    S = np.array([
        [0, 0], [0, 0], [1, 1],
        [0, 0], [0, -2], [0, 0],
    ], dtype=float)
    B = adjoint(trans_inv(M)) @ S
    return M, S, B


def make_rp():
    M = np.eye(4)
    M[:3, 3] = [2, 0, 0]
    S = np.array([
        [0, 0], [0, 0], [1, 0],
        [0, 1], [0, 0], [0, 0],
    ], dtype=float)
    B = adjoint(trans_inv(M)) @ S
    return M, S, B


def assert_close(actual, expected, label, atol=1e-8):
    np.testing.assert_allclose(actual, expected, atol=atol, rtol=0, err_msg=label)
    print(f"  [PASS] {label}")


def test_2r_home():
    M, S, B = make_2r()
    q = np.zeros(2)
    assert_close(jacobian_space(S, q), S, "2R: J_s(0) = S")
    assert_close(jacobian_body(B, q), B, "2R: J_b(0) = B")


def test_2r_specific():
    M, S, B = make_2r()
    q = np.array([np.pi / 2, 0.0])
    Js_expected = np.array([
        [0, 0], [0, 0], [1, 1],
        [0, 2], [0, 0], [0, 0],
    ], dtype=float)
    assert_close(jacobian_space(S, q), Js_expected, "2R: J_s(pi/2, 0)")
    assert_close(jacobian_body(B, q), B, "2R: J_b(pi/2, 0) = B")
    # q1 should not affect J_b in this particular 2R robot.
    assert_close(jacobian_body(B, [0.31, 0.7]),
                 jacobian_body(B, [-1.22, 0.7]), "2R: J_b independent of q1")


def test_rp():
    M, S, B = make_rp()
    q = np.array([np.pi / 2, 1.0])
    Js_expected = np.array([
        [0, 0], [0, 0], [1, 0],
        [0, 0], [0, 1], [0, 0],
    ], dtype=float)
    assert_close(jacobian_space(S, q), Js_expected, "RP: prismatic axis rotates with q1")
    T = fk_space(M, S, q)
    assert_close(jacobian_body(B, q), adjoint(trans_inv(T)) @ Js_expected,
                 "RP: body-space conversion")


def test_random_consistency():
    rng = np.random.default_rng(20261013)
    for n in (1, 2, 4, 6):
        for _ in range(10):
            S = np.empty((6, n))
            for i in range(n):
                if i % 3 == 2:  # prismatic
                    v = rng.normal(size=3)
                    S[:, i] = np.r_[np.zeros(3), v / np.linalg.norm(v)]
                else:  # revolute
                    w = rng.normal(size=3)
                    w /= np.linalg.norm(w)
                    r = rng.normal(size=3)
                    S[:, i] = np.r_[w, -np.cross(w, r)]
            M = np.eye(4)
            axis = rng.normal(size=3)
            axis /= np.linalg.norm(axis)
            M[:3, :3] = so3_exp(axis, rng.uniform(-2, 2))
            M[:3, 3] = rng.normal(size=3)
            B = adjoint(trans_inv(M)) @ S
            q = rng.uniform(-1.5, 1.5, size=n)
            T = fk_space(M, S, q)
            assert_close_silent(fk_body(M, B, q), T)
            Js = jacobian_space(S, q)
            Jb = jacobian_body(B, q)
            assert_close_silent(Jb, adjoint(trans_inv(T)) @ Js)
            assert np.linalg.matrix_rank(Js, tol=1e-8) == np.linalg.matrix_rank(Jb, tol=1e-8)
    print("  [PASS] 40 random robots: FK, J_b = Ad(T^-1) J_s, ranks")


def assert_close_silent(a, b, atol=1e-8):
    np.testing.assert_allclose(a, b, atol=atol, rtol=0)


def test_finite_difference():
    """Compare J_s qdot with numerical derivative of T(q(t))."""
    M, S, B = make_2r()
    q = np.array([0.6, -0.8])
    qdot = np.array([0.4, -0.7])
    h = 1e-6
    T = fk_space(M, S, q)
    T_plus = fk_space(M, S, q + h * qdot)
    T_minus = fk_space(M, S, q - h * qdot)
    Tdot = (T_plus - T_minus) / (2 * h)
    Vs_hat = Tdot @ trans_inv(T)
    Vb_hat = trans_inv(T) @ Tdot
    Vs_numeric = np.r_[np.array([Vs_hat[2, 1], Vs_hat[0, 2], Vs_hat[1, 0]]), Vs_hat[:3, 3]]
    Vb_numeric = np.r_[np.array([Vb_hat[2, 1], Vb_hat[0, 2], Vb_hat[1, 0]]), Vb_hat[:3, 3]]
    assert_close(jacobian_space(S, q) @ qdot, Vs_numeric, "Finite diff: space twist", atol=1e-7)
    assert_close(jacobian_body(B, q) @ qdot, Vb_numeric, "Finite diff: body twist", atol=1e-7)
    # Position Jacobian derived from space twist: p_dot = omega x p + v
    Js = jacobian_space(S, q)
    Jp = -skew(T[:3, 3]) @ Js[:3, :] + Js[3:, :]
    p_dot_numeric = Tdot[:3, 3]
    assert_close(Jp @ qdot, p_dot_numeric, "Finite diff: position Jacobian", atol=1e-7)


def main():
    tests = [test_2r_home, test_2r_specific, test_rp,
             test_random_consistency, test_finite_difference]
    print("Day 6 Jacobian practice — run all tests\n")
    for test in tests:
        print(test.__name__)
        try:
            test()
        except Exception as exc:
            print(f"  [FAIL] {type(exc).__name__}: {exc}")
            raise
    print("\nALL TESTS PASSED!")


if __name__ == "__main__":
    main()
