import numpy as np


# ============================================================
# Basic utilities (Day 2 / Day 3)
# ============================================================

def skew(w):
    """R^3 -> so(3)"""
    x, y, z = w
    return np.array([
        [0, -z, y],
        [z, 0, -x],
        [-y, x, 0]
    ], dtype=float)


def vee(W):
    """so(3) -> R^3"""
    return np.array([
        W[2, 1],
        W[0, 2],
        W[1, 0]
    ], dtype=float)


def so3_exp(phi):
    """
    SO(3) exponential map.

    phi: rotation vector, shape (3,)
    """
    theta = np.linalg.norm(phi)

    if theta < 1e-12:
        return np.eye(3)

    K = skew(phi / theta)

    return (
        np.eye(3)
        + np.sin(theta) * K
        + (1 - np.cos(theta)) * (K @ K)
    )


# ============================================================
# Exercise 1: Twist hat / vee
# ============================================================

def se3_hat(V):
    """
    Convert a 6D twist to a 4x4 se(3) matrix.

    V = [omega_x, omega_y, omega_z, v_x, v_y, v_z]
    """

    w = V[:3] 
    v = V[3:] 
    W = skew(w) 
    M = np.zeros((4, 4)) 
    M[:3, :3] = W 
    M[:3, 3] = v 
    return M 


def se3_vee(M):
    """
    Convert a 4x4 se(3) matrix to a 6D twist.
    """
    w = vee(M[:3, :3]) 
    v = M[:3, 3] 
    V = np.concatenate([w, v]) 
    return V 

# ============================================================
# Exercise 2: Adjoint transformation
# ============================================================

def adjoint(T):
    """
    Compute the 6x6 adjoint matrix of T in SE(3).

    T = [[R, p],
         [0, 1]]

    Convention: V = [omega, v]
    """
    R = T[:3, :3] 
    p = T[:3, 3]  
    Ad = np.zeros((6, 6)) 
    Ad[:3, :3] = R 
    Ad[3:, 3:] = R 
    Ad[3:, :3] = skew(p) @ R 
    return Ad 


# ============================================================
# Exercise 3: SE(3) exponential map
# ============================================================

def se3_exp(V, t):
    """
    Compute exp([V] * t).

    V: constant twist, shape (6,)
    t: time

    Return:
        T: 4x4 homogeneous transformation matrix
    """
    omega = V[:3] 
    v = V[3:] 
    alpha = np.linalg.norm(omega) 
    if alpha < 1e-8: 
        R = np.eye(3) 
        p = t * v 
        
    else: 
        K = omega / alpha 
        K = skew(K) 
        R = so3_exp(omega * t)
        p = t * np.eye(3) + (1 - np.cos(alpha * t)) / alpha * K + (t - np.sin(alpha * t) / alpha) * (K @ K) 
        p = p @ v 

    T = np.zeros((4, 4)) 
    T[:3, :3] = R 
    T[:3, 3] = p 
    T[3, 3] = 1 

    return T 

# ============================================================
# Tests
# ============================================================

def test_exercise_1():
    print("\n===== Exercise 1: se3_hat / se3_vee =====")

    V = np.array([0, 0, 2, 1, 2, 3], dtype=float)

    M = se3_hat(V)

    expected = np.array([
        [0, -2, 0, 1],
        [2,  0, 0, 2],
        [0,  0, 0, 3],
        [0,  0, 0, 0]
    ], dtype=float)

    assert np.allclose(M, expected)
    assert np.allclose(se3_vee(M), V)

    print("Exercise 1 passed!")


def test_exercise_2():
    print("\n===== Exercise 2: adjoint =====")

    T = np.eye(4)
    T[:3, 3] = [1, 0, 0]

    Vb = np.array([0, 0, 2, 0, 2, 0], dtype=float)

    Vs = adjoint(T) @ Vb

    expected = np.array([0, 0, 2, 0, 0, 0], dtype=float)

    assert np.allclose(Vs, expected)

    # Adjoint of identity should be identity
    assert np.allclose(adjoint(np.eye(4)), np.eye(6))

    print("Exercise 2 passed!")


def test_exercise_3():
    print("\n===== Exercise 3: se3_exp =====")

    # Case 1: Pure translation
    V1 = np.array([0, 0, 0, 1, 2, 3], dtype=float)
    T1 = se3_exp(V1, 2.0)

    expected1 = np.eye(4)
    expected1[:3, 3] = [2, 4, 6]

    assert np.allclose(T1, expected1)

    # Case 2: Pure rotation about the z-axis
    V2 = np.array([0, 0, 2, 0, 0, 0], dtype=float)
    T2 = se3_exp(V2, np.pi / 2)

    expected2 = np.eye(4)
    expected2[:3, :3] = np.diag([-1, -1, 1])

    assert np.allclose(T2, expected2, atol=1e-10)

    # Case 3: Rotation about an offset axis
    # omega = [0, 0, 1], q = [1, 0, 0]
    # v = -omega x q = [0, -1, 0]
    V3 = np.array([0, 0, 1, 0, -1, 0], dtype=float)
    T3 = se3_exp(V3, np.pi / 2)

    expected3 = np.array([
        [0, -1, 0,  1],
        [1,  0, 0, -1],
        [0,  0, 1,  0],
        [0,  0, 0,  1]
    ], dtype=float)

    assert np.allclose(T3, expected3, atol=1e-10)

    # Case 4: Screw motion along the z-axis
    V4 = np.array([0, 0, 2, 0, 0, 1], dtype=float)
    T4 = se3_exp(V4, np.pi)

    expected4 = np.eye(4)
    expected4[:3, 3] = [0, 0, np.pi]

    assert np.allclose(T4, expected4, atol=1e-10)

    print("Exercise 3 passed!")

def test_consistency():
    T = np.eye(4)
    T[:3, :3] = so3_exp(np.array([0.2, -0.3, 0.4]))
    T[:3, 3] = [1.0, 2.0, -1.0]

    Vb = np.array([0.3, -0.2, 0.5, 1.0, 0.5, -0.3])
    t = 0.7

    Vs = adjoint(T) @ Vb

    lhs = T @ se3_exp(Vb, t) @ np.linalg.inv(T)
    rhs = se3_exp(Vs, t)

    assert np.allclose(lhs, rhs, atol=1e-10)

    print("SE(3) consistency test passed!")


if __name__ == "__main__":

    # Uncomment exercises after implementing them.

    test_exercise_1()
    test_exercise_2()
    test_exercise_3()

    test_consistency() 