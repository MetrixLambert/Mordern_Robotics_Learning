"""
Day 8 - Robot Statics: Wrench, Jacobian Transpose, and Duality

Reference:
    Modern Robotics, Chapter 5.2

Exercises:
    1. Wrench coordinate transformation
    2. Joint generalized forces from an external Wrench
    3. Space/Body static duality verification

Dependencies:
    numpy

Conventions:
    Twist  V = [omega, v]
    Wrench F = [moment, force]
    T_sb maps coordinates from Body frame to Space frame.

All three exercise functions contain TODOs.
"""

import numpy as np


# ============================================================
# Part 1. Auxiliary functions (already implemented)
# ============================================================

def skew(v):
    """
    Convert a 3D vector into a skew-symmetric matrix.

    Parameters
    ----------
    v : ndarray, shape (3,)

    Returns
    -------
    v_hat : ndarray, shape (3, 3)
    """
    v = np.asarray(v, dtype=float).reshape(3)

    return np.array([
        [0.0, -v[2], v[1]],
        [v[2], 0.0, -v[0]],
        [-v[1], v[0], 0.0]
    ])


def so3_exp(omega, theta):
    """
    Exponential map for SO(3).

    R = exp([omega] * theta)

    Assumes omega is a unit axis.
    """
    omega = np.asarray(omega, dtype=float).reshape(3)

    omega_hat = skew(omega)

    return (
        np.eye(3)
        + np.sin(theta) * omega_hat
        + (1.0 - np.cos(theta)) * (omega_hat @ omega_hat)
    )


def adjoint(T):
    """
    Adjoint representation of SE(3).

    Ad_T = [[R,             0],
            [[p] R,         R]]

    Parameters
    ----------
    T : ndarray, shape (4, 4)

    Returns
    -------
    Ad_T : ndarray, shape (6, 6)
    """
    T = np.asarray(T, dtype=float)

    R = T[:3, :3]
    p = T[:3, 3]

    Ad_T = np.zeros((6, 6))

    Ad_T[:3, :3] = R
    Ad_T[3:, :3] = skew(p) @ R
    Ad_T[3:, 3:] = R

    return Ad_T


def make_transform(R, p):
    """
    Construct a homogeneous transformation matrix.
    """
    T = np.eye(4)

    T[:3, :3] = R
    T[:3, 3] = p

    return T


def planar_2r_position_jacobian(q, L1=2.0, L2=1.0):
    """
    Position Jacobian of a planar 2R robot.

    Parameters
    ----------
    q : ndarray, shape (2,)
        Joint angles [q1, q2].

    Returns
    -------
    J_p : ndarray, shape (2, 2)
    """
    q1, q2 = q

    return np.array([
        [
            -L1 * np.sin(q1) - L2 * np.sin(q1 + q2),
            -L2 * np.sin(q1 + q2)
        ],
        [
            L1 * np.cos(q1) + L2 * np.cos(q1 + q2),
            L2 * np.cos(q1 + q2)
        ]
    ])


# ============================================================
# Part 2. Exercise functions (complete the TODOs)
# ============================================================

def wrench_transform(T, F_b):
    """
    Exercise 1: Body Wrench -> Space Wrench.

    Parameters
    ----------
    T : ndarray, shape (4, 4)
        Body frame pose in Space frame.

    F_b : ndarray, shape (6,)
        Body Wrench [moment, force].

    Returns
    -------
    F_s : ndarray, shape (6,)
        Space Wrench.
    """

    # TODO 1:
    AdT = adjoint(T) 

    # TODO 2:
    # Compute the dual transformation:
    #
    #     F_s = Ad_T^(-T) @ F_b
    #
    # Hint:
    # Use np.linalg.solve() rather than np.linalg.inv().
    # Think about which matrix should be transposed.
    F_s = np.linalg.solve(AdT.T, F_b) 

    return F_s


def joint_torque(J, F):
    """
    Exercise 2: Calculate joint generalized forces.

    Parameters
    ----------
    J : ndarray, shape (6, n)
        Space or Body Jacobian.

    F : ndarray, shape (6,)
        Wrench in the matching coordinate frame.

    Returns
    -------
    tau : ndarray, shape (n,)
        Generalized joint forces caused by external loading.
    """

    # TODO:
    # Implement the Jacobian transpose mapping.
    tau = J.T @ F 

    return tau


def verify_static_duality(T, J_b, F_b):
    """
    Exercise 3: Verify static duality.

    Given:
        J_s = Ad_T @ J_b
        F_s = Ad_T^(-T) @ F_b

    Verify:
        J_s.T @ F_s == J_b.T @ F_b

    Parameters
    ----------
    T : ndarray, shape (4, 4)
    J_b : ndarray, shape (6, n)
    F_b : ndarray, shape (6,)

    Returns
    -------
    error : float
        Euclidean norm of the torque difference.
    """

    # TODO 1:
    # Compute J_s from J_b.
    J_s = adjoint(T) @ J_b 

    # TODO 2:
    # Compute F_s using wrench_transform().
    F_s = wrench_transform(T, F_b)

    # TODO 3:
    # Compute tau_s and tau_b using joint_torque().
    tau_s = joint_torque(J_s, F_s)
    tau_b = joint_torque(J_b, F_b)


    # TODO 4:
    # Return the Euclidean norm of their difference.

    return np.linalg.norm(tau_s - tau_b)


# ============================================================
# Part 3. Tests
# ============================================================

def test_exercise_1():
    print("\n" + "=" * 60)
    print("Exercise 1: Wrench Transformation")
    print("=" * 60)

    R = np.eye(3)
    p = np.array([2.0, 0.0, 0.0])

    T = make_transform(R, p)

    F_b = np.array([
        0.0, 0.0, 5.0,     # moment
        0.0, 10.0, 0.0     # force
    ])

    F_s = wrench_transform(T, F_b)

    expected = np.array([
        0.0, 0.0, 25.0,
        0.0, 10.0, 0.0
    ])

    print("Body Wrench:", F_b)
    print("Space Wrench:", F_s)
    print("Expected:", expected)

    np.testing.assert_allclose(F_s, expected, atol=1e-10)

    print("Exercise 1 passed!")


def test_exercise_2():
    print("\n" + "=" * 60)
    print("Exercise 2: Joint Generalized Forces")
    print("=" * 60)

    # Planar 2R robot:
    # L1 = 2, L2 = 1
    # q1 = 0, q2 = pi/2

    q = np.array([0.0, np.pi / 2])

    J_p = planar_2r_position_jacobian(q)

    # Embed the planar position Jacobian in a 6D
    # spatial Jacobian for the end-effector point.
    #
    # This Jacobian maps joint rates to:
    # [angular velocity; linear velocity of the endpoint].
    #
    # Its linear part is the actual endpoint velocity,
    # not the spatial Twist linear component.

    J_endpoint = np.zeros((6, 2))

    J_endpoint[2, :] = [1.0, 1.0]
    J_endpoint[3:5, :] = J_p

    # Force applied at the endpoint.
    # Moment about the endpoint is zero.

    F_endpoint = np.array([
        0.0, 0.0, 0.0,
        10.0, 0.0, 0.0
    ])

    tau = joint_torque(J_endpoint, F_endpoint)

    expected = np.array([-10.0, -10.0])

    print("Position Jacobian:\n", J_p)
    print("Joint generalized forces:", tau)
    print("Expected:", expected)

    np.testing.assert_allclose(tau, expected, atol=1e-10)

    print("Exercise 2 passed!")


def test_exercise_3():
    print("\n" + "=" * 60)
    print("Exercise 3: Static Duality")
    print("=" * 60)

    rng = np.random.default_rng(42)

    # Nontrivial rotation and translation.
    axis = np.array([1.0, 2.0, 3.0])
    axis /= np.linalg.norm(axis)

    R = so3_exp(axis, theta=0.7)

    p = np.array([1.0, 2.0, 3.0])

    T = make_transform(R, p)

    # Random Body Jacobian and Wrench.
    J_b = rng.standard_normal((6, 3))
    F_b = rng.standard_normal(6)

    error = verify_static_duality(T, J_b, F_b)

    print(f"Static duality error: {error:.3e}")

    assert error < 1e-10

    print("Exercise 3 passed!")


def test_pure_couple():
    print("\n" + "=" * 60)
    print("Bonus: Pure Couple Under Translation")
    print("=" * 60)

    F_b = np.array([
        0.0, 0.0, 5.0,
        0.0, 0.0, 0.0
    ])

    T1 = make_transform(
        np.eye(3),
        np.array([2.0, 0.0, 0.0])
    )

    T2 = make_transform(
        np.eye(3),
        np.array([5.0, -3.0, 4.0])
    )

    F_s1 = wrench_transform(T1, F_b)
    F_s2 = wrench_transform(T2, F_b)

    print("Wrench with translation 1:", F_s1)
    print("Wrench with translation 2:", F_s2)

    np.testing.assert_allclose(F_s1, F_s2, atol=1e-10)
    np.testing.assert_allclose(F_s1, F_b, atol=1e-10)

    print("Pure couple test passed!")


# ============================================================
# Part 4. Main
# ============================================================

if __name__ == "__main__":

    np.set_printoptions(precision=6, suppress=True)

    test_exercise_1()
    test_exercise_2()
    test_exercise_3()
    test_pure_couple()

    print("\n" + "=" * 60)
    print("All Day 8 exercises passed!")
    print("=" * 60)