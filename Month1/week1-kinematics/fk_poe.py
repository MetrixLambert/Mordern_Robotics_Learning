
import numpy as np


# ============================================================
# Day 5: Forward Kinematics using Product of Exponentials
#
# Reference:
# Modern Robotics, Chapter 4 - Forward Kinematics
#
# Convention:
# Twist / Screw Axis = [omega_x, omega_y, omega_z,
#                       v_x,     v_y,     v_z]
#
# Slist / Blist: shape (6, n), one screw axis per column
# ============================================================


# ============================================================
# Part 0: Utilities from Day 3 and Day 4
# These functions are already implemented.
# ============================================================

def skew(w):
    """R^3 -> so(3)."""
    x, y, z = np.asarray(w, dtype=float)

    return np.array([
        [0, -z, y],
        [z, 0, -x],
        [-y, x, 0]
    ], dtype=float)


def vee(W):
    """so(3) -> R^3."""
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


def se3_hat(V):
    """Convert a 6D twist into a 4x4 se(3) matrix."""
    V = np.asarray(V, dtype=float)

    M = np.zeros((4, 4))
    M[:3, :3] = skew(V[:3])
    M[:3, 3] = V[3:]

    return M


def se3_vee(M):
    """Convert a 4x4 se(3) matrix into a 6D twist."""
    return np.concatenate([
        vee(M[:3, :3]),
        M[:3, 3]
    ])


def se3_exp(V, t):
    """
    Compute exp([V] * t).

    V: constant twist / screw axis, shape (6,)
    t: time or joint displacement
    """
    V = np.asarray(V, dtype=float)

    omega = V[:3]
    v = V[3:]
    alpha = np.linalg.norm(omega)

    if alpha < 1e-12:
        R = np.eye(3)
        p = t * v
    else:
        K = skew(omega / alpha)

        R = so3_exp(omega * t)

        G = (
            t * np.eye(3)
            + (1 - np.cos(alpha * t)) / alpha * K
            + (t - np.sin(alpha * t) / alpha) * (K @ K)
        )

        p = G @ v

    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = p

    return T


def adjoint(T):
    """
    Compute Ad_T for angular-first twists.

    Ad_T = [[R,       0],
            [[p]R,    R]]
    """
    R = T[:3, :3]
    p = T[:3, 3]

    Ad = np.zeros((6, 6))
    Ad[:3, :3] = R
    Ad[3:, 3:] = R
    Ad[3:, :3] = skew(p) @ R

    return Ad


def trans_inv(T):
    """Inverse of a homogeneous transformation."""
    R = T[:3, :3]
    p = T[:3, 3]

    T_inv = np.eye(4)
    T_inv[:3, :3] = R.T
    T_inv[:3, 3] = -R.T @ p

    return T_inv


# ============================================================
# Exercise 1: Space PoE Forward Kinematics
# ============================================================

def fk_space(M, Slist, q):
    """
    Forward kinematics using Space PoE.

    T(q) = exp([S1]q1) ... exp([Sn]qn) M

    Parameters
    ----------
    M : ndarray, shape (4, 4)
        Home configuration of the end-effector.

    Slist : ndarray, shape (6, n)
        Space screw axes in the HOME configuration.
        Each column represents one joint.

    q : ndarray, shape (n,)
        Joint coordinates:
          revolute -> radians
          prismatic -> meters

    Returns
    -------
    T : ndarray, shape (4, 4)
        End-effector pose relative to the space frame.
    """
    
    # 1. Initialize the transformation matrix.
    T = M.copy() 

    # 2. Loop over all joints.
    # 3. Use se3_exp(Slist[:, i], q[i]).
    for i in range(len(q))[::-1]: 
        T = se3_exp(Slist[:, i], q[i]) @ T 

    return T 

# ============================================================
# Exercise 2: Body PoE Forward Kinematics
# ============================================================

def fk_body(M, Blist, q):
    """
    Forward kinematics using Body PoE.

    T(q) = M exp([B1]q1) ... exp([Bn]qn)

    Parameters
    ----------
    M : ndarray, shape (4, 4)
        Home configuration.

    Blist : ndarray, shape (6, n)
        Body screw axes expressed in the end-effector
        HOME frame.

    q : ndarray, shape (n,)
        Joint coordinates.

    Returns
    -------
    T : ndarray, shape (4, 4)
        End-effector pose relative to the space frame.
    """

    # 1. Initialize the transformation matrix.
    T = M.copy() 
    # 2. Loop over all joints.
    # 3. Multiply the exponentials in the correct order.
    # 4. Remember that M is on the LEFT in Body PoE.
    for i in range(len(q)): 
        T = T @ se3_exp(Blist[:, i], q[i]) 

    return T 


# ============================================================
# Exercise 3: Convert Space Screw Axes to Body Screw Axes
# ============================================================

def space_to_body_screws(M, Slist):
    """
    Convert Space Screw Axes to Body Screw Axes.

    Bi = Ad_(M^-1) Si

    Parameters
    ----------
    M : ndarray, shape (4, 4)
        Home configuration.

    Slist : ndarray, shape (6, n)
        Space screw axes.

    Returns
    -------
    Blist : ndarray, shape (6, n)
        Body screw axes.
    """

    # 1. Compute M inverse.
    M_inv = trans_inv(M) 

    # 2. Compute Ad_(M^-1).
    AdM_inv = adjoint(M_inv)


    # 3. Apply the adjoint to every screw axis.
    Blist = AdM_inv @ Slist

    return Blist 

# ============================================================
# Tests
# ============================================================

def assert_pose_close(T, T_expected, atol=1e-9):
    assert T.shape == (4, 4)
    assert np.allclose(T, T_expected, atol=atol), (
        "\nComputed:\n"
        f"{T}\n"
        "Expected:\n"
        f"{T_expected}"
    )


def make_pose(R, p):
    T = np.eye(4)
    T[:3, :3] = R
    T[:3, 3] = p
    return T


def test_fk_space_2r():
    """
    Planar 2R manipulator.

    L1 = 2, L2 = 1

    S1 = [0, 0, 1, 0,  0, 0]
    S2 = [0, 0, 1, 0, -2, 0]

    M = translation along x by 3.
    """
    L1 = 2.0
    L2 = 1.0

    M = make_pose(np.eye(3), [L1 + L2, 0, 0])

    Slist = np.array([
        [0, 0],
        [0, 0],
        [1, 1],
        [0, 0],
        [0, -L1],
        [0, 0]
    ], dtype=float)

    q = np.array([np.pi / 2, np.pi / 2])

    T = fk_space(M, Slist, q)

    T_expected = make_pose(
        so3_exp([0, 0, np.pi]),
        [-1, 2, 0]
    )

    assert_pose_close(T, T_expected)
    print("[PASS] Exercise 1A: Space PoE - 2R robot")


def test_fk_space_rp():
    """
    RP manipulator.

    Joint 1: revolute about z.
    Joint 2: prismatic along x in home configuration.

    M: end-effector initially at (2, 0, 0).

    theta = pi/2, d = 1
    Expected p = (0, 3, 0).
    """
    M = make_pose(np.eye(3), [2, 0, 0])

    Slist = np.array([
        [0, 0],
        [0, 0],
        [1, 0],
        [0, 1],
        [0, 0],
        [0, 0]
    ], dtype=float)

    q = np.array([np.pi / 2, 1.0])

    T = fk_space(M, Slist, q)

    T_expected = make_pose(
        so3_exp([0, 0, np.pi / 2]),
        [0, 3, 0]
    )

    assert_pose_close(T, T_expected)
    print("[PASS] Exercise 1B: Space PoE - RP robot")


def test_space_to_body():
    """
    For the planar 2R robot:

    M = translation (3, 0, 0).

    B1 = [0, 0, 1, 0, 3, 0]
    B2 = [0, 0, 1, 0, 1, 0]
    """
    M = make_pose(np.eye(3), [3, 0, 0])

    Slist = np.array([
        [0, 0],
        [0, 0],
        [1, 1],
        [0, 0],
        [0, -2],
        [0, 0]
    ], dtype=float)

    Blist = space_to_body_screws(M, Slist)

    B_expected = np.array([
        [0, 0],
        [0, 0],
        [1, 1],
        [0, 0],
        [3, 1],
        [0, 0]
    ], dtype=float)

    assert np.allclose(Blist, B_expected, atol=1e-9)
    print("[PASS] Exercise 3: Space-to-Body screw axes")


def test_fk_body():
    """Compare Space PoE and Body PoE on the same robot."""
    M = make_pose(np.eye(3), [3, 0, 0])

    Slist = np.array([
        [0, 0],
        [0, 0],
        [1, 1],
        [0, 0],
        [0, -2],
        [0, 0]
    ], dtype=float)

    Blist = np.array([
        [0, 0],
        [0, 0],
        [1, 1],
        [0, 0],
        [3, 1],
        [0, 0]
    ], dtype=float)

    q = np.array([0.4, -0.7])

    T_space = fk_space(M, Slist, q)
    T_body = fk_body(M, Blist, q)

    assert_pose_close(T_body, T_space)
    print("[PASS] Exercise 2: Body PoE matches Space PoE")


def test_random_consistency():
    """
    Randomized consistency test.

    Checks:
      FK_space(M, Slist, q)
      ==
      FK_body(M, Ad_(M^-1) Slist, q)

    Uses mixed rotational/translational screw axes.
    """
    rng = np.random.default_rng(42)

    for _ in range(20):
        phi = rng.normal(size=3) * 0.5
        R = so3_exp(phi)
        p = rng.normal(size=3)
        M = make_pose(R, p)

        Slist = np.zeros((6, 3))

        for i in range(3):
            if i == 1:
                # Prismatic joint
                direction = rng.normal(size=3)
                direction /= np.linalg.norm(direction)
                Slist[3:, i] = direction
            else:
                # Revolute joint
                axis = rng.normal(size=3)
                axis /= np.linalg.norm(axis)
                point = rng.normal(size=3)

                Slist[:3, i] = axis
                Slist[3:, i] = -np.cross(axis, point)

        q = rng.uniform(-1, 1, size=3)

        Blist = space_to_body_screws(M, Slist)

        T_space = fk_space(M, Slist, q)
        T_body = fk_body(M, Blist, q)

        assert_pose_close(T_space, T_body, atol=1e-8)

    print("[PASS] Random Space/Body PoE consistency")


# ============================================================
# Main
# ============================================================

if __name__ == "__main__":

    np.set_printoptions(precision=5, suppress=True)

    print("=" * 60)
    print("Day 5: Forward Kinematics using PoE")
    print("=" * 60)

    # Exercise 1
    test_fk_space_2r()
    test_fk_space_rp()

    # # Exercise 2
    test_fk_body()

    # # Exercise 3
    test_space_to_body()

    # Comprehensive test
    test_random_consistency()

    print("=" * 60)
    print("All Day 5 tests passed!")
