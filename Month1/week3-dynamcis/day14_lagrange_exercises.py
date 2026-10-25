"""ROBO_Learn | Day 14: Euler-Lagrange dynamics for a vertical planar 2R robot.

NumPy only. Run: python day14_lagrange_exercises.py

Tasks
  TODO 1: derive the mass matrix M(q) from kinetic energy.
  TODO 2: construct Christoffel coefficients and C(q,qd).
  TODO 3: compute gravity potential U(q) and g(q)=grad U.
  TODO 4: assemble inverse dynamics M@qdd + C@qd + g.
  TODO 5 (optional): energy balance and forward dynamics simulation.

Model: vertical XY plane; +Y upward, gravity acceleration [0,-9.81,0].
q1 measured from +X; q2 is link-2 angle relative to link 1.
Joint axes +Z. All inertia values are about link COM.
Day 13's RNEA (no gravity) is supplied as a working reference.
A gravity-aware RNEA reference is also supplied for independent testing.
"""
import numpy as np

# ------------------------------------------------------------
# Working Day 13 spatial algebra and gravity-free RNEA
# ------------------------------------------------------------
def skew(x):
    x, y, z = np.asarray(x, dtype=float).reshape(3)
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=float)


def rotz(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)


def transform(R=None, p=None):
    T = np.eye(4)
    T[:3, :3] = np.eye(3) if R is None else R
    T[:3, 3] = np.zeros(3) if p is None else p
    return T


def inv_transform(T):
    R, p = T[:3, :3], T[:3, 3]
    return transform(R.T, -R.T @ p)


def adjoint(T):
    R, p = T[:3, :3], T[:3, 3]
    A = np.zeros((6, 6))
    A[:3, :3] = R
    A[3:, :3] = skew(p) @ R
    A[3:, 3:] = R
    return A


def ad(V):
    w, v = V[:3], V[3:]
    A = np.zeros((6, 6))
    A[:3, :3] = skew(w)
    A[3:, :3] = skew(v)
    A[3:, 3:] = skew(w)
    return A


def spatial_inertia(m, Ic, c):
    C = skew(c)
    G = np.zeros((6, 6))
    G[:3, :3] = Ic - m * C @ C
    G[:3, 3:] = m * C
    G[3:, :3] = -m * C
    G[3:, 3:] = m * np.eye(3)
    return G


class Planar2R:
    def __init__(self):
        self.L1, self.L2 = 2.0, 1.0
        self.m1, self.m2 = 2.0, 1.0
        self.I1, self.I2 = 0.25, 0.10
        self.n = 2
        self.A = np.array([
            [0, 0, 1, 0, self.L1, 0],
            [0, 0, 1, 0, self.L2, 0],
        ], dtype=float)
        self.G = [
            spatial_inertia(self.m1, np.diag([0.1, 0.1, self.I1]),
                            [-self.L1 / 2, 0, 0]),
            spatial_inertia(self.m2, np.diag([0.1, 0.1, self.I2]),
                            [-self.L2 / 2, 0, 0]),
        ]

    def relative_transforms(self, q):
        return [
            transform(rotz(q[0]), rotz(q[0]) @ [self.L1, 0, 0]),
            transform(rotz(q[1]), rotz(q[1]) @ [self.L2, 0, 0]),
        ]


def forward_recursion(robot, q, qd, qdd):
    T = robot.relative_transforms(q)
    V_prev, Vdot_prev = np.zeros(6), np.zeros(6)
    Vs, Vdots = [], []
    for i in range(robot.n):
        H = adjoint(inv_transform(T[i]))
        A = robot.A[i]
        V_i = H @ V_prev + A * qd[i]
        Vdot_i = H @ Vdot_prev + A * qdd[i] + ad(V_i) @ (A * qd[i])
        Vs.append(V_i)
        Vdots.append(Vdot_i)
        V_prev, Vdot_prev = V_i, Vdot_i
    return Vs, Vdots


def backward_recursion(robot, q, Vs, Vdots):
    T = robot.relative_transforms(q)
    tau = np.zeros(robot.n)
    F_child = np.zeros(6)
    for i in range(robot.n - 1, -1, -1):
        F_i = robot.G[i] @ Vdots[i] - ad(Vs[i]).T @ robot.G[i] @ Vs[i]
        if i < robot.n - 1:
            H_child = adjoint(inv_transform(T[i + 1]))
            F_i += H_child.T @ F_child
        tau[i] = robot.A[i] @ F_i
        F_child = F_i
    return tau


def rnea(robot, q, qd, qdd):
    Vs, Vdots = forward_recursion(robot, q, qd, qdd)
    return backward_recursion(robot, q, Vs, Vdots)




# ------------------------------------------------------------
# Working independent reference: RNEA with gravity
# ------------------------------------------------------------

def rnea_with_gravity(robot, q, qd, qdd, gravity=9.81):
    """RNEA gravity via base acceleration -g_world (Modern Robotics convention).

    World gravity is [0, -gravity, 0]. The artificial base acceleration
    is [0, +gravity, 0] in the world/base frame.
    """
    T = robot.relative_transforms(q)
    V_prev = np.zeros(6)
    Vdot_prev = np.array([0., 0., 0., 0., gravity, 0.])
    Vs, Vdots = [], []
    for i in range(robot.n):
        H = adjoint(inv_transform(T[i]))
        A = robot.A[i]
        V_i = H @ V_prev + A * qd[i]
        Vdot_i = H @ Vdot_prev + A * qdd[i] + ad(V_i) @ (A * qd[i])
        Vs.append(V_i)
        Vdots.append(Vdot_i)
        V_prev, Vdot_prev = V_i, Vdot_i
    return backward_recursion(robot, q, Vs, Vdots)


# ------------------------------------------------------------
# Exercises: write your code inside TODO functions
# ------------------------------------------------------------

def lagrange_mass_matrix(robot, q):
    """TODO 1: analytic M(q) from kinetic energy.

    r1=L1/2, r2=L2/2, q=[q1,q2].
    K=0.5*qd.T @ M(q) @ qd.

    M11 = I1+I2+m1*r1**2 + m2*(L1**2+r2**2+2*L1*r2*cos(q2))
    M12 = I2 + m2*(r2**2+L1*r2*cos(q2))
    M22 = I2 + m2*r2**2

    Return symmetric 2x2 matrix.
    """
    M11 = robot.I1 + robot.I2 + robot.m1 * (robot.L1 / 2)**2 + robot.m2 * (
        robot.L1 ** 2 + (robot.L2 / 2) ** 2 + 2 * robot.L1 * robot.L2 / 2 * np.cos(q[1])
    )
    M12 = robot.I2 + robot.m2 * ((robot.L2 / 2) ** 2 + robot.L1 * robot.L2 / 2 * np.cos(q[1]))
    M22 = robot.I2 + robot.m2 * (robot.L2 / 2)**2
    return np.array([[M11, M12], [M12, M22]])


def christoffel_C(robot, q, qd):
    """TODO 2: construct the canonical C(q,qd), shape (2,2).

    Gamma_ijk = 0.5*(dM_ij/dq_k + dM_ik/dq_j - dM_jk/dq_i)
    C_ij = sum_k Gamma_ijk * qd_k

    For this model let h=m2*L1*r2*sin(q2). Then
      C = [[-h*qd2, -h*(qd1+qd2)],
           [ h*qd1, 0             ]]

    Challenge: derive it from Gamma, rather than simply copying.
    """
    h = robot.m2 * robot.L1 * robot.L2 / 2 * np.sin(q[1])
    C11 = -h * qd[1]
    C12 = -h * (qd[0] + qd[1])
    C21 = h * qd[0]
    C22 = 0
    return np.array([[C11, C12], [C21, C22]])


def potential_energy(robot, q, gravity=9.81):
    """TODO 3a: U(q) in the vertical XY plane.

    y_C1 = r1*sin(q1)
    y_C2 = L1*sin(q1) + r2*sin(q1+q2)
    U = m1*gravity*y_C1 + m2*gravity*y_C2

    Zero potential is chosen when both links are horizontal (q=[0,0]).
    """
    y_C1 = robot.L1 / 2 * np.sin(q[0])
    y_C2 = robot.L1 * np.sin(q[0]) + robot.L2 / 2 * np.sin(q[0] + q[1])
    U = robot.m1 * gravity * y_C1 + robot.m2 * gravity * y_C2
    return U


def gravity_vector(robot, q, gravity=9.81):
    """TODO 3b: return grad_q U(q), shape (2,)."""
    grad_U = np.zeros(2)
    grad_U[0] = robot.m1 * gravity * robot.L1 / 2 * np.cos(q[0]) + robot.m2 * gravity * (robot.L1 * np.cos(q[0]) + robot.L2 / 2 * np.cos(q[0] + q[1]))
    grad_U[1] = robot.m2 * gravity * robot.L2 / 2 * np.cos(q[0] + q[1])
    return grad_U


def lagrange_inverse_dynamics(robot, q, qd, qdd, gravity=9.81):
    """TODO 4: tau = M(q)@qdd + C(q,qd)@qd + g(q)."""
    M = lagrange_mass_matrix(robot, q)
    C = christoffel_C(robot, q, qd)
    g = gravity_vector(robot, q)
    tau = M @ qdd + C @ qd + g
    return tau


def lagrange_forward_dynamics(robot, q, qd, tau, gravity=9.81):
    """TODO 5 (optional): solve M(q)@qdd = tau - C(q,qd)@qd - g(q).

    Use np.linalg.solve rather than computing matrix inverse.
    """
    M = lagrange_mass_matrix(robot, q)
    C = christoffel_C(robot, q, qd)
    g = gravity_vector(robot, q)
    qdd = np.linalg.solve(M, tau - C @ qd - g)
    return qdd


# ------------------------------------------------------------
# Tests: RNEA is an independent dynamics reference
# ------------------------------------------------------------

def test_day14():
    robot = Planar2R()
    rng = np.random.default_rng(1402)
    eps = 1e-6

    print("[1] Mass matrix: 100 random states vs Day 13 RNEA")
    for _ in range(100):
        q = rng.uniform(-np.pi, np.pi, 2)
        M = lagrange_mass_matrix(robot, q)
        M_ref = np.column_stack([
            rnea(robot, q, np.zeros(2), np.eye(2)[i])
            for i in range(2)
        ])
        np.testing.assert_allclose(M, M_ref, atol=1e-10, rtol=1e-10)
        np.testing.assert_allclose(M, M.T, atol=1e-12)
        assert np.linalg.eigvalsh(M).min() > 0
    print("  PASS")

    print("[2] Christoffel C: 100 random states and skew-symmetry")
    for _ in range(100):
        q = rng.uniform(-np.pi, np.pi, 2)
        qd = rng.uniform(-2, 2, 2)
        C = christoffel_C(robot, q, qd)
        c_ref = rnea(robot, q, qd, np.zeros(2))
        np.testing.assert_allclose(C @ qd, c_ref, atol=1e-10, rtol=1e-10)
        dMdt = sum(
            (lagrange_mass_matrix(robot, q + eps * np.eye(2)[j])
             - lagrange_mass_matrix(robot, q - eps * np.eye(2)[j]))
            / (2 * eps) * qd[j]
            for j in range(2)
        )
        S = dMdt - 2 * C
        np.testing.assert_allclose(S + S.T, np.zeros((2, 2)),
                                   atol=1e-8, rtol=0)
    print("  PASS")

    print("[3] Potential U and gravity vector: gradient and gravity RNEA")
    for _ in range(100):
        q = rng.uniform(-np.pi, np.pi, 2)
        grad_fd = np.array([
            (potential_energy(robot, q + eps * np.eye(2)[j])
             - potential_energy(robot, q - eps * np.eye(2)[j]))
            / (2 * eps)
            for j in range(2)
        ])
        g = gravity_vector(robot, q)
        np.testing.assert_allclose(g, grad_fd, atol=1e-7, rtol=1e-7)
        np.testing.assert_allclose(
            g, rnea_with_gravity(robot, q, np.zeros(2), np.zeros(2)),
            atol=1e-10, rtol=1e-10)
    print("  PASS")

    print("[4] Full Lagrange inverse dynamics vs gravity-aware RNEA")
    for _ in range(100):
        q = rng.uniform(-np.pi, np.pi, 2)
        qd = rng.uniform(-2, 2, 2)
        qdd = rng.uniform(-3, 3, 2)
        tau = lagrange_inverse_dynamics(robot, q, qd, qdd)
        tau_ref = rnea_with_gravity(robot, q, qd, qdd)
        np.testing.assert_allclose(tau, tau_ref, atol=1e-9, rtol=1e-9)
    print("  PASS")
    print("\n全部通过：4 组测试，每组 100 次随机采样。")


def test_optional_forward_dynamics():
    robot = Planar2R()
    rng = np.random.default_rng(1414)
    for _ in range(100):
        q = rng.uniform(-np.pi, np.pi, 2)
        qd = rng.uniform(-2, 2, 2)
        qdd_target = rng.uniform(-3, 3, 2)
        tau = rnea_with_gravity(robot, q, qd, qdd_target)
        qdd = lagrange_forward_dynamics(robot, q, qd, tau)
        np.testing.assert_allclose(qdd, qdd_target, atol=1e-9, rtol=1e-9)
    print("PASS: 100 forward/inverse dynamics round trips")


if __name__ == "__main__":
    test_day14()
    # Complete TODO 5, then uncomment:
    test_optional_forward_dynamics()
