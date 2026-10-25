"""ROBO_Learn Day 12 — RNEA inverse dynamics exercises (NumPy only).

Run: python day12_rnea_exercises.py
Complete TODO 1–3. No third-party robotics library is needed.

Conventions:
  - angular-first twists/wrenches: [omega_x, omega_y, omega_z, v_x, v_y, v_z]
  - T_parent_child maps coordinates from child frame to parent frame
  - link frames are attached at the distal ends of links
  - each A_i is a CONSTANT screw axis expressed in link i's own frame
  - gravity and external wrenches are omitted in this exercise
"""
import numpy as np


def skew(x):
    x, y, z = np.asarray(x, dtype=float).reshape(3)
    return np.array([[0, -z, y], [z, 0, -x], [-y, x, 0]], dtype=float)


def transform(R=None, p=None):
    T = np.eye(4)
    T[:3, :3] = np.eye(3) if R is None else R
    T[:3, 3] = np.zeros(3) if p is None else p
    return T


def rotz(theta):
    c, s = np.cos(theta), np.sin(theta)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)


def inv_transform(T):
    R, p = T[:3, :3], T[:3, 3]
    return transform(R.T, -R.T @ p)


def adjoint(T):
    R, p = T[:3, :3], T[:3, 3]
    result = np.zeros((6, 6))
    result[:3, :3] = R
    result[3:, :3] = skew(p) @ R
    result[3:, 3:] = R
    return result


def ad(V):
    w, v = V[:3], V[3:]
    result = np.zeros((6, 6))
    result[:3, :3] = skew(w)
    result[3:, :3] = skew(v)
    result[3:, 3:] = skew(w)
    return result


def spatial_inertia(m, Ic, c):
    """Spatial inertia in link frame; c = vector from link origin to COM."""
    C = skew(c)
    G = np.zeros((6, 6))
    G[:3, :3] = Ic - m * (C @ C)
    G[:3, 3:] = m * C
    G[3:, :3] = -m * C
    G[3:, 3:] = m * np.eye(3)
    return G


class Planar2R:
    def __init__(self):
        self.L1, self.L2 = 2.0, 1.0
        self.m1, self.m2 = 2.0, 1.0
        self.I1, self.I2 = 0.25, 0.10  # COM moments about z, kg m^2
        # Frames {1}, {2} are at the distal ends of the two links.
        # A_i = [omega, v] in frame {i}; revolute joint axis is z.
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
        """T_{0,1}(q1), T_{1,2}(q2): child -> parent."""
        return [
            transform(rotz(q[0]), rotz(q[0]) @ [self.L1, 0, 0]),
            transform(rotz(q[1]), rotz(q[1]) @ [self.L2, 0, 0]),
        ]


def forward_recursion(robot, q, qd, qdd):
    """Return lists V[i], Vdot[i] for links i=1,2 (0-based Python indices)."""
    T = robot.relative_transforms(q)
    V_prev = np.zeros(6)       # fixed base
    Vdot_prev = np.zeros(6)    # no gravity for now
    Vs, Vdots = [], []
    for i in range(2):
        H = adjoint(inv_transform(T[i]))
        A = robot.A[i]

        V_i = H @ V_prev + A * qd[i] 
        Vdot_i = H @ Vdot_prev + A * qdd[i] + ad(V_i) @ (A * qd[i]) 

        Vs.append(V_i)
        Vdots.append(Vdot_i)
        V_prev, Vdot_prev = V_i, Vdot_i
    return Vs, Vdots


def backward_recursion(robot, q, Vs, Vdots):
    """Return joint torques tau (shape (2,))."""
    T = robot.relative_transforms(q)
    tau = np.zeros(2)
    F_child = np.zeros(6)  # no external end-effector wrench
    for i in range(1, -1, -1):
        # TODO 2: compute the net inertia wrench of link i, then add the
        # wrench transmitted from its child (if any), and project onto A_i.
        # For i=1 (link 2): no child, transmitted wrench is zero.
        # For i=0 (link 1): child wrench is expressed in frame {2};
        # transform it to frame {1} using adjoint(inv_transform(T[1])).T.
        F_i = robot.G[i] @ Vdots[i] - ad(Vs[i]).T @ robot.G[i] @ Vs[i]
        if i < 1: 
            F_i += adjoint(inv_transform(T[i+1])).T @ F_child
        tau[i] = robot.A[i] @ F_i
        F_child = F_i
    return tau


def rnea(robot, q, qd, qdd):
    # TODO 3: connect forward_recursion and backward_recursion.
    Vs, Vdots = forward_recursion(robot, q, qd, qdd)
    tau = backward_recursion(robot, q, Vs, Vdots)
    return tau

def analytic_tau(robot, q, qd, qdd):
    """Independent planar 2R Lagrangian reference, without gravity.

    COM distances from proximal joint: r1=L1/2, r2=L2/2.
    This function is supplied for validation; it does not use RNEA.
    """
    L1, L2 = robot.L1, robot.L2
    m1, m2 = robot.m1, robot.m2
    I1, I2 = robot.I1, robot.I2
    r1, r2 = L1 / 2, L2 / 2
    c, s = np.cos(q[1]), np.sin(q[1])
    M11 = I1 + I2 + m1 * r1**2 + m2 * (L1**2 + r2**2 + 2 * L1 * r2 * c)
    M12 = I2 + m2 * (r2**2 + L1 * r2 * c)
    M22 = I2 + m2 * r2**2
    M = np.array([[M11, M12], [M12, M22]])
    h = m2 * L1 * r2 * s
    velocity_term = np.array([-h * (2 * qd[0] * qd[1] + qd[1]**2),
                              h * qd[0]**2])
    return M @ qdd + velocity_term


def run_checks():
    robot = Planar2R()
    scenarios = [
        ("静止 + 仅关节2加速", [0.4, 0.7], [0, 0], [0, 1]),
        ("两个关节匀速", [0.4, 0.7], [1.2, -0.8], [0, 0]),
        ("混合运动", [-0.8, 1.1], [0.7, 0.3], [0.5, -0.2]),
        ("完全静止", [0.2, -0.5], [0, 0], [0, 0]),
    ]
    for name, q, qd, qdd in scenarios:
        q, qd, qdd = map(lambda x: np.array(x, dtype=float), (q, qd, qdd))
        actual = rnea(robot, q, qd, qdd)
        expected = analytic_tau(robot, q, qd, qdd)
        err = np.max(np.abs(actual - expected))
        print(f"{name}: RNEA={actual}, 解析={expected}, 最大误差={err:.3e}")
        np.testing.assert_allclose(actual, expected, atol=1e-10, rtol=1e-10)

    rng = np.random.default_rng(2026)
    for _ in range(100):
        q = rng.uniform(-np.pi, np.pi, 2)
        qd = rng.uniform(-2, 2, 2)
        qdd = rng.uniform(-3, 3, 2)
        np.testing.assert_allclose(rnea(robot, q, qd, qdd),
                                   analytic_tau(robot, q, qd, qdd),
                                   atol=1e-9, rtol=1e-9)
    print("全部通过：4 组典型案例 + 100 组随机案例。")


if __name__ == "__main__":
    run_checks()
