"""ROBO_Learn — Day 13: joint-space dynamics from RNEA (NumPy only).

Run: python day13_dynamics_exercises.py

Day 12 RNEA is included as a WORKING foundation, so this file runs independently.
Your tasks:
  TODO 1: mass_matrix_from_rnea(q)
  TODO 2: velocity_bias_from_rnea(q, qd)
  TODO 3: forward_dynamics(q, qd, tau)
  TODO 4 (optional): simulate() using semi-implicit Euler

Conventions: angular-first Twist/Wrench; no gravity, friction, or external loads.
The two-link model is identical to Day 12 (L1=2, L2=1, m1=2, m2=1).
"""
import numpy as np


# ============================================================
# Supplied foundation: spatial algebra + Day 12 RNEA
# ============================================================

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


# ============================================================
# Independent analytic reference: planar 2R, no gravity
# ============================================================

def analytic_mass_matrix(robot, q):
    L1, L2 = robot.L1, robot.L2
    m1, m2 = robot.m1, robot.m2
    I1, I2 = robot.I1, robot.I2
    r1, r2 = L1 / 2, L2 / 2
    c = np.cos(q[1])
    M11 = I1 + I2 + m1 * r1**2 + m2 * (L1**2 + r2**2 + 2 * L1 * r2 * c)
    M12 = I2 + m2 * (r2**2 + L1 * r2 * c)
    M22 = I2 + m2 * r2**2
    return np.array([[M11, M12], [M12, M22]])


def analytic_velocity_bias(robot, q, qd):
    h = robot.m2 * robot.L1 * (robot.L2 / 2) * np.sin(q[1])
    return np.array([
        -h * (2 * qd[0] * qd[1] + qd[1] ** 2),
        h * qd[0] ** 2,
    ])


# ============================================================
# Your exercises
# ============================================================

def mass_matrix_from_rnea(robot, q):
    """TODO 1. Return M(q), shape (n,n).

    At zero joint velocity, the k-th column of M is
    RNEA(q, zero_velocity, e_k), where e_k is a unit acceleration.
    Hint: np.eye(robot.n) and M[:, k].
    """
    qd = np.zeros(robot.n) 
    M = np.zeros((robot.n, robot.n)) 

    for k in range(robot.n): 
        qdd = np.zeros(robot.n) 
        qdd[k] = 1 

        tau = rnea(robot, q, qd, qdd)
        M[:, k] = tau 

    return M  


def velocity_bias_from_rnea(robot, q, qd):
    """TODO 2. Return c(q,qd), shape (n,).

    For this NO-GRAVITY model, use one RNEA call with zero acceleration.
    This function returns the velocity-bias VECTOR, not a unique C matrix.
    """
    qdd = np.zeros(robot.n) 
    c = rnea(robot, q, qd, qdd) 
    return c  


def forward_dynamics(robot, q, qd, tau):
    """TODO 3. Return qdd, shape (n,).

    Solve M(q) qdd = tau - c(q,qd).
    Use np.linalg.solve, NOT np.linalg.inv.
    """
    M = mass_matrix_from_rnea(robot, q)
    c = velocity_bias_from_rnea(robot, q, qd)
    qdd = np.linalg.solve(M, tau - c)
    return qdd  


def simulate(robot, q0, qd0, torque_fn, dt=0.002, steps=1000):
    """TODO 4 (optional). Semi-implicit Euler time integration.

    For k = 0,...,steps-1:
      tau = torque_fn(t, q, qd)
      qdd = forward_dynamics(robot, q, qd, tau)
      qd_new = qd + dt * qdd
      q_new = q + dt * qd_new

    Return (times, qs, qds), shapes (steps+1,), (steps+1,n), (steps+1,n).
    Initial state must be included. Avoid modifying q0 or qd0 in place.
    This is a simple educational integrator; numerical energy drift is expected.
    """
    times = np.linspace(0, steps * dt, steps + 1)
    qs = np.zeros((steps + 1, robot.n))
    qds = np.zeros((steps + 1, robot.n))
    qs[0] = q0
    qds[0] = qd0
    for k in range(steps):
        tau = torque_fn(times[k], qs[k], qds[k])
        qdd = forward_dynamics(robot, qs[k], qds[k], tau)
        qds[k + 1] = qds[k] + dt * qdd
        qs[k + 1] = qs[k] + dt * qds[k + 1]
    return times, qs, qds  


# ============================================================
# Checks: independent formulas and physical properties
# ============================================================

def run_checks():
    robot = Planar2R()
    rng = np.random.default_rng(2026)

    print("\n[1] Inertia matrix: symmetry, positive definiteness, analytic check")
    for _ in range(100):
        q = rng.uniform(-np.pi, np.pi, robot.n)
        M = mass_matrix_from_rnea(robot, q)
        assert M.shape == (robot.n, robot.n)
        np.testing.assert_allclose(M, analytic_mass_matrix(robot, q),
                                   rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(M, M.T, rtol=1e-10, atol=1e-10)
        assert np.min(np.linalg.eigvalsh(M)) > 0
    print("  PASS: 100 random configurations")

    print("\n[2] Velocity bias: analytic check + quadratic speed scaling")
    for _ in range(100):
        q = rng.uniform(-np.pi, np.pi, robot.n)
        qd = rng.uniform(-2, 2, robot.n)
        c = velocity_bias_from_rnea(robot, q, qd)
        np.testing.assert_allclose(c, analytic_velocity_bias(robot, q, qd),
                                   rtol=1e-10, atol=1e-10)
        np.testing.assert_allclose(
            velocity_bias_from_rnea(robot, q, 2 * qd), 4 * c,
            rtol=1e-10, atol=1e-10)
    print("  PASS: 100 random states")

    print("\n[3] Forward dynamics: recover acceleration from RNEA torques")
    for _ in range(100):
        q = rng.uniform(-np.pi, np.pi, robot.n)
        qd = rng.uniform(-2, 2, robot.n)
        qdd_target = rng.uniform(-3, 3, robot.n)
        tau = rnea(robot, q, qd, qdd_target)
        qdd_actual = forward_dynamics(robot, q, qd, tau)
        np.testing.assert_allclose(qdd_actual, qdd_target,
                                   rtol=1e-9, atol=1e-9)
    print("  PASS: 100 random inverse/forward dynamics round trips")

    q = np.array([0.4, 0.7])
    qd = np.array([1.2, -0.8])
    M = mass_matrix_from_rnea(robot, q)
    c = velocity_bias_from_rnea(robot, q, qd)
    print("\nExample at q=[0.4, 0.7], qd=[1.2, -0.8]")
    print("M(q) =\n", M)
    print("c(q,qd) =", c)
    print("eigenvalues(M) =", np.linalg.eigvalsh(M))
    print("\n全部通过：M/c/正动力学各 100 组随机测试。")


def run_simulation_check():
    """Optional check. Run after completing TODO 4."""
    robot = Planar2R()
    q0 = np.array([0.4, 0.7])
    qd0 = np.array([0.0, 0.0])
    zero_torque = lambda t, q, qd: np.zeros(robot.n)
    times, qs, qds = simulate(robot, q0, qd0, zero_torque,
                             dt=0.002, steps=1000)
    assert times.shape == (1001,)
    assert qs.shape == (1001, 2)
    assert qds.shape == (1001, 2)
    np.testing.assert_allclose(qs, np.broadcast_to(q0, qs.shape), atol=1e-10)
    np.testing.assert_allclose(qds, np.broadcast_to(qd0, qds.shape), atol=1e-10)
    print("PASS: zero torque + zero initial velocity -> remains at rest")

    # A second experiment with nonzero initial velocity:
    # For a no-gravity, no-friction robot, kinetic energy is conserved
    # in continuous time. Small numerical drift is normal with Euler.
    qd0 = np.array([1.0, -0.5])
    times, qs, qds = simulate(robot, q0, qd0, zero_torque,
                             dt=0.0005, steps=2000)
    energies = np.array([
        0.5 * qds[k] @ analytic_mass_matrix(robot, qs[k]) @ qds[k]
        for k in range(len(times))
    ])
    rel_drift = abs(energies[-1] - energies[0]) / max(1.0, abs(energies[0]))
    print(f"Free-motion kinetic energy: initial={energies[0]:.6f}, "
          f"final={energies[-1]:.6f}, relative drift={rel_drift:.3e}")
    assert rel_drift < 0.02, "Energy drift too large: check integration"
    print("PASS: optional integration tests")


if __name__ == "__main__":
    run_checks()
    # After completing TODO 4, uncomment the next line:
    run_simulation_check()
