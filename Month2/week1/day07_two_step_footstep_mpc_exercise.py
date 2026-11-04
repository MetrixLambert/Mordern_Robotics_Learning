"""Day 7 (Part 2): Two-step 1D Footstep MPC -- TODO exercise.

Run: python day07_two_step_footstep_mpc_exercise.py
Dependencies: numpy, matplotlib, cvxpy, osqp

Model: two consecutive single-support intervals, each T=N*dt.
- During interval 0: ZMP in foot p0's support interval.
- At t=T: ideal instantaneous support switch to foot p1.
- During interval 1: ZMP in foot p1's support interval.
- At t=2T: foot p2 lands; enforce xi(2T)-p2=b_terminal.
- No double support, swing dynamics, impacts, or closed-loop walking.
- b_terminal is a chosen terminal target, not a proof of periodic walking.

Complete TODO 1-5. Keep CVXPY inequalities separate (no a <= x <= b).
"""
from dataclasses import dataclass
import numpy as np
import matplotlib.pyplot as plt
import cvxpy as cp


@dataclass(frozen=True)
class Config:
    g: float = 9.81
    h: float = 1.0
    dt: float = 0.05
    N: int = 10                 # samples per support phase
    x0: float = 0.0
    v0: float = 0.4
    p0: float = 0.0
    half_width: float = 0.10   # foot support interval +/- half_width
    step_min: float = 0.05     # each step relative to preceding foot
    step_max: float = 0.45
    L_ref: float = 0.30
    r_u: float = 1.0
    r_du: float = 1.0
    r_p: float = 10.0

    @property
    def omega(self):
        return np.sqrt(self.g / self.h)

    @property
    def T(self):
        return self.N * self.dt

    @property
    def a(self):
        return np.exp(self.omega * self.dt)

    @property
    def xi0(self):
        return self.x0 + self.v0 / self.omega

    @property
    def b_terminal(self):
        # Point-support periodic reference, used only as an illustrative target.
        return self.L_ref / (np.exp(self.omega * self.T) - 1.0)


@dataclass
class Result:
    name: str
    status: str
    U: np.ndarray | None = None
    p1: float | None = None
    p2: float | None = None
    objective: float | None = None


def predict_dcm(xi0: float, U: np.ndarray, a: float) -> np.ndarray:
    """Given a numeric ZMP sequence, return xi[0], ..., xi[len(U)]."""
    U = np.asarray(U, dtype=float).reshape(-1)
    xi = np.empty(len(U) + 1)
    xi[0] = xi0
    for j, uj in enumerate(U):
        xi[j + 1] = a * xi[j] + (1 - a) * uj
    return xi


def terminal_expr(xi_start, U: cp.Expression, a: float, N: int):
    """Affine CVXPY expression for terminal DCM after N controls."""
    # TODO 1: Implement exact N-step prediction:
    # xi_end = a**N * xi_start + (1-a) * sum(a**(N-1-j) * U[j]).
    # xi_start may itself be a CVXPY affine expression.
    weights = np.array([a**(N-1-j) for j in range(N)])
    return a**N*xi_start + (1-a)*(weights @ U)


def solve_two_step(cfg: Config, fixed_footsteps: bool = False) -> Result:
    """Two-step QP. If fixed_footsteps=True, p1/p2 equal nominal positions."""
    U0 = cp.Variable(cfg.N)
    U1 = cp.Variable(cfg.N)
    p1 = cfg.p0 + cfg.L_ref if fixed_footsteps else cp.Variable()
    p2 = cfg.p0 + 2 * cfg.L_ref if fixed_footsteps else cp.Variable()

    # Provided: phase boundary predictions (TODO 1 must be completed first).
    xi_T = terminal_expr(cfg.xi0, U0, cfg.a, cfg.N)
    xi_2T = terminal_expr(xi_T, U1, cfg.a, cfg.N)

    # TODO 2: Write convex objective using:
    #   r_u * (||U0-p0||^2 + ||U1-p1||^2)
    #   r_du * (||diff(U0)||^2 + ||diff(U1)||^2)
    #   r_p * ((p1-(p0+L_ref))^2 + (p2-(p1+L_ref))^2)
    # The final term is zero when fixed_footsteps=True; it may be retained.
    # Use cp.sum_squares, cp.square. You may optionally penalize the
    # support-switch ZMP jump, but do NOT add it for the initial experiment.
    obj = cfg.r_u * cp.sum_squares(U0 - cfg.p0) + cfg.r_u * cp.sum_squares(U1 - p1) 
    obj = obj + cfg.r_du * cp.sum_squares(U0[1:] - U0[:-1]) + cfg.r_du * cp.sum_squares(U1[1:] - U1[:-1]) 
    if not fixed_footsteps:
        obj = obj + cfg.r_p * cp.sum_squares(p1 - (cfg.p0 + cfg.L_ref)) + cfg.r_p * cp.sum_squares(p2 - (p1 + cfg.L_ref))

    # TODO 3: Write linear constraints:
    # phase 0: p0-half_width <= U0 <= p0+half_width
    # phase 1: p1-half_width <= U1 <= p1+half_width
    # step reach: step_min <= p1-p0 <= step_max
    #             step_min <= p2-p1 <= step_max
    # terminal: xi_2T-p2 == b_terminal
    # IMPORTANT: Split each two-sided inequality into two CVXPY constraints.
    constraints = [
        U0 >= cfg.p0 - cfg.half_width,
        U0 <= cfg.p0 + cfg.half_width,
        U1 >= p1 - cfg.half_width,
        U1 <= p1 + cfg.half_width,
    ]
    if not fixed_footsteps:
        constraints.extend([
            p1 >= cfg.step_min + cfg.p0,
            p1 <= cfg.step_max + cfg.p0,
            p2 >= cfg.step_min + p1,
            p2 <= cfg.step_max + p1,
        ])
    constraints.extend([
        xi_2T - p2 == cfg.b_terminal,
    ])

    # TODO 4: Solve using cp.Problem(cp.Minimize(objective), constraints),
    # solver=cp.OSQP. Return status safely when infeasible; otherwise
    # return Result(name, status, U=np.r_[U0.value, U1.value],
    #               p1=float(...), p2=float(...), objective=float(...)).
    problem = cp.Problem(cp.Minimize(obj), constraints)
    problem.solve()
    return Result(
        name="Two-step Footstep MPC",
        status=problem.status,
        U=np.r_[U0.value, U1.value],
        p1=float(p1.value) if not fixed_footsteps else cfg.p0 + cfg.L_ref,
        p2=float(p2.value) if not fixed_footsteps else cfg.p0 + 2 * cfg.L_ref,
        objective=problem.value,
    )


def validate_solution(cfg: Config, sol: Result):
    """Check dynamics, contact-dependent ZMP bounds and terminal constraint."""
    if sol.U is None:
        print(f"[{sol.name}] status={sol.status} (no trajectory)")
        return
    xi = predict_dcm(cfg.xi0, sol.U, cfg.a)
    U0, U1 = sol.U[:cfg.N], sol.U[cfg.N:]
    print(f"\n[{sol.name}] status={sol.status}")
    print(f"  p0, p1, p2 = {cfg.p0:.5f}, {sol.p1:.5f}, {sol.p2:.5f} m")
    print(f"  steps      = {sol.p1-cfg.p0:.5f}, {sol.p2-sol.p1:.5f} m")
    print(f"  xi(T)      = {xi[cfg.N]:.5f} m")
    print(f"  xi(2T)     = {xi[-1]:.5f} m")
    print(f"  terminal b = {xi[-1]-sol.p2:.5f} m; target={cfg.b_terminal:.5f} m")
    print(f"  phase-0 ZMP relative to p0: [{(U0-cfg.p0).min():.5f}, {(U0-cfg.p0).max():.5f}]")
    print(f"  phase-1 ZMP relative to p1: [{(U1-sol.p1).min():.5f}, {(U1-sol.p1).max():.5f}]")
    print(f"  objective  = {sol.objective:.6f}")
    # TODO 5: Add assertions (tolerance ~1e-4) for:
    #   (a) both ZMP sequences in their corresponding support intervals;
    #   (b) both step lengths in [step_min, step_max];
    #   (c) abs(xi[-1]-sol.p2-b_terminal) < 1e-4.
    # Optionally verify DCM continuity at T and coordinate reset:
    #   (xi[T]-p0) - (xi[T]-p1) == p1-p0.
        # TODO 5: Validate optimization results
    tol = 1e-4

    # (a) ZMP constraints for both support phases
    assert np.all(U0 >= cfg.p0 - cfg.half_width - tol)
    assert np.all(U0 <= cfg.p0 + cfg.half_width + tol)

    assert np.all(U1 >= sol.p1 - cfg.half_width - tol)
    assert np.all(U1 <= sol.p1 + cfg.half_width + tol)

    # (b) Step length constraints
    L0 = sol.p1 - cfg.p0
    L1 = sol.p2 - sol.p1

    assert cfg.step_min - tol <= L0 <= cfg.step_max + tol
    assert cfg.step_min - tol <= L1 <= cfg.step_max + tol

    # (c) Terminal DCM offset constraint
    terminal_error = xi[-1] - sol.p2 - cfg.b_terminal
    assert abs(terminal_error) < tol

    # (d) Optional: DCM coordinate reset at support switch
    b_old = xi[cfg.N] - cfg.p0
    b_new = xi[cfg.N] - sol.p1

    assert abs((b_old - b_new) - L0) < tol

    print("  All constraint checks passed.")


def plot_results(cfg: Config, results: list[Result]):
    good = [s for s in results if s.U is not None]
    if not good:
        return
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    t_xi = np.arange(2 * cfg.N + 1) * cfg.dt
    t_u = np.arange(2 * cfg.N) * cfg.dt
    for sol in good:
        xi = predict_dcm(cfg.xi0, sol.U, cfg.a)
        ax0.plot(t_xi, xi, marker="o", ms=2.5, label=sol.name)
        ax0.scatter([cfg.T, 2 * cfg.T], [sol.p1, sol.p2], marker="x", s=45)
        ax1.step(t_u, sol.U, where="post", label=sol.name)
        ax1.hlines([cfg.p0-cfg.half_width, cfg.p0+cfg.half_width], 0, cfg.T,
                   linestyles="dotted", linewidth=1)
        ax1.hlines([sol.p1-cfg.half_width, sol.p1+cfg.half_width], cfg.T, 2*cfg.T,
                   linestyles="dashed", linewidth=1)
    for ax in (ax0, ax1):
        ax.axvline(cfg.T, color="gray", linestyle=":", linewidth=1)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    ax0.set_ylabel("DCM xi (m)")
    ax1.set_ylabel("ZMP u (m)")
    ax1.set_xlabel("Time (s)")
    ax0.set_title("Two-step Footstep MPC (ideal instantaneous support switch)")
    plt.tight_layout()
    plt.show()


def self_check_prediction(cfg: Config):
    U_test = np.linspace(-0.03, 0.09, 2*cfg.N)
    expr1 = terminal_expr(cfg.xi0, cp.Constant(U_test[:cfg.N]), cfg.a, cfg.N)
    expr2 = terminal_expr(expr1, cp.Constant(U_test[cfg.N:]), cfg.a, cfg.N)
    expected = predict_dcm(cfg.xi0, U_test, cfg.a)
    np.testing.assert_allclose(float(expr1.value), expected[cfg.N], atol=1e-10)
    np.testing.assert_allclose(float(expr2.value), expected[-1], atol=1e-10)
    print("Two-phase DCM prediction checks passed.")


def main():
    cfg = Config()
    print("=== Day 7 Part 2: Two-step Footstep MPC ===")
    print(f"omega={cfg.omega:.5f}, T={cfg.T:.3f}s per phase, xi0={cfg.xi0:.5f}m")
    print(f"L_ref={cfg.L_ref:.3f}m, b_terminal={cfg.b_terminal:.5f}m")
    try:
        self_check_prediction(cfg)
        solutions = [solve_two_step(cfg, fixed_footsteps=True),
                     solve_two_step(cfg, fixed_footsteps=False)]
        for sol in solutions:
            validate_solution(cfg, sol)
        plot_results(cfg, solutions)
    except NotImplementedError as e:
        print(f"Exercise incomplete: {e}")
        print("Complete TODO 1-5, then rerun.")


if __name__ == "__main__":
    main()
