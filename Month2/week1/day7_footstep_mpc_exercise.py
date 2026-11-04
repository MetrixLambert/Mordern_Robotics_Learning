"""Day 7 — 1D Footstep MPC exercise.

Dependencies: numpy, matplotlib, cvxpy, osqp
Run: python day7_footstep_mpc_exercise.py

TODOs are intentionally left for the learner. The script runs before completion
and prints guidance; after completion it compares three planning strategies.

Model assumptions:
- Single support during [0,T]; fixed current foot centered at 0.
- ZMP remains in the current foot during this interval.
- New foot lands at T. No double-support transfer is modeled yet.
- Terminal DCM offset is xi(T) - p_next = b_des.
- This is a one-step planning demonstration, NOT a closed-loop walking simulation.
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
    N: int = 10
    x0: float = 0.0
    v0: float = 0.4
    foot_current: float = 0.0
    zmp_half_width: float = 0.1
    step_min: float = 0.05
    step_max: float = 0.45
    desired_step_length: float = 0.30
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
    def p_ref(self):
        return self.foot_current + self.desired_step_length

    @property
    def b_des(self):
        # Ideal periodic point-support formula; reference design only.
        return self.desired_step_length / (np.exp(self.omega * self.T) - 1.0)


@dataclass
class Solution:
    name: str
    status: str
    U: np.ndarray | None = None
    p: float | None = None
    cost: float | None = None


def predict_dcm(xi0: float, U: np.ndarray, a: float) -> np.ndarray:
    """Return [xi_0, ..., xi_N] by exact DCM discretization."""
    U = np.asarray(U, dtype=float).reshape(-1)
    xi = np.empty(len(U) + 1)
    xi[0] = xi0
    # TODO 1: Implement xi[j+1] = a*xi[j] + (1-a)*U[j].
    for j in range(len(U)):
        xi[j+1] = a*xi[j] + (1-a)*U[j]
    return xi


def terminal_dcm_expression(xi0: float, U: cp.Expression, a: float, N: int):
    """Construct affine CVXPY expression for xi_N (no Python state loop required)."""
    # TODO 2: Compute weights w_j=a**(N-1-j), then return
    #         a**N*xi0 + (1-a)*(weights @ U).
    weights = np.array([a**(N-1-j) for j in range(N)])
    return a**N*xi0 + (1-a)*(weights @ U)


def solve_footstep_mpc(cfg: Config, fixed_p: float | None = None) -> Solution:
    """Solve a QP; fixed_p=None means jointly optimize U and p.

    With fixed_p set, the footstep is prescribed and only U is optimized.
    """
    U = cp.Variable(cfg.N)
    p = cp.Variable() if fixed_p is None else float(fixed_p)
    u_ref = cfg.foot_current

    # TODO 3: Call terminal_dcm_expression to define xi_N.
    xi_N = terminal_dcm_expression(cfg.xi0, U, cfg.a, cfg.N)

    # TODO 4: Construct objective with:
    #   cfg.r_u * sum_squares(U - u_ref)
    #   cfg.r_du * sum_squares(U[1:] - U[:-1])
    #   cfg.r_p * square(p - cfg.p_ref) [only needed when p is variable]
    obj = cfg.r_u * cp.sum_squares(U - u_ref) + cfg.r_du * cp.sum_squares(U[1:] - U[:-1]) + cfg.r_p * cp.square(p - cfg.p_ref)

    # TODO 5: Construct constraints:
    #   current-foot ZMP limits (foot_current +/- zmp_half_width)
    #   xi_N - p == cfg.b_des
    #   cfg.step_min <= p <= cfg.step_max (when p is variable)
    constraints = [
        U >= cfg.foot_current - cfg.zmp_half_width, 
        U <= cfg.foot_current + cfg.zmp_half_width,
        xi_N - p == cfg.b_des,
    ]
    if fixed_p is None:
        constraints.extend([
            p >= cfg.step_min, 
            p <= cfg.step_max,
        ])

    # TODO 6: Solve with OSQP, check status before reading .value,
    #         return Solution with U, p, cost; for infeasible return status.
    problem = cp.Problem(cp.Minimize(obj), constraints)
    problem.solve()
    return Solution(
        name="Optimized ZMP + " + ("fixed foot" if fixed_p is not None else "optimized foot"),
        status=problem.status,
        U=U.value,
        p=float(p.value) if fixed_p is None else cfg.p_ref,
        cost=problem.value,
    )


def fixed_zmp_baseline(cfg: Config) -> Solution:
    """No optimization: ZMP fixed at current foot center; choose p=p_ref.

    This baseline need not satisfy terminal offset; it is used for comparison.
    """
    U = np.full(cfg.N, cfg.foot_current)
    return Solution("Fixed ZMP + fixed foot", "baseline", U, cfg.p_ref)


def print_result(cfg: Config, sol: Solution):
    print(f"\n[{sol.name}] status={sol.status}")
    if sol.U is None:
        return
    xi = predict_dcm(cfg.xi0, sol.U, cfg.a)
    print(f"  footstep p       = {sol.p:.5f} m")
    print(f"  terminal xi_N    = {xi[-1]:.5f} m")
    print(f"  terminal offset  = {xi[-1] - sol.p:.5f} m")
    print(f"  target offset    = {cfg.b_des:.5f} m")
    print(f"  ZMP range        = [{np.min(sol.U):.5f}, {np.max(sol.U):.5f}] m")
    if sol.cost is not None:
        print(f"  QP objective     = {sol.cost:.6f}")


def plot_comparison(cfg: Config, solutions: list[Solution]):
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
    t_xi = np.arange(cfg.N + 1) * cfg.dt
    t_u = np.arange(cfg.N) * cfg.dt
    for sol in solutions:
        if sol.U is None:
            continue
        xi = predict_dcm(cfg.xi0, sol.U, cfg.a)
        axes[0].plot(t_xi, xi, marker="o", ms=3, label=sol.name)
        axes[0].scatter([cfg.T], [sol.p + cfg.b_des], marker="x", s=55)
        axes[1].step(t_u, sol.U, where="post", label=sol.name)
    axes[0].axhline(cfg.p_ref + cfg.b_des, ls="--", lw=1, label="p_ref + b_des")
    axes[0].set_ylabel("DCM xi (m)")
    axes[0].legend(fontsize=8)
    axes[0].grid(alpha=0.3)
    axes[1].axhline(cfg.foot_current + cfg.zmp_half_width, ls="--", lw=1)
    axes[1].axhline(cfg.foot_current - cfg.zmp_half_width, ls="--", lw=1)
    axes[1].set_xlabel("Time (s)")
    axes[1].set_ylabel("ZMP u (m)")
    axes[1].legend(fontsize=8)
    axes[1].grid(alpha=0.3)
    plt.tight_layout()
    plt.show()


def self_checks(cfg: Config):
    """Basic correctness checks after TODOs are filled."""
    U0 = np.zeros(cfg.N)
    xi = predict_dcm(cfg.xi0, U0, cfg.a)
    np.testing.assert_allclose(xi[-1], cfg.a**cfg.N * cfg.xi0, atol=1e-10)
    U_test = np.linspace(-0.05, 0.08, cfg.N)
    expr = terminal_dcm_expression(cfg.xi0, cp.Constant(U_test), cfg.a, cfg.N)
    np.testing.assert_allclose(float(expr.value), predict_dcm(cfg.xi0, U_test, cfg.a)[-1], atol=1e-10)
    print("DCM prediction checks passed.")


def main():
    cfg = Config()
    print("=== Day 7: Footstep MPC ===")
    print(f"omega={cfg.omega:.5f}, T={cfg.T:.3f}s, xi0={cfg.xi0:.5f}m")
    print(f"p_ref={cfg.p_ref:.5f}m, b_des={cfg.b_des:.5f}m")
    print("Assumption: current single support; no double support modeled.")
    try:
        self_checks(cfg)
        solutions = [
            fixed_zmp_baseline(cfg),
            solve_footstep_mpc(cfg, fixed_p=cfg.p_ref),
            solve_footstep_mpc(cfg, fixed_p=None),
        ]
        solutions[1].name = "Optimized ZMP + fixed foot"
        solutions[2].name = "Optimized ZMP + optimized foot"
        for sol in solutions:
            print_result(cfg, sol)
        plot_comparison(cfg, solutions)
    except NotImplementedError as exc:
        print(f"\nExercise not complete: {exc}")
        print("Fill TODO 1-6 in order, then rerun this file.")


if __name__ == "__main__":
    main()
