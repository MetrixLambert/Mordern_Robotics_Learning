"""Day 12.2: event-driven hybrid simulation of a rimless wheel.

Dependencies: numpy, scipy, matplotlib
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

G = 9.81
L = 1.0
ALPHA = np.deg2rad(15.0)
GAMMA = np.deg2rad(5.0)
C = np.cos(2.0 * ALPHA)
K = 4.0 * G / L * np.sin(GAMMA) * np.sin(ALPHA)
OMEGA_STAR = np.sqrt(C**2 * K / (1.0 - C**2))
N_STEPS = 20
OMEGA0 = 2.0


def dynamics(t, state):
    theta, omega = state
    return [omega, G / L * np.sin(theta + GAMMA)]


def touchdown(t, state):
    return state[0] - ALPHA


touchdown.terminal = True
touchdown.direction = 1


def step_is_possible(omega):
    # For 0 < gamma < alpha, the wheel must cross a potential-energy barrier.
    if 0 < GAMMA < ALPHA:
        threshold = np.sqrt(2.0 * G / L * (1.0 - np.cos(ALPHA - GAMMA)))
        return omega > threshold, threshold
    return omega > 0, 0.0


def main():
    ok, threshold = step_is_possible(OMEGA0)
    if not ok:
        raise ValueError(f"Initial speed {OMEGA0:.4f} must exceed {threshold:.4f} rad/s")

    t0 = 0.0
    state = np.array([-ALPHA, OMEGA0], dtype=float)
    segments = []
    impact_times = []
    post_speeds = [OMEGA0]
    map_errors = []

    for k in range(N_STEPS):
        sol = solve_ivp(
            dynamics, (t0, t0 + 20.0), state,
            events=touchdown, max_step=0.01,
            rtol=1e-10, atol=1e-12,
        )
        if not sol.success or len(sol.t_events[0]) == 0:
            raise RuntimeError(f"Step {k + 1}: touchdown not reached; wheel may have reversed")

        segments.append((sol.t, sol.y))
        t_impact = float(sol.t_events[0][0])
        state_minus = np.array(sol.y_events[0][0])
        omega_plus = C * state_minus[1]
        predicted = C * np.sqrt(state[1] ** 2 + K)
        map_errors.append(abs(omega_plus - predicted))
        impact_times.append(t_impact)
        post_speeds.append(omega_plus)

        state = np.array([-ALPHA, omega_plus])
        t0 = t_impact
        ok, _ = step_is_possible(omega_plus)
        if not ok and k < N_STEPS - 1:
            raise RuntimeError(f"Step {k + 1}: not enough speed for next step")

    print(f"Fixed-point post-impact speed: {OMEGA_STAR:.8f} rad/s")
    print(f"Initial post-impact speed:     {OMEGA0:.8f} rad/s")
    print(f"Final post-impact speed:       {post_speeds[-1]:.8f} rad/s")
    print(f"Stability multiplier:          {C**2:.8f}")
    print(f"Max ODE-vs-map discrepancy:    {max(map_errors):.3e} rad/s")
    print(f"Simulated {N_STEPS} steps; final time: {t0:.4f} s")

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), constrained_layout=True)
    for ts, ys in segments:
        axes[0].plot(ts, np.rad2deg(ys[0]), color='C0', lw=1.5)
        axes[1].plot(ts, ys[1], color='C1', lw=1.5)
    for t in impact_times:
        axes[0].axvline(t, color='0.75', lw=0.5, ls=':')
        axes[1].axvline(t, color='0.75', lw=0.5, ls=':')
    axes[0].set_ylabel('Support-leg angle (deg)')
    axes[0].set_title('Hybrid rimless-wheel motion (each segment ends at impact)')
    axes[1].set_ylabel('Angular speed (rad/s)')
    axes[1].set_xlabel('Time (s)')
    axes[2].plot(range(N_STEPS + 1), post_speeds, 'o-', ms=4)
    axes[2].axhline(OMEGA_STAR, ls='--', color='black', label='Fixed point')
    axes[2].set_xlabel('Impact count k')
    axes[2].set_ylabel('Post-impact speed (rad/s)')
    axes[2].legend()
    for ax in axes:
        ax.grid(True, alpha=0.35)
    output = 'day12_2_hybrid.png'
    fig.savefig(output, dpi=180)
    print(f"Plot saved: {output}")
    plt.show()


if __name__ == '__main__':
    main()
