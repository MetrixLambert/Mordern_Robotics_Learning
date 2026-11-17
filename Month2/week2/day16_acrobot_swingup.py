"""ROBO_Learn | Month 2 Day 16: Acrobot energy shaping / swing-up.

Standalone example. Dependencies: numpy, scipy, matplotlib.
Run:
    python month2_day16_acrobot_swingup.py
    python month2_day16_acrobot_swingup.py --mode compare
    python month2_day16_acrobot_swingup.py --mode pfl --no-show

The PFL law is a teaching/experimental controller, not a guaranteed swing-up
or global stabilization law. Torque saturation invalidates exact PFL.
"""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp


@dataclass(frozen=True)
class Params:
    m1: float = 1.0
    m2: float = 1.0
    l1: float = 1.0
    lc1: float = 0.5
    lc2: float = 0.5
    I1: float = 1.0  # about the proximal joint (Day 15 convention)
    I2: float = 1.0  # about the proximal joint (Day 15 convention)
    g: float = 9.81


P = Params()


def dynamics_terms(q: np.ndarray, dq: np.ndarray, p: Params = P):
    q1, q2 = q
    dq1, dq2 = dq
    a = p.m2 * p.l1 * p.lc2
    c2, s2 = np.cos(q2), np.sin(q2)
    M = np.array([
        [p.I1 + p.I2 + p.m2 * p.l1**2 + 2 * a * c2, p.I2 + a * c2],
        [p.I2 + a * c2, p.I2],
    ])
    c_dq = np.array([-a * s2 * (2 * dq1 * dq2 + dq2**2), a * s2 * dq1**2])
    G = np.array([
        (p.m1 * p.lc1 + p.m2 * p.l1) * p.g * np.sin(q1)
        + p.m2 * p.g * p.lc2 * np.sin(q1 + q2),
        p.m2 * p.g * p.lc2 * np.sin(q1 + q2),
    ])
    return M, c_dq, G


def acceleration(q: np.ndarray, dq: np.ndarray, tau: float, p: Params = P):
    M, c_dq, G = dynamics_terms(q, dq, p)
    return np.linalg.solve(M, np.array([0.0, tau]) - c_dq - G)


def total_energy(q: np.ndarray, dq: np.ndarray, p: Params = P) -> float:
    M, _, _ = dynamics_terms(q, dq, p)
    q1, q2 = q
    kinetic = 0.5 * dq @ M @ dq
    potential = -(p.m1 * p.lc1 + p.m2 * p.l1) * p.g * np.cos(q1)
    potential -= p.m2 * p.g * p.lc2 * np.cos(q1 + q2)
    return float(kinetic + potential)


def desired_energy(p: Params = P) -> float:
    return total_energy(np.array([np.pi, 0.0]), np.zeros(2), p)


def pfl_collocated(q: np.ndarray, dq: np.ndarray, v: float, p: Params = P) -> float:
    """Compute actual elbow torque that realizes ddq2=v without saturation."""
    M, c_dq, G = dynamics_terms(q, dq, p)
    h = c_dq + G
    m11, m12, m21, m22 = M[0, 0], M[0, 1], M[1, 0], M[1, 1]
    return float((m22 - m21 * m12 / m11) * v + h[1] - m21 * h[0] / m11)


def control(x: np.ndarray, mode: str, p: Params = P,
            k1: float = 10.0, k2: float = 5.0, k3: float = 1.0,
            k_energy: float = 1.0, tau_max: float | None = 50.0) -> tuple[float, float]:
    q, dq = x[:2], x[2:]
    e_err = total_energy(q, dq, p) - desired_energy(p)
    if mode == 'pfl':
        # Virtual acceleration, not actual motor torque.
        v = -k1 * q[1] - k2 * dq[1] + k3 * dq[0] * e_err
        tau_raw = pfl_collocated(q, dq, v, p)
    elif mode == 'direct':
        # Direct motor-power-based energy feedback: E_dot = dq2*tau.
        tau_raw = -k_energy * dq[1] * e_err
    else:
        raise ValueError(f'Unknown mode: {mode}')
    tau = float(np.clip(tau_raw, -tau_max, tau_max)) if tau_max is not None else float(tau_raw)
    return tau, float(tau_raw)


def simulate(mode: str, duration: float, dt: float, angle0: float,
             tau_max: float | None, p: Params = P):
    x0 = np.array([np.deg2rad(angle0), 0.0, 0.0, 0.0])
    t_eval = np.linspace(0, duration, int(round(duration / dt)) + 1)

    def rhs(_t, x):
        tau, _ = control(x, mode, p=p, tau_max=tau_max)
        return np.r_[x[2:], acceleration(x[:2], x[2:], tau, p)]

    sol = solve_ivp(rhs, (0.0, duration), x0, t_eval=t_eval,
                    rtol=1e-8, atol=1e-10, max_step=0.01)
    if not sol.success:
        raise RuntimeError(sol.message)
    energies = np.array([total_energy(sol.y[:2, i], sol.y[2:, i], p)
                         for i in range(sol.t.size)])
    torques = np.array([control(sol.y[:, i], mode, p=p, tau_max=tau_max)[0]
                        for i in range(sol.t.size)])
    raw_torques = np.array([control(sol.y[:, i], mode, p=p, tau_max=tau_max)[1]
                            for i in range(sol.t.size)])
    # The exact continuous-time identity is E_dot = tau * dq2.
    power = torques * sol.y[3]
    # A centered finite difference gives a numerical consistency check.
    energy_dot_numeric = np.gradient(energies, sol.t)
    return sol, energies, torques, raw_torques, power, energy_dot_numeric


def plot_result(axs, result, label: str, p: Params = P):
    sol, energies, torques, _, power, energy_dot_numeric = result
    axs[0].plot(sol.t, energies, label=label)
    axs[1].plot(sol.t, np.rad2deg(sol.y[0]), label=f'{label}: q1')
    axs[1].plot(sol.t, np.rad2deg(sol.y[1]), '--', label=f'{label}: q2')
    axs[2].plot(sol.t, torques, label=label)
    axs[3].plot(sol.t, energy_dot_numeric - power, label=label)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mode', choices=['pfl', 'direct', 'compare'], default='compare')
    parser.add_argument('--duration', type=float, default=15.0)
    parser.add_argument('--dt', type=float, default=0.005)
    parser.add_argument('--angle0', type=float, default=2.865,
                        help='Initial q1 angle in degrees; 2.865 deg is about 0.05 rad')
    parser.add_argument('--tau-max', type=float, default=50.0,
                        help='Motor torque limit in Nm; use --no-limit to disable')
    parser.add_argument('--no-limit', action='store_true')
    parser.add_argument('--no-show', action='store_true', help='Save plot without opening a window')
    parser.add_argument('--output', type=Path, default=Path('month2_day16_swingup.png'))
    args = parser.parse_args()
    if args.duration <= 0 or args.dt <= 0 or args.tau_max <= 0:
        parser.error('duration, dt, and tau-max must be positive')
    tau_max = None if args.no_limit else args.tau_max
    modes = ['pfl', 'direct'] if args.mode == 'compare' else [args.mode]
    results = {}
    print(f'Desired upright energy: {desired_energy():.4f} J')
    print(f'Initial energy: {total_energy(np.deg2rad([args.angle0, 0.0]), np.zeros(2)):.4f} J')
    for mode in modes:
        res = simulate(mode, args.duration, args.dt, args.angle0, tau_max)
        results[mode] = res
        sol, E, tau, raw, power, edot = res
        print(f'\n[{mode}] final energy={E[-1]:.4f} J, energy error={E[-1]-desired_energy():+.4f} J')
        print(f'  max |tau|={np.max(np.abs(tau)):.3f} Nm; saturation samples={np.count_nonzero(np.abs(raw-tau)>1e-8)}')
        print(f'  power identity: RMS |dE/dt - tau*dq2|={np.sqrt(np.mean((edot-power)**2)):.4g} W')
        # Wrapped distance to upright is only a diagnostic, not a success certificate.
        q1err = (sol.y[0, -1] - np.pi + np.pi) % (2 * np.pi) - np.pi
        q2err = (sol.y[1, -1] + np.pi) % (2 * np.pi) - np.pi
        print(f'  final wrapped q errors: q1={np.rad2deg(q1err):+.2f} deg, q2={np.rad2deg(q2err):+.2f} deg')

    fig, axs = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    for mode, result in results.items():
        plot_result(axs, result, mode)
    axs[0].axhline(desired_energy(), color='gray', ls=':', label='target E')
    axs[0].set_ylabel('Total energy (J)')
    axs[1].axhline(180, color='gray', ls=':', label='upright q1 = 180 deg')
    axs[1].set_ylabel('Angles (deg, unwrapped)')
    axs[2].set_ylabel('Elbow torque (Nm)')
    axs[3].set_ylabel('dE/dt - tau*dq2 (W)')
    axs[3].set_xlabel('Time (s)')
    for ax in axs:
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, ncol=2)
    fig.suptitle('ROBO_Learn — Month 2 Day 16: Acrobot swing-up experiments')
    fig.tight_layout()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, dpi=160)
    print(f'\nFigure saved: {args.output.resolve()}')
    if not args.no_show:
        plt.show()


if __name__ == '__main__':
    main()
