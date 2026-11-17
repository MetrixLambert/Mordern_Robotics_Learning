"""Day 12.1: Rimless wheel Poincare map and local stability.
Run: python run_12_1.py
"""
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

g = 9.81
l = 1.0
alpha = np.deg2rad(15.0)
gamma = np.deg2rad(5.0)
c = np.cos(2.0 * alpha)
K = 4.0 * g / l * np.sin(gamma) * np.sin(alpha)


def poincare_map(omega):
    """One step, provided the next touchdown is physically reachable."""
    return c * np.sqrt(omega**2 + K)


# Starting at theta=-alpha, the wheel must pass the unstable equilibrium
# theta=-gamma before reaching theta=+alpha (when alpha > gamma).
barrier = np.sqrt(2 * g / l * (1 - np.cos(alpha - gamma)))
omega_star = np.sqrt(c**2 * K / (1 - c**2))
multiplier = c**2
print(f"c = cos(2 alpha): {c:.6f}")
print(f"K: {K:.6f} (rad/s)^2")
print(f"Minimum post-impact speed to clear potential barrier: {barrier:.6f} rad/s")
print(f"Fixed point omega*: {omega_star:.6f} rad/s")
print(f"P(omega*) - omega*: {poincare_map(omega_star) - omega_star:.3e}")
print(f"Local stability multiplier P'(omega*): {multiplier:.6f}")

n_steps = 30
initial_speeds = [0.65, 0.95, 1.40, 2.00]
fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
for omega0 in initial_speeds:
    if omega0 <= barrier:
        raise ValueError(f"omega0={omega0} cannot clear the potential barrier")
    history = [omega0]
    for _ in range(n_steps):
        history.append(poincare_map(history[-1]))
    history = np.array(history)
    axes[0].plot(range(n_steps + 1), history, marker='.', label=f'omega0={omega0:.2f}')
    axes[1].semilogy(range(n_steps + 1), np.maximum(abs(history - omega_star), 1e-15), marker='.')
    print(f"omega0={omega0:.2f} -> omega30={history[-1]:.6f}, error={history[-1]-omega_star:+.3e}")

axes[0].axhline(omega_star, linestyle='--', color='black', label='Fixed point')
axes[0].set_ylabel('Post-impact angular speed (rad/s)')
axes[0].set_title('Rimless wheel: Poincare map convergence')
axes[0].legend()
axes[0].grid(True)
axes[1].set_xlabel('Step number k')
axes[1].set_ylabel('Absolute error |omega_k - omega*|')
axes[1].grid(True, which='both')
fig.tight_layout()
output = Path(__file__).with_name('day12_1_poincare.png')
fig.savefig(output, dpi=160)
print(f'Plot saved: {output}')
plt.show()
