"""Day 9: Jacobian manipulability, velocity and force ellipses.

Dependencies: numpy, matplotlib
Run: python day09_manipulability.py

Complete the three TODO functions before running all tests.
"""
import numpy as np
import matplotlib.pyplot as plt

L1, L2 = 2.0, 1.0


def planar_2r_fk(q, L1=L1, L2=L2):
    q1, q2 = np.asarray(q, dtype=float)
    return np.array([L1*np.cos(q1)+L2*np.cos(q1+q2),
                     L1*np.sin(q1)+L2*np.sin(q1+q2)])


def planar_2r_jacobian(q, L1=L1, L2=L2):
    q1, q2 = np.asarray(q, dtype=float)
    return np.array([[-L1*np.sin(q1)-L2*np.sin(q1+q2), -L2*np.sin(q1+q2)],
                     [ L1*np.cos(q1)+L2*np.cos(q1+q2),  L2*np.cos(q1+q2)]])


# ===================== EXERCISES =====================
def ellipse_axes(J):
    """Return U, velocity semi-axes, force semi-axes.

    U columns are principal directions; lengths are descending singular values.
    Force ellipse is defined only for full-rank square J here.
    """
    # TODO 1: SVD of J (np.linalg.svd); get U and singular values s.
    U, s, _ = np.linalg.svd(J)
    # TODO 2: velocity_axes = s.
    velocity_axes = s 

    # TODO 3: if min(s) <= 1e-10, set force_axes to None;
    #         otherwise force_axes = 1 / s.
    if min(s) <= 1e-10:
        force_axes = None
    else:
        force_axes = 1 / s
    
    return U, velocity_axes, force_axes


def yoshikawa_measure(J):
    """w = product of singular values for a full-row-rank task Jacobian.

    Here J is 2x2. At a singular configuration w = 0.
    """
    # TODO: compute w using singular values (not determinant).
    _, s, _ = np.linalg.svd(J)
    w = np.prod(s)
    return w


def jacobian_condition_number(J):
    """Return sigma_max / sigma_min; inf if singular."""
    # TODO: compute singular values; if sigma_min <= 1e-10 return np.inf.
    _, s, _ = np.linalg.svd(J)
    if min(s) <= 1e-10:
        return np.inf
    else:
        return max(s) / min(s)


# ===================== PLOTTING HELPERS =====================
def ellipse_points(U, axes, n=361):
    """Boundary points of a 2D ellipse, shape (2, n)."""
    theta = np.linspace(0, 2*np.pi, n)
    unit_circle = np.vstack((np.cos(theta), np.sin(theta)))
    return U @ (np.asarray(axes)[:, None] * unit_circle)


def plot_ellipse_pair(q, filename):
    J = planar_2r_jacobian(q)
    U, velocity_axes, force_axes = ellipse_axes(J)
    fig, axs = plt.subplots(1, 2, figsize=(10, 4))
    for ax, axes, name in zip(axs, [velocity_axes, force_axes],
                              ['Velocity ellipse', 'Force ellipse']):
        if axes is None:
            ax.text(0.5, 0.5, 'Singular: force set unbounded',
                    transform=ax.transAxes, ha='center', va='center')
        else:
            pts = ellipse_points(U, axes)
            ax.plot(pts[0], pts[1], linewidth=2)
            for i in range(2):
                direction = U[:, i] * axes[i]
                ax.plot([-direction[0], direction[0]],
                        [-direction[1], direction[1]], '--', alpha=0.7)
        ax.axhline(0, linewidth=0.7, alpha=0.5)
        ax.axvline(0, linewidth=0.7, alpha=0.5)
        ax.set_aspect('equal', adjustable='datalim')
        ax.set_title(name)
        ax.set_xlabel('x component')
        ax.set_ylabel('y component')
        ax.grid(True, alpha=0.25)
    fig.suptitle(f'2R robot: q1={q[0]:.3f}, q2={q[1]:.3f} rad')
    fig.tight_layout()
    fig.savefig(filename, dpi=160)
    plt.close(fig)
    print('Saved:', filename)


def plot_metrics(filename):
    q2_values = np.linspace(-np.pi, np.pi, 721)
    w_values = []
    kappa_values = []
    for q2 in q2_values:
        J = planar_2r_jacobian([0.0, q2])
        w_values.append(yoshikawa_measure(J))
        kappa_values.append(jacobian_condition_number(J))
    w_values = np.asarray(w_values)
    kappa_values = np.asarray(kappa_values)
    fig, axs = plt.subplots(2, 1, figsize=(8, 7), sharex=True)
    axs[0].plot(q2_values, w_values, label='Numerical Yoshikawa measure')
    axs[0].plot(q2_values, 2*np.abs(np.sin(q2_values)), '--',
                label='Theory: 2|sin(q2)|')
    axs[0].set_ylabel('w(q)')
    axs[0].legend()
    axs[0].grid(True, alpha=0.25)
    # Hide singularities in the plot (they are infinite, not zero).
    finite_kappa = np.where(np.isfinite(kappa_values), kappa_values, np.nan)
    axs[1].plot(q2_values, finite_kappa)
    axs[1].set_ylim(0, min(30, max(10, np.nanpercentile(finite_kappa, 90))))
    axs[1].set_ylabel('Condition number')
    axs[1].set_xlabel('q2 (rad)')
    axs[1].grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(filename, dpi=160)
    plt.close(fig)
    print('Saved:', filename)
    return q2_values, w_values, kappa_values


# ===================== TESTS =====================
def test_exercise_1():
    print('\nExercise 1: ellipse axes')
    J = np.diag([4.0, 0.5])
    U, velocity_axes, force_axes = ellipse_axes(J)
    np.testing.assert_allclose(velocity_axes, [4.0, 0.5], atol=1e-12)
    np.testing.assert_allclose(force_axes, [0.25, 2.0], atol=1e-12)
    np.testing.assert_allclose(U @ U.T, np.eye(2), atol=1e-12)
    _, _, singular_force = ellipse_axes(np.array([[0., 0.], [3., 1.]]))
    assert singular_force is None
    print('Passed!')


def test_exercise_2():
    print('\nExercise 2: Yoshikawa measure')
    for q2 in [0., np.pi/6, np.pi/2, np.pi, -np.pi/3]:
        J = planar_2r_jacobian([0.4, q2])
        expected = 2 * abs(np.sin(q2))
        np.testing.assert_allclose(yoshikawa_measure(J), expected, atol=1e-10)
    print('Passed!')


def test_exercise_3():
    print('\nExercise 3: condition number')
    np.testing.assert_allclose(jacobian_condition_number(np.diag([4., 0.5])), 8.)
    np.testing.assert_allclose(jacobian_condition_number(np.eye(2)), 1.)
    assert np.isinf(jacobian_condition_number(np.array([[0., 0.], [3., 1.]])))
    print('Passed!')


if __name__ == '__main__':
    np.set_printoptions(precision=6, suppress=True)
    test_exercise_1()
    test_exercise_2()
    test_exercise_3()
    plot_ellipse_pair(np.array([0., np.pi/2]), 'day09_ellipses.png')
    q2, w, kappa = plot_metrics('day09_metrics.png')
    idx_w = np.argmax(w)
    idx_k = np.argmin(kappa)
    print('\nGrid-search results over [-pi, pi]:')
    print(f'Max w: q2={q2[idx_w]:.5f}, w={w[idx_w]:.5f}')
    print(f'Min kappa: q2={q2[idx_k]:.5f}, kappa={kappa[idx_k]:.5f}')
    print('Note: symmetric equivalent optima may exist.')
    print('\nAll Day 9 tests passed!')
