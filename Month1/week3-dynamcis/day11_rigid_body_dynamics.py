
import numpy as np


def vec_to_so3(w):
    """Convert a 3D vector to its skew-symmetric matrix."""
    w = np.asarray(w, dtype=float).reshape(3)

    return np.array([
        [0, -w[2], w[1]],
        [w[2], 0, -w[0]],
        [-w[1], w[0], 0]
    ])


def spatial_inertia(m, I_C):
    """Spatial inertia at the center of mass, angular-first."""
    I_C = np.asarray(I_C, dtype=float).reshape(3, 3)

    # TODO: Construct the 6x6 spatial inertia matrix
    G = np.zeros((6, 6))
    G[:3, :3] = I_C
    G[3:, 3:] = m * np.eye(3)
    return G


def ad_twist(V):
    """Lie algebra adjoint ad_V, angular-first convention."""
    V = np.asarray(V, dtype=float).reshape(6)

    w = V[:3] 
    v = V[3:] 
    Ad = np.zeros((6, 6))
    Ad[:3, :3] = vec_to_so3(w)
    Ad[3:, :3] = vec_to_so3(v)
    Ad[3:, 3:] = vec_to_so3(w)
    return Ad


def newton_euler(G, V, V_dot):
    """Compute the body wrench using Newton-Euler dynamics."""
    G = np.asarray(G, dtype=float).reshape(6, 6)
    V = np.asarray(V, dtype=float).reshape(6)
    V_dot = np.asarray(V_dot, dtype=float).reshape(6)

    F = G @ V_dot - ad_twist(V).T @ G @ V 

    return F

def main():
    m = 3.0
    I_C = np.diag([1.0, 2.0, 3.0])

    V = np.array([1.0, 1.0, 0.0, 1.0, 0.0, 0.0])
    V_dot = np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0])

    G = spatial_inertia(m, I_C)
    F = newton_euler(G, V, V_dot)

    omega = V[:3]
    v = V[3:]
    omega_dot = V_dot[:3]
    v_dot = V_dot[3:]

    # Independent 3D Newton-Euler equations
    moment_3d = I_C @ omega_dot + np.cross(
        omega, I_C @ omega
    )
    force_3d = m * (v_dot + np.cross(omega, v))

    F_expected = np.concatenate([moment_3d, force_3d])

    print("Spatial inertia G:\n", G)
    print("\nad_V:\n", ad_twist(V))
    print("\n6D Newton-Euler result:", F)
    print("3D equations result:   ", F_expected)

    assert np.allclose(F, F_expected, atol=1e-10)
    assert np.allclose(G, G.T)
    assert np.allclose(vec_to_so3(omega) @ omega, 0.0)

    print("\nAll tests passed!")


if __name__ == "__main__":
    main()
