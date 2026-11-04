
import numpy as np
import matplotlib.pyplot as plt

# Physical parameters
g = 9.81
h = 1.0
dt = 0.05
omega = np.sqrt(g / h)

# MPC parameters
N = 20
u_min = -0.10
u_max = 0.10

# Initial and target states
x0 = np.array([0.0, 0.2])
x_ref = np.zeros(2)

# TODO 1: Exact ZOH discretization
def get_discrete_lipm(omega, dt):
    # Return Ad (2, 2) and Bd (2, 1)
    Ad = np.array([
        [np.cosh(omega * dt), np.sinh(omega * dt) / omega],
        [omega * np.sinh(omega * dt), np.cosh(omega * dt)]
    ])
    Bd = np.array([
        [1 - np.cosh(omega * dt)], 
        [-omega * np.sinh(omega * dt)] 
    ])
    return Ad, Bd

Ad, Bd = get_discrete_lipm(omega, dt)

print("Ad =\n", Ad)
print("Bd =\n", Bd)
print("Eigenvalues =", np.linalg.eigvals(Ad))


def build_prediction_matrices(Ad, Bd, N):
    nx = Ad.shape[0]
    nu = Bd.shape[1]

    A_bar = np.zeros((N * nx, nx))
    B_bar = np.zeros((N * nx, N * nu))

    # TODO 2:
    # Fill A_bar and B_bar
    for i in range(N):
        A_bar[i * nx: (i+1) * nx, :] = np.linalg.matrix_power(Ad, i+1)
        for j in range(i + 1): 
            B_bar[i * nx: (i+1) * nx, j * nu: (j+1) * nu] = np.linalg.matrix_power(Ad, i - j) @ Bd
    
    return A_bar, B_bar


A_bar, B_bar = build_prediction_matrices(Ad, Bd, N)

assert A_bar.shape == (2 * N, 2)
assert B_bar.shape == (2 * N, N)


rng = np.random.default_rng(42)
U = rng.uniform(u_min, u_max, size=N)

# Method 1: Prediction matrices
X_pred = A_bar @ x0 + B_bar @ U

# Method 2: Step-by-step simulation
x = x0.copy()
states = []

for k in range(N):
    x = Ad @ x + Bd.flatten() * U[k]
    states.append(x.copy())

X_sim = np.concatenate(states)

error = np.max(np.abs(X_pred - X_sim))
print("Maximum prediction error:", error)

assert np.allclose(X_pred, X_sim)
