import numpy as np 
import cvxpy as cp
from day6_lipm_mpc_sys import *

Q = np.diag([10.0, 1.0])
r = 0.01

Q_bar = np.kron(np.eye(N), Q)
R_bar = r * np.eye(N)


def solve_mpc(x_current, A_bar, B_bar, Q_bar, R_bar,
              u_min, u_max, omega=None, terminal_dcm=False):

    N = B_bar.shape[1]
    U = cp.Variable(N)

    X_pred = A_bar @ x_current + B_bar @ U

    objective = (
        cp.quad_form(X_pred, Q_bar)
        + cp.quad_form(U, R_bar)
    )

    constraints = [
        U >= u_min,
        U <= u_max
    ]

    if terminal_dcm:
        X_terminal = X_pred[-2:]

        xi_terminal = (
            X_terminal[0]
            + X_terminal[1] / omega
        )

        constraints += [
            xi_terminal >= u_min,
            xi_terminal <= u_max
        ]

    problem = cp.Problem(
        cp.Minimize(objective),
        constraints
    )

    problem.solve(solver=cp.OSQP)

    if problem.status not in ("optimal", "optimal_inaccurate"):
        return None, None, problem.status

    return float(U.value[0]), U.value.copy(), problem.status



sim_time = 5.0
sim_steps = int(sim_time / dt)

x = x0.copy()

state_history = [x.copy()]
zmp_history = []

for k in range(sim_steps):

    # TODO 5: Solve MPC
    u, U_opt, status = solve_mpc(
        x, A_bar, B_bar, Q_bar, R_bar,
        u_min, u_max, omega,
        terminal_dcm=True
    )

    if u is None:
        print(f"MPC infeasible at step {k}: {status}")
        break

    x = Ad @ x + Bd.flatten() * u

    state_history.append(x.copy())
    zmp_history.append(u)

state_history = np.array(state_history)
zmp_history = np.array(zmp_history)


# ============================
# Plot MPC simulation results
# ============================

time_state = np.arange(sim_steps + 1) * dt
time_control = np.arange(sim_steps) * dt

# Compute DCM
dcm_history = (
    state_history[:, 0]
    + state_history[:, 1] / omega
)

fig, axes = plt.subplots(
    4, 1, figsize=(10, 11), sharex=True
)

# 1. CoM position
axes[0].plot(
    time_state, state_history[:, 0],
    label="CoM position"
)
axes[0].axhline(0, color="gray", linestyle="--")
axes[0].set_ylabel("Position (m)")
axes[0].legend()
axes[0].grid(True)

# 2. CoM velocity
axes[1].plot(
    time_state, state_history[:, 1],
    label="CoM velocity"
)
axes[1].axhline(0, color="gray", linestyle="--")
axes[1].set_ylabel("Velocity (m/s)")
axes[1].legend()
axes[1].grid(True)

# 3. ZMP control input
axes[2].step(
    time_control, zmp_history,
    where="post", label="ZMP"
)
axes[2].axhline(
    u_max, color="red", linestyle="--",
    label="ZMP limits"
)
axes[2].axhline(
    u_min, color="red", linestyle="--"
)
axes[2].set_ylabel("ZMP (m)")
axes[2].legend()
axes[2].grid(True)

# 4. DCM
axes[3].plot(
    time_state, dcm_history,
    label="DCM"
)
axes[3].axhline(
    u_max, color="red", linestyle="--",
    label="Support limits"
)
axes[3].axhline(
    u_min, color="red", linestyle="--"
)
axes[3].axhline(0, color="gray", linestyle=":")
axes[3].set_ylabel("DCM (m)")
axes[3].set_xlabel("Time (s)")
axes[3].legend()
axes[3].grid(True)

fig.suptitle("LIPM-MPC Balance Control", fontsize=14)
plt.tight_layout()
plt.show()


# ============================
# Effect of MPC Prediction Horizon
# ============================

N_list = [5, 10, 20, 40]

results = {}

for N_test in N_list:
    A_test, B_test = build_prediction_matrices(
        Ad, Bd, N_test
    )

    Q_test = np.kron(np.eye(N_test), Q)
    R_test = r * np.eye(N_test)

    x = x0.copy()
    states = [x.copy()]
    controls = []

    for k in range(sim_steps):
        u, _, status = solve_mpc(
            x, A_test, B_test,
            Q_test, R_test,
            u_min, u_max
        )

        if u is None:
            print(f"MPC infeasible at step {k}: {status}")
            break

        x = Ad @ x + Bd.flatten() * u

        states.append(x.copy())
        controls.append(u)

    states = np.array(states)
    controls = np.array(controls)

    results[N_test] = {
        "states": states,
        "controls": controls,
        "dcm": states[:, 0] + states[:, 1] / omega
    }

# Compare DCM trajectories
plt.figure(figsize=(10, 5))

for N_test in N_list:
    plt.plot(
        time_state,
        results[N_test]["dcm"],
        label=f"N = {N_test}"
    )

plt.axhline(u_max, linestyle="--", color="red")
plt.axhline(u_min, linestyle="--", color="red")
plt.axhline(0, linestyle=":", color="gray")

plt.xlabel("Time (s)")
plt.ylabel("DCM (m)")
plt.title("Effect of MPC Prediction Horizon")
plt.legend()
plt.grid(True)
plt.show()



# ======================================
# Experiment: Different initial velocities
# ======================================

v0_list = [0.1, 0.2, 0.3]
N_test = 20

A_test, B_test = build_prediction_matrices(
    Ad, Bd, N_test
)

Q_test = np.kron(np.eye(N_test), Q)
R_test = r * np.eye(N_test)

velocity_results = {}

for v0 in v0_list:
    x = np.array([0.0, v0])

    states = [x.copy()]
    controls = []

    for k in range(sim_steps):
        u, _, status = solve_mpc(
            x, A_test, B_test,
            Q_test, R_test,
            u_min, u_max
        )

        if u is None:
            print(f"MPC infeasible at step {k}: {status}")
            break

        x = Ad @ x + Bd.flatten() * u

        states.append(x.copy())
        controls.append(u)

    states = np.array(states)
    controls = np.array(controls)

    dcm = states[:, 0] + states[:, 1] / omega

    velocity_results[v0] = {
        "states": states,
        "controls": controls,
        "dcm": dcm
    }

fig, axes = plt.subplots(
    3, 1, figsize=(10, 10), sharex=True
)

for v0 in v0_list:
    result = velocity_results[v0]
    label = f"v0 = {v0:.1f} m/s"

    axes[0].plot(
        time_state,
        result["states"][:, 0],
        label=label
    )

    axes[1].plot(
        time_state,
        result["states"][:, 1],
        label=label
    )

    axes[2].plot(
        time_state,
        result["dcm"],
        label=label
    )

axes[0].set_ylabel("CoM position (m)")
axes[0].set_title("Effect of Initial Velocity")

axes[1].set_ylabel("CoM velocity (m/s)")

axes[2].axhline(
    u_max, color="red", linestyle="--",
    label="Support limits"
)
axes[2].axhline(
    u_min, color="red", linestyle="--"
)
axes[2].set_ylabel("DCM (m)")
axes[2].set_xlabel("Time (s)")

for ax in axes:
    ax.grid(True)
    ax.legend()

plt.tight_layout()
plt.show()


# ============================
# Effect of Terminal DCM constraints
# ============================

for v0 in [0.2, 0.3]:
    x_init = np.array([0.0, v0])

    u, _, status = solve_mpc(
        x_init,
        A_bar, B_bar,
        Q_bar, R_bar,
        u_min, u_max,
        omega,
        terminal_dcm=True
    )

    if u is None:
        print(f"MPC infeasible at step {k}: {status}")
        break

    print(
        f"v0 = {v0:.1f} m/s, "
        f"status = {status}, "
        f"u0 = {u}"
    )


# ========================================
# Compare baseline MPC vs terminal DCM MPC
# ========================================

v0 = 0.3
N_test = 20
sim_time = 5.0
sim_steps = int(sim_time / dt)

A_test, B_test = build_prediction_matrices(
    Ad, Bd, N_test
)

Q_test = np.kron(np.eye(N_test), Q)
R_test = r * np.eye(N_test)

results = {}

for terminal_dcm in [False, True]:
    x = np.array([0.0, v0])

    states = [x.copy()]
    controls = []
    terminal_dcm_predictions = []
    statuses = []

    for k in range(sim_steps):
        u, U_opt, status = solve_mpc(
            x, A_test, B_test,
            Q_test, R_test,
            u_min, u_max,
            omega,
            terminal_dcm=terminal_dcm
        )

        statuses.append(status)

        if u is None:
            print(
                f"terminal_dcm={terminal_dcm}, "
                f"step={k}, status={status}"
            )
            break

        # Predicted terminal state at current MPC step
        X_pred = A_test @ x + B_test @ U_opt
        x_terminal = X_pred[-2:]

        xi_terminal = (
            x_terminal[0]
            + x_terminal[1] / omega
        )
        terminal_dcm_predictions.append(xi_terminal)

        # Apply only the first optimal control
        x = Ad @ x + Bd.flatten() * u

        states.append(x.copy())
        controls.append(u)

    states = np.array(states)
    controls = np.array(controls)

    results[terminal_dcm] = {
        "states": states,
        "controls": controls,
        "dcm": states[:, 0] + states[:, 1] / omega,
        "terminal_dcm_predictions": np.array(
            terminal_dcm_predictions
        ),
        "statuses": statuses
    }

# ========================================
# Plot
# ========================================

fig, axes = plt.subplots(
    4, 1, figsize=(11, 12), sharex=True
)

for terminal_dcm in [False, True]:
    result = results[terminal_dcm]

    label = (
        "MPC + DCM terminal constraint"
        if terminal_dcm else "Baseline MPC"
    )

    ts = np.arange(len(result["states"])) * dt
    tu = np.arange(len(result["controls"])) * dt

    axes[0].plot(
        ts, result["states"][:, 0], label=label
    )
    axes[1].plot(
        ts, result["states"][:, 1], label=label
    )
    axes[2].step(
        tu, result["controls"], where="post",
        label=label
    )
    axes[3].plot(
        ts, result["dcm"], label=label
    )

axes[0].set_ylabel("CoM position (m)")
axes[1].set_ylabel("CoM velocity (m/s)")
axes[2].set_ylabel("ZMP (m)")
axes[3].set_ylabel("DCM (m)")
axes[3].set_xlabel("Time (s)")

for ax in [axes[2], axes[3]]:
    ax.axhline(u_max, color="red", linestyle="--")
    ax.axhline(u_min, color="red", linestyle="--")

for ax in axes:
    ax.grid(True)
    ax.legend()

plt.tight_layout()


plt.figure(figsize=(10, 5))

for terminal_dcm in [False, True]:
    result = results[terminal_dcm]
    xi_pred = result["terminal_dcm_predictions"]

    t = np.arange(len(xi_pred)) * dt

    label = (
        "MPC + DCM terminal constraint"
        if terminal_dcm else "Baseline MPC"
    )

    plt.plot(t, xi_pred, label=label)

plt.axhline(u_max, color="red", linestyle="--")
plt.axhline(u_min, color="red", linestyle="--")
plt.axhline(0, color="gray", linestyle=":")

plt.xlabel("Current MPC time (s)")
plt.ylabel("Predicted terminal DCM (m)")
plt.title("Predicted Terminal DCM Comparison")
plt.legend()
plt.grid(True)
plt.show()
