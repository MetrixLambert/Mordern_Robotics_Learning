"""Day 20 — ISS and robust stability (NumPy only).

Complete TODO 1–4, then run:
    python day20_iss_robustness_exercises.py

Model:
    e_dd + kd*e_d + kp*e = delta_p*e + delta_d*e_d + d(t)
Equivalently:
    e_dd + (kd-delta_d)*e_d + (kp-delta_p)*e = d(t).

This script intentionally uses no scipy and no plotting dependencies.
"""
import numpy as np


# TODO 1 --------------------------------------------------------------------
def closed_loop_coefficients(kp, kd, delta_p, delta_d):
    """Return (kp_eff, kd_eff) after moving state-dependent terms to LHS."""
    return (kp - delta_p, kd - delta_d)


# TODO 2 --------------------------------------------------------------------
def is_iss(kp, kd, delta_p, delta_d):
    """Return bool: is the scalar LTI error system ISS w.r.t. d(t)?

    Apply the second-order Hurwitz criterion (strict inequalities).
    """
    return (kp - delta_p) > 0 and (kd - delta_d) > 0


# TODO 3 --------------------------------------------------------------------
def constant_disturbance_equilibrium(kp, kd, delta_p, delta_d, d0):
    """Return steady-state position error under constant d0.

    Raise ValueError when the closed loop is not ISS: a unique attracting
    equilibrium is not guaranteed then.
    """
    if not is_iss(kp, kd, delta_p, delta_d):
        raise ValueError("The closed loop is not ISS")
    return d0 / (kp - delta_p)


def disturbance(t, kind="constant", amplitude=0.2, omega=3.0):
    """Provided disturbance signals (scalar acceleration units)."""
    if kind == "zero":
        return 0.0
    if kind == "constant":
        return float(amplitude)
    if kind == "sinusoid":
        return float(amplitude * np.sin(omega * t))
    raise ValueError("kind must be zero, constant, or sinusoid")


def simulate(kp, kd, delta_p, delta_d, *, kind="constant",
             amplitude=0.2, omega=3.0, initial_error=0.3,
             initial_velocity=0.0, duration=12.0, dt=0.001):
    """Provided semi-implicit Euler simulation of the true error dynamics.

    Caution: unstable cases can grow large; duration is deliberately finite.
    """
    assert kp > 0 and kd > 0 and dt > 0 and duration > 0
    n = int(round(duration / dt))
    t = np.arange(n + 1) * dt
    e = np.zeros(n + 1)
    ed = np.zeros(n + 1)
    d = np.zeros(n + 1)
    e[0], ed[0] = initial_error, initial_velocity
    for k in range(n):
        d[k] = disturbance(t[k], kind, amplitude, omega)
        # True dynamics, independent of TODO implementations:
        edd = -(kd - delta_d) * ed[k] - (kp - delta_p) * e[k] + d[k]
        ed[k + 1] = ed[k] + dt * edd
        e[k + 1] = e[k] + dt * ed[k + 1]
    d[-1] = disturbance(t[-1], kind, amplitude, omega)
    return {"t": t, "e": e, "ed": ed, "d": d}


# TODO 4 --------------------------------------------------------------------
def compare_robustness_cases(kp=25.0, kd=10.0):
    """Return a dict with exactly these six simulation outputs:

    'nominal'       delta_p=0,  delta_d=0,  constant d0=0.2
    'near_p_limit'  delta_p=24, delta_d=0,  constant d0=0.2
    'p_unstable'    delta_p=30, delta_d=0,  zero disturbance
    'near_d_limit'  delta_p=0,  delta_d=9,  zero disturbance
    'd_unstable'    delta_p=0,  delta_d=12, zero disturbance
    'periodic'      delta_p=5,  delta_d=2,  sinusoidal d(t)=0.2*sin(3t)

    For every case, use simulate() defaults except the listed parameters.
    """
    results = {}
    results['nominal'] = simulate(kp, kd, 0, 0, amplitude=0.2)
    results['near_p_limit'] = simulate(kp, kd, 24, 0, amplitude=0.2)
    results['p_unstable'] = simulate(kp, kd, 30, 0, amplitude=0.2)
    results['near_d_limit'] = simulate(kp, kd, 0, 9, amplitude=0.2)
    results['d_unstable'] = simulate(kp, kd, 0, 12, amplitude=0.2)
    results['periodic'] = simulate(kp, kd, 5, 2, kind="sinusoid", amplitude=0.2, omega=3.0)
    return results


def _close(actual, expected, tol=1e-10):
    assert np.isclose(actual, expected, atol=tol, rtol=tol), (actual, expected)


def run_tests():
    kp, kd = 25.0, 10.0
    print("[1] Effective stiffness and damping")
    _close(closed_loop_coefficients(kp, kd, 4.0, 2.0)[0], 21.0)
    _close(closed_loop_coefficients(kp, kd, 4.0, 2.0)[1], 8.0)
    print("  PASS")

    print("[2] ISS / Hurwitz boundary")
    for dp, dd, expected in [(0, 0, True), (24, 9, True),
                             (25, 0, False), (0, 10, False),
                             (30, 0, False), (0, 12, False),
                             (-10, -10, True)]:
        assert bool(is_iss(kp, kd, dp, dd)) == expected, (dp, dd)
    print("  PASS")

    print("[3] Constant-disturbance equilibrium")
    _close(constant_disturbance_equilibrium(kp, kd, 0, 0, 0.2), 0.008)
    _close(constant_disturbance_equilibrium(kp, kd, 24, 0, 0.2), 0.2)
    try:
        constant_disturbance_equilibrium(kp, kd, 25, 0, 0.2)
    except ValueError:
        pass
    else:
        raise AssertionError("Non-ISS case should raise ValueError")
    print("  PASS")

    print("[4] Six robustness experiments")
    results = compare_robustness_cases(kp, kd)
    keys = {'nominal', 'near_p_limit', 'p_unstable',
            'near_d_limit', 'd_unstable', 'periodic'}
    assert set(results) == keys
    for name, data in results.items():
        assert set(data) == {'t', 'e', 'ed', 'd'}
        assert len(data['t']) == 12001
        assert np.all(np.isfinite(data['e']))
        print(f"  {name:14s} final_e={data['e'][-1]:+.6f}, "
              f"max_abs_e={np.max(np.abs(data['e'])):.6f}")
    _close(results['nominal']['e'][-1], 0.008, tol=2e-4)
    # Near stiffness boundary, slow convergence; equilibrium is 0.2,
    # but at t=12 the transient has not fully vanished.
    assert 0.15 < results['near_p_limit']['e'][-1] < 0.3
    assert abs(results['p_unstable']['e'][-1]) > 1.0
    assert np.max(np.abs(results['d_unstable']['e'])) > 1.0
    print("  PASS")

    print("[5] Bounded sinusoidal input: finite response")
    periodic = results['periodic']
    assert np.max(np.abs(periodic['d'])) <= 0.2 + 1e-12
    assert np.max(np.abs(periodic['e'])) < 1.0
    print("  PASS")
    print("\nDay 20: all tests passed.")


if __name__ == '__main__':
    run_tests()
