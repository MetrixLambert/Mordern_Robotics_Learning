"""ROBO_Learn · Month 2 Day 20: Swing foot trajectory + planar leg IK.

Complete TODO 1–4, then run:
    python day20_swing_foot_ik_practice.py

Dependencies: numpy, matplotlib

Coordinates: world x = forward, z = upward; ground z=0.
Hip stays fixed at (hip_x, hip_z). The planar two-link leg's joint
angles follow the standard 2R arm convention, NOT anatomical angles.
This is a kinematic exercise: no ground contact forces or dynamics.
"""

from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np


@dataclass(frozen=True)
class Config:
    duration: float = 0.50  # seconds
    dt: float = 0.002       # seconds
    foot_start: tuple = (-0.12, 0.0)  # (x,z) in world [m]
    foot_end: tuple = (0.18, 0.0)
    clearance: float = 0.10  # [m]
    hip: tuple = (0.0, 0.75) # fixed hip, world [m]
    l1: float = 0.42         # thigh [m]
    l2: float = 0.42         # shank [m]


def quintic(u):
    """TODO 1: normalized quintic interpolation.

    Input u in [0,1], scalar or ndarray. Return alpha(u) satisfying
    alpha(0)=0, alpha(1)=1, alpha'(0)=alpha'(1)=0,
    alpha''(0)=alpha''(1)=0.
    Hint: 10*u**3 - 15*u**4 + 6*u**5.
    """
    return 10*u**3 - 15*u**4 + 6*u**5


def swing_foot(t, cfg):
    """TODO 2: return [x,z] in WORLD frame, for 0 <= t <= duration.

    Horizontal: interpolate foot_start[0] to foot_end[0] with quintic(t/T).
    Vertical: first half rises from ground to clearance with quintic;
              second half falls back to ground with another quintic.
    At t=T/2, z=clearance; endpoint velocities AND accelerations vanish.
    """
    x = cfg.foot_start[0] + (cfg.foot_end[0] - cfg.foot_start[0]) * quintic(t/cfg.duration)
    z = cfg.clearance * quintic(t/(cfg.duration/2)) if t <= cfg.duration/2 else cfg.clearance - cfg.clearance * quintic((t - cfg.duration/2)/(cfg.duration/2))
    return np.array([x, z])


def leg_ik(foot_world, cfg):
    """TODO 3: analytic 2R planar IK, returning [q1,q2] in radians.

    1. Convert foot_world to hip-relative vector [x,z].
    2. Solve cos(q2) by cosine law; use q2 = -arccos(cos_q2)
       consistently to choose one knee branch.
    3. Solve q1 using atan2 as in Month 1.

    IMPORTANT: Check reachability BEFORE arccos. Raise ValueError if
    distance > l1+l2 or < abs(l1-l2), within small tolerance.
    Do not hide an unreachable target by blindly clipping cos_q2.
    """
    x, z = foot_world - cfg.hip 
    r = np.sqrt(x**2 + z**2)
    if r > cfg.l1 + cfg.l2 or r < abs(cfg.l1 - cfg.l2):
        raise ValueError("Foot is out of reach")
    cos_q2 = (x**2 + z**2 - cfg.l1**2 - cfg.l2**2) / (2 * cfg.l1 * cfg.l2)
    if cos_q2 < -1 or cos_q2 > 1:
        raise ValueError("Invalid cosine value")
    q2 = -np.arccos(cos_q2)
    q1 = np.atan2(z, x) - np.atan2(cfg.l2 * np.sin(q2), cfg.l1 + cfg.l2 * np.cos(q2))
    return np.array([q1, q2])


def leg_fk(q, cfg):
    """TODO 4: return foot [x,z] in WORLD frame from [q1,q2].

    Standard planar 2R forward kinematics, then add world hip position.
    Use this INDEPENDENTLY of IK to validate the solution.
    """
    x = cfg.hip[0] + cfg.l1 * np.cos(q[0]) + cfg.l2 * np.cos(q[0] + q[1])
    z = cfg.hip[1] + cfg.l1 * np.sin(q[0]) + cfg.l2 * np.sin(q[0] + q[1])
    return np.array([x, z])


def main():
    cfg = Config()
    t = np.linspace(0, cfg.duration, round(cfg.duration / cfg.dt) + 1)

    # Validate interpolation boundary conditions.
    np.testing.assert_allclose([quintic(0.0), quintic(1.0)], [0, 1], atol=1e-12)
    foot = np.array([swing_foot(float(ti), cfg) for ti in t], dtype=float)
    assert foot.shape == (len(t), 2)
    assert np.all(np.isfinite(foot))
    np.testing.assert_allclose(foot[0], cfg.foot_start, atol=1e-10)
    np.testing.assert_allclose(foot[-1], cfg.foot_end, atol=1e-10)
    np.testing.assert_allclose(foot[len(t)//2, 1], cfg.clearance, atol=1e-10)
    assert np.min(foot[:, 1]) >= -1e-10, "Foot penetrates the ground."

    vel = np.gradient(foot, t, axis=0, edge_order=2)
    acc = np.gradient(vel, t, axis=0, edge_order=2)
    # Finite differences near endpoints are approximate; tolerate O(dt^2).
    assert np.max(np.abs(vel[[0, -1]])) < 2e-3, "Check endpoint velocity."
    # A finite-difference estimate of a zero endpoint acceleration has
    # O(dt) error when jerk is nonzero, especially for piecewise quintics.
    assert np.max(np.abs(acc[[0, -1]])) < 1.5, "Check endpoint acceleration."

    q = np.array([leg_ik(p, cfg) for p in foot], dtype=float)
    assert q.shape == (len(t), 2)
    foot_check = np.array([leg_fk(qi, cfg) for qi in q])
    err = np.linalg.norm(foot_check - foot, axis=1)
    print(f"Max FK(IK) position error: {np.max(err):.3e} m")
    np.testing.assert_allclose(foot_check, foot, atol=1e-9)

    fig, ax = plt.subplots(2, 2, figsize=(11, 7))
    ax[0, 0].plot(foot[:, 0], foot[:, 1], label="desired foot path")
    ax[0, 0].plot(foot_check[:, 0], foot_check[:, 1], "--", label="FK(IK)")
    ax[0, 0].set(xlabel="world x [m]", ylabel="world z [m]", title="Foot path")
    ax[0, 0].axis("equal")
    ax[0, 0].legend()

    for k, label in enumerate(["x", "z"]):
        ax[0, 1].plot(t, foot[:, k], label=f"{label} [m]")
        ax[1, 0].plot(t, vel[:, k], label=f"{label} velocity")
        ax[1, 1].plot(t, acc[:, k], label=f"{label} acceleration")
    ax[0, 1].set_title("Position")
    ax[1, 0].set_title("Velocity")
    ax[1, 1].set_title("Acceleration")
    for a in ax.flat:
        a.grid(alpha=0.3)
        a.legend()
    for a in ax.flat[1:]:
        a.set_xlabel("time [s]")
    fig.tight_layout()

    plt.figure(figsize=(8, 3))
    plt.plot(t, q[:, 0], label="hip angle q1")
    plt.plot(t, q[:, 1], label="knee angle q2")
    plt.xlabel("time [s]")
    plt.ylabel("angle [rad]")
    plt.title("Planar 2R inverse kinematics")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.show()

    print("\nQuestions:")
    print("1. Why should we check FK(IK(q)) instead of only plotting joint angles?")
    print("2. What changes if the hip moves during swing? Which frame is affected?")
    print("3. If the target is unreachable, why is clipping cos(q2) insufficient?")
    print("4. Why does smooth foot position not guarantee smooth joint torques?")


if __name__ == "__main__":
    main()
