"""2R planar robot: true plant and simplified nominal model."""
import numpy as np

L1, L2 = 2.0, 1.0

def mass_matrix(q):
    q2 = float(q[1])
    return np.array([[6.60 + 2*np.cos(q2), 0.35 + np.cos(q2)],
                     [0.35 + np.cos(q2), 0.35]])

def coriolis_vector(q, qd):
    v1, v2 = qd
    return np.sin(q[1])*np.array([-2*v1*v2-v2*v2, v1*v1])

def gravity_vector(q):
    q1, q2 = q
    return 9.81*np.array([4*np.cos(q1)+0.5*np.cos(q1+q2),
                          0.5*np.cos(q1+q2)])

def model_terms(q, qd, scale=(1.,1.,1.)):
    a,b,d = scale
    return a*mass_matrix(q), b*coriolis_vector(q,qd), d*gravity_vector(q)

def joint_points(q):
    """Return (3,2) positions: base, elbow, end effector."""
    q1, q2 = q
    elbow = np.array([L1*np.cos(q1), L1*np.sin(q1)])
    tip = elbow + L2*np.array([np.cos(q1+q2), np.sin(q1+q2)])
    return np.vstack((np.zeros(2), elbow, tip))

def jacobian(q):
    q1,q2 = q
    return np.array([[-L1*np.sin(q1)-L2*np.sin(q1+q2), -L2*np.sin(q1+q2)],
                     [ L1*np.cos(q1)+L2*np.cos(q1+q2),  L2*np.cos(q1+q2)]])

def forward_dynamics(q, qd, tau):
    return np.linalg.solve(mass_matrix(q), tau-coriolis_vector(q,qd)-gravity_vector(q))
