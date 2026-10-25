"""Implement four controllers and componentwise total-torque saturation."""
import numpy as np
from robot import model_terms

def controller_torque(kind, q, qd, qr, vr, ar, kp=25., kd=10., model_scale=(1.,1.,1.)):
    """TODO 1: PD, GPD, FFPD, CTC. Same formulas as Day 22.

    Use model_terms(q,qd,model_scale) for current-state dynamics;
    use model_terms(qr,vr,model_scale) for reference-state feedforward.
    Raise ValueError for unknown kind.
    """
    e = q - qr 
    ed = qd - vr 
    if kind == "pd": 
        return - kp * e - kd * ed 
    elif kind == 'gpd': 
        M_hat, c_hat, g_hat = model_terms(q, qd, model_scale)
        return g_hat - kp * e - kd * ed  
    elif kind == 'ffpd': 
        M_hat, c_hat, g_hat = model_terms(q, qd, model_scale)
        return M_hat @ ar + c_hat + g_hat - kp * e - kd * ed   
    elif kind == 'ctc': 
        M_hat, c_hat, g_hat = model_terms(q, qd, model_scale) 
        v = ar - kd * ed - kp * e 
        return M_hat @ v + c_hat + g_hat  
    else: 
        raise ValueError(f"Unknown controller kind: {kind}")

def saturate_torque(tau_cmd, tau_limit):
    """TODO 2: Return (tau_actual, saturation_mask), not delta_tau.

    tau_limit may be scalar, vector, or np.inf. Never mutate tau_cmd.
    """
    tau_limit = np.asarray(tau_limit) 
    tau_actual = np.clip(tau_cmd, -tau_limit, tau_limit) 
    saturation_mask = tau_cmd != tau_actual  
    return tau_actual, saturation_mask 
