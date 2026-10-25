"""Semi-implicit Euler simulator: the plant always uses true dynamics."""
import numpy as np
from robot import forward_dynamics, joint_points
from controllers import controller_torque, saturate_torque
from trajectories import desired_state

def simulate(kind='ctc', task='tracking', model_scale=(1.,1.,1.),
             tau_limit=np.inf, kp=25., kd=10., duration=4., dt=0.002,
             initial_error=(0.15,-0.10)):
    """TODO 3: Integrate plant and collect histories.

    Return dict with keys t,q,qd,q_des,error,tau_cmd,tau_actual,
    saturated,rmse,final_error,peak_actual_torque,saturation_fraction,
    ee,ee_des. Shapes: t (n+1,), q/qd/q_des/error/ee/ee_des (n+1,2),
    tau_cmd/tau_actual/saturated (n,2).

    Steps at k: desired_state -> controller -> saturation ->
    forward_dynamics -> qd[k+1]=qd[k]+dt*qdd -> q[k+1]=q[k]+dt*qd[k+1].
    Start q[0]=reference(0)+initial_error; qd[0]=reference_velocity(0).
    RMSE = sqrt(mean(sum(error**2, axis=1))).
    Saturation fraction = mean(any(saturated,axis=1)).
    End effector = joint_points(q)[-1].
    """
    k = int(duration / dt) 
    t = np.linspace(0, duration, k+1) 

    qs = np.zeros((k+1, 2)) 
    qds = np.zeros((k+1, 2)) 
    qdds = np.zeros((k+1, 2)) 
    q_des_s = np.zeros((k+1, 2)) 

    # initial states 
    q_des_s[0], qds[0], qdds[0] = desired_state(0, task)
    qs[0] = q_des_s[0] + initial_error

    # control signals 
    tau_cmds = np.zeros((k, 2)) 
    tau_actuals = np.zeros((k, 2)) 
    saturations = np.zeros((k, 2)) 

    # end effector 
    ee_s = np.zeros((k+1, 2)) 
    ee_des_s = np.zeros((k+1, 2)) 

    # simulation loop 
    for i in range(k): 
        # desired state 
        qr, vr, ar = desired_state(t[i], task) 
        q_des_s[i+1] = qr 
        qds[i+1] = vr 
        qdds[i+1] = ar 

        # controller 
        tau_cmd = controller_torque(kind, qs[i], qds[i], qr, vr, ar, model_scale=model_scale)
        tau_actual, saturation_mask = saturate_torque(tau_cmd, tau_limit)
        tau_cmds[i] = tau_cmd 
        tau_actuals[i] = tau_actual 
        saturations[i] = saturation_mask 

        # forward dynamics 
        qdds[i+1] = forward_dynamics(qs[i], qds[i], tau_actual) 

        # update states 
        qds[i+1] = qds[i] + dt * qdds[i+1] 
        qs[i+1] = qs[i] + dt * qds[i+1] 

    # calculate metrics 
    error = qs - q_des_s 
    rmse = np.sqrt(np.mean(np.sum(error**2, axis=1)))
    final_error = error[-1] 
    peak_actual_torque = np.max(np.linalg.norm(tau_actuals, axis=1)) 
    saturation_fraction = np.mean(np.any(saturations, axis=1)) 

    
    # cal positions of end effector 
    for q in qs: 
        ee_s[i] = joint_points(q)[-1]
    for q in q_des_s: 
        ee_des_s[i] = joint_points(q)[-1]

    # return results 
    return {
        't': t, 
        'q': qs, 
        'qd': qds, 
        'q_des': q_des_s, 
        'error': qs - q_des_s, 
        'tau_cmd': tau_cmds, 
        'tau_actual': tau_actuals, 
        'saturated': saturations, 
        'rmse': rmse, 
        'final_error': final_error, 
        'peak_actual_torque': peak_actual_torque, 
        'saturation_fraction': saturation_fraction, 
        'ee': ee_s, 
        'ee_des': ee_des_s, 
    }