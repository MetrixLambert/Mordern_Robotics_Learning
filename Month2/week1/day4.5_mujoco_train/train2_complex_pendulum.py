import numpy as np 
import time 
import mujoco 
import mujoco.viewer 
import matplotlib.pyplot as plt 


model_XML = """
<mujoco model="complex_pendulum"> 
    <option timestep="0.002" gravity="0 0 -9.81" /> 
    <worldbody> 
        <body name="link1" pos="0 0 0"> 
            <joint  name="joint1" 
                    type="hinge" 
                    axis="0 1 0" /> 
            <geom   type="capsule" 
                    fromto="0 0 0 0 0 -1" 
                    size="0.04" 
                    mass = "1" /> 
            <body name="link2" pos="0 0 -1"> 
                <joint  name="joint2"
                        type="hinge" 
                        axis="0 1 0" /> 
                <geom   type="capsule" 
                        fromto="0 0 0 0 0 -0.8" 
                        size="0.04" 
                        mass="0.8" /> 
            </body> 
        </body> 
    </worldbody> 
    <actuator> 
        <motor name="motor1" joint="joint1"/> 
        <motor name="motor2" joint="joint2"/> 
    </actuator> 
</mujoco> 
"""

model = mujoco.MjModel.from_xml_string(model_XML)  
data = mujoco.MjData(model) 


# 仿真参数
sim_time = 5 
qpos_init = np.array([0.5, -0.3])
qvel_init = np.array([0.0, 0.0])
data.qpos[:] = qpos_init 
data.qvel[:] = qvel_init 

qpos_des = np.array([0, 0.5]) 

Kp = np.array([30.0, 30.0]) 
Kd = np.array([8.0, 8.0]) 

times = []
angles = [] 
vels = [] 
torques = [] 


# 仿真循环
with mujoco.viewer.launch_passive(model, data) as viewer: 
    while data.time < sim_time and viewer.is_running(): 
        start_time = time.time() 

        tau = Kp * (qpos_des - data.qpos) - Kd * data.qvel 
        data.ctrl[:] = tau 

        times.append(data.time) 
        angles.append(data.qpos.copy()) 
        vels.append(data.qvel.copy()) 
        torques.append(data.ctrl.copy()) 

        mujoco.mj_step(model, data) 

        viewer.sync() 

        elapsed_time = time.time() - start_time 
        remaining = model.opt.timestep - elapsed_time 
        if remaining > 0: 
            time.sleep(remaining) 

    while viewer.is_running(): 
        time.sleep(0.05) 

# 绘制图像
times = np.array(times) 
angles = np.array(angles) 
vels = np.array(vels) 
torques = np.array(torques) 

plt.figure(figsize=(9, 5)) 
plt.plot(times, angles[:, 0], label="Angle 1 (rad)") 
plt.plot(times, angles[:, 1], label="Angle 2 (rad)")
plt.xlabel("Time (s)")
plt.ylabel("Angle (rad)")
plt.grid() 
plt.legend() 
plt.tight_layout() 

plt.figure(figsize=(9, 5)) 
plt.plot(times, vels[:, 0], label="Velocity 1 (rad/s)") 
plt.plot(times, vels[:, 1], label="Velocity 2 (rad/s)") 
plt.xlabel("Time (s)")
plt.ylabel("Velocity (rad/s)")
plt.grid() 
plt.legend() 
plt.tight_layout() 

plt.figure(figsize=(9, 5)) 
plt.plot(times, torques[:, 0], label="Torque 1 (Nm)") 
plt.plot(times, torques[:, 1], label="Torque 2 (Nm)") 
plt.xlabel("Time (s)")
plt.ylabel("Torque (Nm)") 
plt.grid() 
plt.legend() 
plt.tight_layout() 

print("nq = ", model.nq) 
print("nv = ", model.nv) 
print("nu = ", model.nu) 
print("Final angles = ", data.qpos) 
print("Final velocities = ", data.qvel) 

plt.show() 

