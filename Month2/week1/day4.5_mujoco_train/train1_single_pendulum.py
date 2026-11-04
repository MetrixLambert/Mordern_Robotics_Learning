# 这个文件是用于熟悉/自主完成mujoco框架的模型搭建和仿真

import numpy as np 
import mujoco 
import matplotlib.pyplot as plt 


model_XML = """
<mujoco model="simple_pendulum"> 
    <option timestep="0.002" gravity="0 0 -9.81" /> 
    <worldbody> 
        <body name="pendulum" pos="0 0 0"> 
            <joint  name="joint1" 
                    type="hinge"
                    axis="0 1 0" /> 
            <geom   type="capsule" 
                    fromto="0 0 0 0 0 -1"
                    size="0.04"
                    mass="1" /> 
        </body> 
    </worldbody> 
    <actuator> 
        <motor name="motor1" joint="joint1" /> 
    </actuator> 
</mujoco> 
"""


model = mujoco.MjModel.from_xml_string(model_XML) 
data = mujoco.MjData(model) 


# 仿真参数
sim_time = 5  # simulate for 5 seconds 
qpos_init = np.array([0.8]) 
q_des = np.array([0.5])
qvel_init = np.array([0.0]) 

Kp = 20 
Kd = 5 


data.qpos[:] = qpos_init 
data.qvel[:] = qvel_init  

times = [] 
angles = [] 
vels = [] 
torques = [] 


# 仿真循环
while data.time < sim_time: 
    tau = Kp * (q_des - data.qpos[0]) - Kd * data.qvel[0]  + 9.81 * 0.5 * np.sin(data.qpos[0]) 
    data.ctrl[0] = tau[0]  

    times.append(data.time) 
    angles.append(data.qpos[0]) 
    vels.append(data.qvel[0]) 
    torques.append(tau) 

    mujoco.mj_step(model, data) 


# 绘制图像
plt.figure(figsize=(9, 5)) 
plt.plot(times, angles, label="Angle (rad)")
plt.plot(times, vels, label="Velocity (rad/s)")
plt.xlabel("Time (s)")
plt.ylabel("State") 
plt.grid() 
plt.legend() 
plt.tight_layout() 
plt.show() 


print("nq = ", model.nq) 
print("nv = ", model.nv)  
print("nu = ", model.nu) 
print("Final angle = ", data.qpos[0])

