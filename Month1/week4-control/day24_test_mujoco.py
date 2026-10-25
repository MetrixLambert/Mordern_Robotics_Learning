import numpy as np
import mujoco


# ========================================
# 1. 创建机器人模型
# ========================================

XML = """
<mujoco model="two_link_arm">
    <option timestep="0.002" gravity="0 0 -9.81"/>

    <worldbody>
        <body name="link1" pos="0 0 0">
            <joint name="joint1" type="hinge" axis="0 1 0"/>

            <geom type="capsule"
                  fromto="0 0 0 0 0 -1"
                  size="0.05"
                  mass="2"/>

            <body name="link2" pos="0 0 -1">
                <joint name="joint2" type="hinge" axis="0 1 0"/>

                <geom type="capsule"
                      fromto="0 0 0 0 0 -1"
                      size="0.05"
                      mass="1"/>
            </body>
        </body>
    </worldbody>

    <actuator>
        <motor joint="joint1" gear="1"/>
        <motor joint="joint2" gear="1"/>
    </actuator>
</mujoco>
"""

model = mujoco.MjModel.from_xml_string(XML)
data = mujoco.MjData(model)


print(mujoco.__version__)
print([name for name in dir(data) if "mass" in name.lower() or "qm" in name.lower()])

# ========================================
# 2. 设置初始状态
# ========================================

data.qpos[:] = [0.3, -0.5]
data.qvel[:] = [0.2, -0.1]
data.ctrl[:] = [5.0, -2.0]

# 更新当前状态下的动力学量
mujoco.mj_forward(model, data)


# ========================================
# 3. 读取动力学参数
# ========================================

M = np.zeros((model.nv, model.nv), dtype=np.float64)
mujoco.mj_fullM(model, data, M)
print("M =")
print(M)


# TODO 2:
# 读取偏置广义力 b = c + g
b = data.qfrc_bias.copy()


# TODO 3:
# 读取执行器产生的广义力 tau
# 提示：这里可以使用 data.qfrc_actuator
tau = data.qfrc_actuator.copy()


# ========================================
# 4. NumPy 正向动力学
# ========================================

# TODO 4:
# 根据 M qdd + b = tau 计算加速度
qdd_numpy = np.linalg.solve(M, tau - b) 


# ========================================
# 5. MuJoCo 正向动力学
# ========================================

qdd_mujoco = data.qacc.copy()


# ========================================
# 6. 比较结果
# ========================================

print("M =")
print(M)

print("\nb =", b)
print("tau =", tau)

print("\nqdd_numpy =", qdd_numpy)
print("qdd_mujoco =", qdd_mujoco)

# TODO 5:
# 计算两种加速度的无穷范数误差
error = np.linalg.norm(qdd_numpy - qdd_mujoco, np.inf)

print("\nerror =", error)

assert error < 1e-8
print("PASS: forward dynamics consistency")

def computed_torque_control(
    model, data,
    q_des, qd_des, qdd_des,
    kp, kd
):
    # TODO 1: 获取完整惯量矩阵
    M = np.zeros((model.nv, model.nv), dtype=np.float64)
    mujoco.mj_fullM(model, data, M) 

    # TODO 2: 计算位置误差和速度误差
    e = data.qpos - q_des 
    ed = data.qvel - qd_des

    # TODO 3: 计算期望的闭环加速度 v
    v = qdd_des - kp * e - kd * ed

    # TODO 4: 计算关节力矩 tau
    tau = M @ v + data.qfrc_bias 

    return tau

# ========================================
# Closed-loop CTC simulation
# ========================================

mujoco.mj_resetData(model, data)

data.qpos[:] = [0.0, 0.0]
data.qvel[:] = [0.0, 0.0]

q_des = np.array([0.5, -0.3])
qd_des = np.zeros(2)
qdd_des = np.zeros(2)

kp = 25.0
kd = 10.0

dt = model.opt.timestep
T = 3.0
steps = int(T / dt)

time_hist = []
error_hist = []

for k in range(steps):
    # 更新当前状态下的动力学量
    mujoco.mj_forward(model, data)

    # TODO 1: 调用你实现的 CTC 控制器
    tau = computed_torque_control(model, data, q_des, qd_des, qdd_des, kp, kd)

    # TODO 2: 将 tau 写入 MuJoCo 执行器
    data.ctrl[:] = tau

    # TODO 3: 推进一步仿真
    mujoco.mj_step(model, data) 

    # 记录位置误差
    e = data.qpos - q_des

    time_hist.append(data.time)
    error_hist.append(e.copy())

error_hist = np.array(error_hist)

print("Final position error:", error_hist[-1])
print("Final error norm:", np.linalg.norm(error_hist[-1]))

import matplotlib.pyplot as plt

time_hist = np.asarray(time_hist)
error_hist = np.asarray(error_hist)

e0 = np.array([-0.5, 0.3])

# TODO 1：根据临界阻尼解析解计算整个理论误差轨迹
error_theory = ((1 + 5 * time_hist) * np.exp(-5 * time_hist))[:, None] * e0[None, :]

# TODO 2：计算整个时间序列的最大误差范数
max_deviation = np.max(np.linalg.norm(error_hist - error_theory, axis=1))

fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)

for i in range(2):
    axes[i].plot(time_hist, error_hist[:, i], label="MuJoCo")
    axes[i].plot(
        time_hist, error_theory[:, i],
        "--", label="Analytical"
    )
    axes[i].set_ylabel(f"e{i+1} (rad)")
    axes[i].grid(True)
    axes[i].legend()

axes[-1].set_xlabel("Time (s)")
plt.tight_layout()
plt.show()

print("Maximum trajectory deviation:", max_deviation)