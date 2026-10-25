import numpy as np

def rot_z(theta):
    """返回绕 z 轴旋转 theta 的旋转矩阵"""
    # TODO
    return np.array([
        [np.cos(theta), -np.sin(theta), 0],
        [np.sin(theta), np.cos(theta), 0],
        [0, 0, 1]
    ])


def skew(w):
    """将三维向量 w 转换成反对称矩阵"""
    # TODO
    wx, wy, wz = w 
    return np.array([
        [0, -wz, wy],
        [wz, 0, -wx],
        [-wy, wx, 0]
    ])


def vee(W):
    """将 3x3 反对称矩阵恢复成三维向量"""
    # TODO
    wx, wy, wz = W[2,1], W[0,2], W[1,0]
    return np.array([wx, wy, wz])


omega = 2.0
t = 0.5
h = 1e-5

R = rot_z(omega * t)

# 中心差分
R_plus = rot_z(omega * (t + h))
R_minus = rot_z(omega * (t - h))
R_dot = (R_plus - R_minus) / (2 * h)

# 世界坐标系下的角速度
Omega_s = R_dot @ R.T
omega_s = vee(Omega_s)

# 刚体坐标系下的角速度
Omega_b = R.T @ R_dot
omega_b = vee(Omega_b)

print("omega_s =", omega_s)
print("omega_b =", omega_b)
print("Expected =", np.array([0, 0, omega]))