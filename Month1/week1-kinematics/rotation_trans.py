import numpy as np

def skew(w):
    wx, wy, wz = w
    return np.array([
        [0, -wz, wy],
        [wz, 0, -wx],
        [-wy, wx, 0]
    ])


def so3_exp(phi):
    """
    SO(3) 指数映射

    Parameters
    ----------
    phi : ndarray, shape (3,)
        旋转向量，模长为旋转角度（rad）

    Returns
    -------
    R : ndarray, shape (3, 3)
        旋转矩阵
    """
    # TODO: 计算 theta
    theta = np.linalg.norm(phi) 

    # TODO: 处理 theta = 0
    if theta == 0:
        return np.eye(3)
    
    # TODO: 计算单位旋转轴
    omega = phi / theta 


    # TODO: 利用 Rodrigues 公式计算 R
    K = skew(omega) 
    R = np.eye(3) + np.sin(theta) * K + (1 - np.cos(theta)) * K @ K 

    return R 

def vee(W):
    """将反对称矩阵转换成三维向量"""
    # TODO
    wx, wy, wz = W[2,1], W[0,2], W[1,0]
    return np.array([wx, wy, wz])


def so3_log(R):
    """
    SO(3) 对数映射（一般情况）

    Parameters
    ----------
    R : ndarray, shape (3, 3)
        旋转矩阵，假设 0 < theta < pi

    Returns
    -------
    phi : ndarray, shape (3,)
        主值旋转向量
    """
    # TODO 1: 利用 trace(R) 计算 cos(theta)
    cos_theta = (np.trace(R) - 1) / 2 
    cos_theta = np.clip(cos_theta, -1, 1)


    # TODO 2: 计算 theta
    theta = np.arccos(cos_theta) 
    if theta < 1e-8: 
        return np.zeros(3) 
    elif np.isclose(theta, np.pi):
        # 求解 R 对应特征值 1 的特征向量
        # 将其归一化
        # 返回 theta * omega
        eig_vals, eig_vecs = np.linalg.eig(R)
        idx = np.argmin(np.abs(eig_vals - 1))
        omega = eig_vecs[:, idx]
        omega = omega / np.linalg.norm(omega)
        return theta * omega 

    # TODO 3: 利用 R - R.T 恢复旋转向量
    phi = vee(R - R.T) / (2 * np.sin(theta)) 

    return theta * phi 


# Test 1: 绕 z 轴旋转 90°
phi1 = np.array([0.0, 0.0, np.pi / 2])
R1 = so3_exp(phi1)

R_expected = np.array([
    [0, -1, 0],
    [1,  0, 0],
    [0,  0, 1]
])

print("Test 1:", np.allclose(R1, R_expected))

# Test 2: 任意旋转轴
phi2 = np.array([0.3, -0.5, 0.8])
R2 = so3_exp(phi2)

print("Test 2 (orthogonality):",
      np.allclose(R2.T @ R2, np.eye(3)))

print("Test 2 (determinant):",
      np.linalg.det(R2))

# Test 3: 零旋转
phi3 = np.zeros(3)
R3 = so3_exp(phi3)

print("Test 3:", np.allclose(R3, np.eye(3)))


# Test 4: 验证对数映射
phi_original = np.array([0.3, -0.5, 0.8])

R = so3_exp(phi_original)
phi_recovered = so3_log(R)

print("Original:", phi_original)
print("Recovered:", phi_recovered)
print("Error:", np.linalg.norm(phi_original - phi_recovered))

# Test 5: 验证对数映射（特殊情况）
test_phis = [
    np.array([0.0, 0.0, 0.0]),
    np.array([0.3, -0.5, 0.8]),
    np.array([0.0, 0.0, np.pi]),
    np.array([np.pi / np.sqrt(2),
              np.pi / np.sqrt(2),
              0.0])
]

for phi in test_phis:
    R = so3_exp(phi)
    phi_recovered = so3_log(R)
    R_recovered = so3_exp(phi_recovered)

    print("Original phi:", phi)
    print("Recovered phi:", phi_recovered)
    print("Rotation error:",
          np.linalg.norm(R - R_recovered))
    print()