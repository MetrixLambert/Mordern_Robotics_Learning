
import numpy as np

def rot2d(theta):
    # TODO: 根据旋转矩阵定义实现
    return np.array([
        [np.cos(theta), -np.sin(theta)],
        [np.sin(theta), np.cos(theta)]
    ])

theta = np.pi / 2
R = rot2d(theta)
print("R = \n", R)

# 1. 验证正交性
print("R * R^T = \n", R @ R.T)

# 2. 验证行列式
print("det(R) = ", np.linalg.det(R)) 

# 3. a -> b 坐标变换
p_a = np.array([1.0, 0.0])
p_b = R @ p_a 
print("p_b = \n", p_b)


# 4. b -> a 坐标变换
p_a2 = R.T @ p_b 
print("p_a2 = \n", p_a2)
