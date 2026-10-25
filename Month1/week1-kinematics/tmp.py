import numpy as np 

def skew(w):
    """R^3 -> so(3)"""
    x, y, z = w
    return np.array([
        [0, -z, y],
        [z, 0, -x],
        [-y, x, 0]
    ], dtype=float)


def vee(W):
    """so(3) -> R^3"""
    return np.array([
        W[2, 1],
        W[0, 2],
        W[1, 0]
    ], dtype=float)

def adjoint(T):
    """
    Compute the 6x6 adjoint matrix of T in SE(3).

    T = [[R, p],
         [0, 1]]

    Convention: V = [omega, v]
    """
    R = T[:3, :3] 
    p = T[:3, 3]  
    Ad = np.zeros((6, 6)) 
    Ad[:3, :3] = R 
    Ad[3:, 3:] = R 
    Ad[3:, :3] = skew(p) @ R 
    return Ad 


R = np.eye(3) 
p = np.array([2, 0, 0])
T = np.eye(4)
T[:3, :3] = R
T[:3, 3] = p
print(T)
print(adjoint(T))


F_b = np.array([0, 0, 5, 0, 10, 0]) 
F_s = np.linalg.inv(adjoint(T)).T @ F_b 

print(F_s)









