import numpy as np 

Mq = np.array([
    [3, 1], 
    [1, 2]
])

q_dd = np.array([0, 1]) 

tau = Mq @ q_dd 
print(tau) 

tau = np.array([0, 2]) 
q_dd = np.linalg.pinv(Mq) @ tau 
print(q_dd) 









