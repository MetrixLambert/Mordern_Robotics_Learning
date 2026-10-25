"""Run: python test_robot.py"""
import numpy as np
from robot import mass_matrix, jacobian, joint_points, forward_dynamics
from controllers import controller_torque, saturate_torque
from simulator import simulate
from trajectories import desired_state

def test_kinematics():
    assert np.allclose(joint_points([0.,0.])[-1],[3.,0.])
    q=np.array([.3,-.4]); h=1e-6
    for i in range(2):
        d=np.zeros(2); d[i]=h
        numerical=(joint_points(q+d)[-1]-joint_points(q-d)[-1])/(2*h)
        assert np.allclose(numerical,jacobian(q)[:,i],atol=1e-8)
    print('[1] Kinematics PASS')

def test_dynamics():
    q=np.array([.2,-.3]); M=mass_matrix(q)
    assert np.allclose(M,M.T)
    assert np.all(np.linalg.eigvalsh(M)>0)
    tau=np.array([2.,3.]); qd=np.array([.1,-.2])
    from robot import coriolis_vector, gravity_vector
    a=forward_dynamics(q,qd,tau)
    assert np.allclose(M@a+coriolis_vector(q,qd)+gravity_vector(q),tau)
    print('[2] Dynamics PASS')

def test_controllers():
    q=np.array([.2,-.3]); qd=np.array([.1,-.2]); qr=np.array([.4,-.5])
    vr=np.array([.2,.1]); ar=np.array([.3,-.1]); e=q-qr; ed=qd-vr
    from robot import model_terms
    M,c,g=model_terms(q,qd)
    assert np.allclose(controller_torque('pd',q,qd,qr,vr,ar),-25*e-10*ed)
    assert np.allclose(controller_torque('gpd',q,qd,qr,vr,ar),g-25*e-10*ed)
    assert np.allclose(controller_torque('ctc',q,qd,qr,vr,ar),M@(ar-10*ed-25*e)+c+g)
    actual,mask=saturate_torque(np.array([60.,-4.]),[45.,8.])
    assert np.allclose(actual,[45.,-4.])
    assert np.array_equal(mask,[True,False])
    print('[3] Controllers and saturation PASS')

def test_simulation():
    r=simulate('ctc')
    assert r['q'].shape==(2001,2)
    assert r['tau_actual'].shape==(2000,2)
    assert r['ee'].shape==(2001,2)
    assert np.isfinite(r['rmse'])
    assert abs(r['rmse']-.04506)<.003
    s=simulate('ctc',tau_limit=[45.,8.])
    assert abs(s['rmse']-1.09585)<.03
    assert abs(s['saturation_fraction']-.54)<.03
    print('[4] Simulator and Day 22 regression PASS')

if __name__=='__main__':
    test_kinematics();test_dynamics();test_controllers();test_simulation()
    print('Day 23: all tests passed.')
