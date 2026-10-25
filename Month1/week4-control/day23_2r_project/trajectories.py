"""Desired joint-space motion."""
import numpy as np

def desired_state(t, task='tracking'):
    goal=np.array([0.8,-0.6])
    if task == 'regulation':
        return goal, np.zeros(2), np.zeros(2)
    if task != 'tracking':
        raise ValueError(task)
    T=2.0
    s=np.clip(t/T,0.,1.)
    h=10*s**3-15*s**4+6*s**5
    hd=(30*s**2-60*s**3+30*s**4)/T if t<T else 0.
    hdd=(60*s-180*s**2+120*s**3)/T**2 if t<T else 0.
    return goal*h, goal*hd, goal*hdd
