"""Run after completing TODO 1-3. python main.py"""
import numpy as np 
import matplotlib.pyplot as plt
from simulator import simulate
from visualization import plot_results, animate_robot

if __name__=='__main__':
    r=simulate(kind='ctc',tau_limit=[45.,8.])
    print(f'RMSE={r["rmse"]:.5f}, final={np.linalg.norm(r["final_error"]):.5f}, sat={r["saturation_fraction"]:.1%}')
    plot_results(r)
    fig,anim=animate_robot(r)
    plt.show()
