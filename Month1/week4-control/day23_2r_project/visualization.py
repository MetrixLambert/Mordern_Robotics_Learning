"""Plot and animate a 2R simulation (matplotlib required)."""
import numpy as np
from robot import joint_points

def plot_results(result):
    import matplotlib.pyplot as plt
    t=result['t']; q=result['q']; qr=result['q_des']; e=result['error']
    fig, axes=plt.subplots(2,2,figsize=(11,7))
    for j in range(2):
        axes[0,0].plot(t,q[:,j],label=f'q{j+1}')
        axes[0,0].plot(t,qr[:,j],'--',label=f'q{j+1} desired')
        axes[0,1].plot(t,e[:,j],label=f'e{j+1}')
        axes[1,0].plot(t[:-1],result['tau_actual'][:,j],label=f'tau{j+1}')
    axes[1,1].plot(*result['ee'].T,label='actual')
    axes[1,1].plot(*result['ee_des'].T,'--',label='desired')
    for ax in axes.flat: ax.grid(True); ax.legend()
    axes[0,0].set_title('Joint positions'); axes[0,1].set_title('Tracking errors')
    axes[1,0].set_title('Actual torque (Nm)'); axes[1,1].set_title('End-effector path')
    axes[1,1].set_aspect('equal',adjustable='box')
    fig.tight_layout(); return fig

def animate_robot(result, stride=20, interval=40):
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    fig,ax=plt.subplots(figsize=(6,6))
    ax.set(xlim=(-3.2,3.2),ylim=(-3.2,3.2),xlabel='x (m)',ylabel='y (m)')
    ax.set_aspect('equal'); ax.grid(True)
    line,=ax.plot([],[],'o-',lw=3,label='actual')
    desired,=ax.plot([],[],'o--',alpha=0.5,label='desired')
    ax.legend()
    ids=np.arange(0,len(result['t']),stride)
    def update(i):
        k=ids[i]
        a=joint_points(result['q'][k]); b=joint_points(result['q_des'][k])
        line.set_data(a[:,0],a[:,1]); desired.set_data(b[:,0],b[:,1])
        ax.set_title(f't={result["t"][k]:.2f}s')
        return line,desired
    anim=FuncAnimation(fig,update,frames=len(ids),interval=interval,blit=False)
    anim.save('animation.gif', writer='pillow', fps=20)
    return fig,anim
