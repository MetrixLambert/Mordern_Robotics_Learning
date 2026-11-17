"""Day 16 Part 3: Acrobot swing-up + LQR, tuning and comparisons.
Standalone: numpy scipy matplotlib.
Default tuned PFL gains (6,3,1) were tested at dt=0.02s, 18s, 50 Nm.
These gains are empirical and can be sensitive to sampling and initial state.
"""
import argparse
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.linalg import solve_continuous_are

G=9.81
ED=19.62

def wrap(a): return (a+np.pi)%(2*np.pi)-np.pi

def terms(x):
 q1,q2,d1,d2=x
 c=np.cos(q2);s=np.sin(q2)
 M=np.array([[3+np.cos(q2),1+.5*c],[1+.5*c,1.]])
 h=np.array([-.5*s*(2*d1*d2+d2*d2),.5*s*d1*d1])+G*np.array([1.5*np.sin(q1)+.5*np.sin(q1+q2),.5*np.sin(q1+q2)])
 return M,h

def energy(x):
 M,_=terms(x);q1,q2,d1,d2=x
 return .5*np.array([d1,d2])@M@np.array([d1,d2])-G*(1.5*np.cos(q1)+.5*np.cos(q1+q2))

def dynamics(x,tau):
 M,h=terms(x)
 dd=np.linalg.solve(M,np.array([0.,tau])-h)
 return np.array([x[2],x[3],*dd])

def lqr_gain():
 x=np.array([np.pi,0.,0.,0.]);eps=1e-6
 A=np.column_stack([(dynamics(x+np.eye(4)[i]*eps,0)-dynamics(x-np.eye(4)[i]*eps,0))/(2*eps) for i in range(4)])
 B=((dynamics(x,eps)-dynamics(x,-eps))/(2*eps)).reshape(-1,1)
 Q=np.diag([100.,100.,10.,10.]);R=np.array([[1.]])
 P=solve_continuous_are(A,B,Q,R)
 return np.linalg.solve(R,B.T@P).ravel()
K=lqr_gain()

def control(x,mode,gains,tau_max):
 M,h=terms(x)
 if mode=='pfl':
  k1,k2,k3=gains
  # wrapped q2 avoids unbounded position feedback over multiple revolutions
  v=-k1*wrap(x[1])-k2*x[3]+k3*x[2]*(energy(x)-ED)
  tau=(M[1,1]-M[1,0]*M[0,1]/M[0,0])*v+h[1]-M[1,0]/M[0,0]*h[0]
 elif mode=='direct':
  tau=-gains[0]*x[3]*(energy(x)-ED)
 elif mode=='lqr':
  err=np.array([wrap(x[0]-np.pi),wrap(x[1]),x[2],x[3]])
  tau=-K@err
 else: raise ValueError(mode)
 return float(np.clip(tau,-tau_max,tau_max))

def within(x,angle_deg,vel):
 return abs(wrap(x[0]-np.pi))<np.deg2rad(angle_deg) and abs(wrap(x[1]))<np.deg2rad(angle_deg) and max(abs(x[2]),abs(x[3]))<vel

def run(mode,gains,duration=18.,tau_max=50.,angle0=2.865,dt=.02,enter_angle=22.,enter_vel=3.,exit_angle=40.,exit_vel=5.):
 # Sample-and-hold control, mode updated only at sample boundaries (no solver-side effects).
 x=np.array([np.deg2rad(angle0),0.,0.,0.]);ts=[0.];xs=[x.copy()];us=[];ms=[];switch=[];state='swing'
 n=int(round(duration/dt))
 for i in range(n):
  t=i*dt
  if state=='swing' and within(x,enter_angle,enter_vel):
   state='lqr';switch.append((t,'enter'))
  elif state=='lqr' and not within(x,exit_angle,exit_vel):
   state='swing';switch.append((t,'exit'))
  u=control(x,'lqr' if state=='lqr' else mode,K if state=='lqr' else gains,tau_max)
  # Fixed-step RK4, constant torque during each controller interval.
  k1=dynamics(x,u);k2=dynamics(x+.5*dt*k1,u);k3=dynamics(x+.5*dt*k2,u);k4=dynamics(x+dt*k3,u)
  x=x+dt*(k1+2*k2+2*k3+k4)/6
  ts.append((i+1)*dt);xs.append(x.copy());us.append(u);ms.append(state)
 xs=np.array(xs);ts=np.array(ts);us=np.array(us)
 err=np.array([wrap(xs[-1,0]-np.pi),wrap(xs[-1,1]),xs[-1,2],xs[-1,3]])
 success=within(xs[-1],5.,.5) and all(within(z,7.,.7) for z in xs[-min(50,len(xs)):])
 return dict(t=ts,x=xs,u=us,m=ms,switch=switch,success=success,err=err,energies=np.array([energy(z) for z in xs]),lqr_time=sum(m=='lqr' for m in ms)*dt)

def score(res):
 x=res['x'];err=np.column_stack((wrap(x[:,0]-np.pi),wrap(x[:,1]),x[:,2],x[:,3]))
 # focus on late trajectory and successful capture
 last=err[int(.75*len(err)):]
 return float(np.mean(last[:,0]**2+last[:,1]**2+.07*(last[:,2]**2+last[:,3]**2)) + (0 if res['success'] else 5))

def tune(tau_max,angle0,dt,duration):
 candidates=[]
 # Coarse grid plus near-baseline settings; bounded search is empirical, not a proof.
 for k1 in [0.,1.,3.,6.,10.,15.]:
  for k2 in [0.,1.,3.,5.,8.]:
   for k3 in [.05,.15,.3,.6,1.,2.,4.]:
    gains=(k1,k2,k3)
    r=run('pfl',gains,duration,tau_max,angle0,dt)
    candidates.append((score(r),gains,r['success'],len(r['switch'])))
 candidates.sort(key=lambda v:v[0])
 return candidates

def plot(results,output):
 fig,axs=plt.subplots(4,1,figsize=(11,11),sharex=True)
 for name,r in results.items():
  t=r['t'];x=r['x']
  axs[0].plot(t,r['energies'],label=name)
  axs[1].plot(t,np.rad2deg(wrap(x[:,0]-np.pi)),label=name+' q1')
  axs[2].plot(t,np.rad2deg(wrap(x[:,1])),label=name+' q2')
  axs[3].plot(t[:-1],r['u'],label=name)
  for tt,_ in r['switch']:
   axs[1].axvline(tt,color='gray',alpha=.12)
 axs[0].axhline(ED,color='gray',ls='--',label='target')
 for ax in axs:ax.grid(alpha=.3);ax.legend(fontsize=8)
 axs[0].set_ylabel('Energy (J)');axs[1].set_ylabel('q1 error (deg)');axs[2].set_ylabel('q2 error (deg)');axs[3].set_ylabel('Torque (Nm)');axs[3].set_xlabel('Time (s)')
 fig.suptitle('Day 16 Part 3: hybrid swing-up / LQR');fig.tight_layout();fig.savefig(output,dpi=150);plt.close(fig)

def main():
 ap=argparse.ArgumentParser()
 ap.add_argument('--tune',action='store_true',help='Coarse search for PFL gains; may take a while')
 ap.add_argument('--duration',type=float,default=18.)
 ap.add_argument('--dt',type=float,default=.02)
 ap.add_argument('--tau-max',type=float,default=50.)
 ap.add_argument('--angle0',type=float,default=2.865)
 ap.add_argument('--pfl-gains',nargs=3,type=float,default=[6.,3.,1.],metavar=('K1','K2','K3'))
 ap.add_argument('--direct-gain',type=float,default=1.)
 ap.add_argument('--output',type=Path,default=Path('day16_part3_hybrid.png'))
 args=ap.parse_args()
 if args.duration<=0 or args.dt<=0 or args.tau_max<=0:
  ap.error('duration, dt and tau-max must be positive')
 if args.tune:
  top=tune(args.tau_max,args.angle0,args.dt,args.duration)
  print('Top PFL candidates (score, gains, captured, switches):')
  for item in top[:10]:print(item)
  gains=top[0][1]
 else:gains=tuple(args.pfl_gains)
 print('LQR K:',K,'PFL gains:',gains)
 results={
  'PFL + LQR':run('pfl',gains,args.duration,args.tau_max,args.angle0,args.dt),
  'Direct + LQR':run('direct',(args.direct_gain,),args.duration,args.tau_max,args.angle0,args.dt)
 }
 for name,r in results.items():
  print(f'{name}: capture={r["success"]}, switches={r["switch"][:15]}, final E={r["energies"][-1]:.3f}, final error={r["err"]}, LQR time={r["lqr_time"]:.2f}s')
 args.output.parent.mkdir(parents=True,exist_ok=True)
 plot(results,args.output)
 print('Saved',args.output.resolve())
if __name__=='__main__':main()
