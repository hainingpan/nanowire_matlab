import matplotlib
matplotlib.use('Agg')
import numpy as np
import tinyarray
import matplotlib.pyplot as plt
from scipy.sparse import spdiags
from scipy.sparse import eye
from scipy.sparse import kron
from scipy.sparse.linalg import inv
from scipy.sparse import csr_matrix
import adaptive
from functools import partial
from scipy.interpolate import griddata
import sys
from mpi4py.futures import MPIPoolExecutor


s0 = tinyarray.array([[1, 0], [0, 1]]);
sx = tinyarray.array([[0, 1], [1, 0]]);
sy = tinyarray.array([[0, -1j], [1j, 0]]);
sz = tinyarray.array([[1, 0], [0, -1]]);

def hdis(a,mu,delta,vz,alpha_R,dim,vimp):
    t=25/a**2
    alpha=alpha_R/(2*a)
    band11sm=spdiags(np.vstack([np.ones(dim),np.ones(dim)]),np.array([-1,1]),dim,dim,format = 'csr')
    band1m1sm=spdiags(np.vstack([np.ones(dim),-np.ones(dim)]),np.array([-1,1]),dim,dim,format = 'csr')
    eyesm=eye(dim)
    mulist=mu*np.ones(dim)-vimp
    diagmulist=spdiags(mulist,0,dim,dim)
    return kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1j*alpha*band1m1sm)))\
    +kron(eye(2),kron(sz,vz*eyesm))+kron(sx,kron(eye(2),delta*eyesm))

def ldosall_dis(a,mu,Delta,Vz,alpha_R,mulist,dim,omega,delta):
    ham=hdis(a,mu,Delta,Vz,alpha_R,dim,mulist);
    hh=csr_matrix((omega+1j*delta)*eye(4*dim)-ham)
    G=inv(hh)
    Gdiag=(G).diagonal()
    return -np.sum((np.reshape(Gdiag,(4,-1))),0).imag/np.pi

def LDOS_dis(p,a,mu,Delta,alpha_R,mulist,dim,delta):
    Vz,energy=p
    z=ldosall_dis(a,mu,Delta,Vz,alpha_R,mulist,dim,energy,delta)
    return np.array([z.mean(),z[0],z[int(dim/2)],z[-1]])

def main():
    vars=len(sys.argv)
    if vars>1:
        loss=float(sys.argv[1])
    else:
        loss=0.1
    fname='savloss'+str(loss)
    learner = adaptive.Learner2D(partial(LDOS_dis,a=1,mu=1,Delta=0.2,alpha_R=5,mulist=0,dim=1000,delta=1e-3),\
                                 bounds=[(0., 2.), (-0.3, 0.3)])
    learner.load(fname)
    runner = adaptive.Runner(learner, executor=MPIPoolExecutor(),shutdown_executor=True,\
    goal=lambda l: l.loss() < loss
    runner.start_periodic_saving(dict(fname=fname), interval=600)
    runner.ioloop.run_until_complete(runner.task)
    learner.save(fname)

    dd=np.array(list(learner.data.items()))
    dz=dd[:,1]
    dx=np.empty(dd.shape[0])
    dy=np.empty(dd.shape[0])
    for i in range(dd.shape[0]):
        dx[i],dy[i]=dd[i,0]
    dz=np.vstack(dz)
    dxx, dyy = np.meshgrid(np.linspace(0,2,400),np.linspace(-.3,.3,400))
    dzz = griddata((dx,dy),dz[:,0],(dxx,dyy), method='linear')
    fig,ax=plt.subplots()
    ax.pcolormesh(dxx,dyy,dzz)
    fig.savefig('loss'+str(loss)+'.png')
    np.savetxt('loss'+str(loss)+'x.dat',dxx)
    np.savetxt('loss'+str(loss)+'y.dat',dyy)
    np.savetxt('loss'+str(loss)+'z.dat',dzz)

if __name__=="__main__":
        main()
