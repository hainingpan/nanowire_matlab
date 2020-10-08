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
import argparse


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

    parser = argparse.ArgumentParser()
    parser.add_argument('--loss', default=0.1);
    parser.add_argument('--dim', default=100);
    parser.add_argument('--mu', default=1);
    parser.add_argument('--Delta', default=0.2);
    parser.add_argument('--alpha_R', default=5);
    parser.add_argument('--muVar', default=0);
    args = parser.parse_args();

    print("loss = %s" % args.loss)
    print("dim = %s" % args.dim)
    print("mu = %s" % args.mu)
    print("Delta = %s" % args.Delta)
    print("alpha_R = %s" % args.alpha_R)
    print("muVar = %s" % args.muVar)

    loss=float(args.loss)
    dim=int(args.dim)
    mu=float(args.mu)
    Delta=float(args.Delta)
    alpha_R=float(args.alpha_R)
    muVar=float(args.muVar)

    fn='loss'+str(loss)+'m'+str(mu)+'D'+str(Delta)+'muVar'+str(muVar)+'L'+str(dim)
    fname=fn+'.sav'
    learner = adaptive.Learner2D(partial(LDOS_dis,a=1,mu=1,Delta=0.2,alpha_R=5,mulist=0,dim=1000,delta=1e-3),\
                                 bounds=[(0., 2.), (-0.3, 0.3)])
    learner.load(fname)
    runner = adaptive.Runner(learner, executor=MPIPoolExecutor(),shutdown_executor=True,\
        goal=lambda l: l.loss() < loss)
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
    dxx, dyy = np.meshgrid(np.linspace(0,2.048,401),np.linspace(-.3,.3,401))
    dzz0 = griddata((dx,dy),dz[:,0],(dxx,dyy), method='linear')
    dzz1 = griddata((dx,dy),dz[:,1],(dxx,dyy), method='linear')
    dzz2 = griddata((dx,dy),dz[:,2],(dxx,dyy), method='linear')
    dzz3 = griddata((dx,dy),dz[:,3],(dxx,dyy), method='linear')
    fig,ax=plt.subplots()
    ax.pcolormesh(dxx,dyy,dzz0)
    fig.savefig(fn+'_DOS.png')

    fig,ax=plt.subplots()
    ax.pcolormesh(dxx,dyy,dzz1)
    fig.savefig(fn+'_LDOS_L.png')

    fig,ax=plt.subplots()
    ax.pcolormesh(dxx,dyy,dzz2)
    fig.savefig(fn+'_LDOS_M.png')

    fig,ax=plt.subplots()
    ax.pcolormesh(dxx,dyy,dzz3)
    fig.savefig(fn+'_LDOS_R.png')

    np.savetxt(fn+'_Vz.dat',dxx)
    np.savetxt(fn+'_Vbias.dat',dyy)
    np.savetxt(fn+'_DOS.dat',dzz0)
    np.savetxt(fn+'_LDOS_L.dat',dzz1)
    np.savetxt(fn+'_LDOS_M.dat',dzz2)
    np.savetxt(fn+'_LDOS_R.dat',dzz3)
    
    scatterpts=np.vstack([dx,dy,dz.T]).T
    np.savetxt(fn+'_s.dat',scatterpts)

if __name__=="__main__":
        main()
