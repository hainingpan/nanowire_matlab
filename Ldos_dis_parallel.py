import matplotlib
matplotlib.use('Agg')
import numpy as np
import tinyarray
import matplotlib.pyplot as plt
from scipy.sparse import spdiags
from scipy.sparse import eye
from scipy.sparse import kron
from scipy.sparse.linalg import inv
from scipy.sparse import csc_matrix
import adaptive
from functools import partial
from scipy.interpolate import griddata
import sys
from mpi4py import MPI
import argparse


s0 = tinyarray.array([[1, 0], [0, 1]]);
sx = tinyarray.array([[0, 1], [1, 0]]);
sy = tinyarray.array([[0, -1j], [1j, 0]]);
sz = tinyarray.array([[1, 0], [0, -1]]);

comm = MPI.COMM_WORLD;
rank = comm.Get_rank();
size=comm.Get_size();
# print('size:'+str(size)+'rank:'+str(rank))
def hdis(a,mu,delta,vz,alpha_R,dim,vimp):
    t=25/a**2
    alpha=alpha_R/(2*a)
    band11sm=spdiags(np.vstack([np.ones(dim),np.ones(dim)]),np.array([-1,1]),dim,dim,format = 'csc')
    band1m1sm=spdiags(np.vstack([np.ones(dim),-np.ones(dim)]),np.array([-1,1]),dim,dim,format = 'csc')
    eyesm=eye(dim)
    mulist=mu*np.ones(dim)-vimp
    diagmulist=spdiags(mulist,0,dim,dim)
    return kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1j*alpha*band1m1sm)))\
    +kron(eye(2),kron(sz,vz*eyesm))+kron(sx,kron(eye(2),delta*eyesm))

def ldosall_dis(a,mu,Delta,Vz,alpha_R,mulist,dim,omega,delta):
    ham=hdis(a,mu,Delta,Vz,alpha_R,dim,mulist);
    hh=csc_matrix((omega+1j*delta)*eye(4*dim)-ham)
    G=inv(hh)
    Gdiag=(G).diagonal()
    return -np.sum((np.reshape(Gdiag,(4,-1))),0).imag/np.pi

def LDOS_dis(p,a,mu,Delta,alpha_R,mulist,dim,delta):
    Vz,energy=p
    z=ldosall_dis(a,mu,Delta,Vz,alpha_R,mulist,dim,energy,delta)
    return np.array([z.mean(),z[0],z[int(dim/2)],z[-1]])

def main():
    parameters={'error':0}
    if rank==0:
        parser = argparse.ArgumentParser()
        parser.add_argument('--dim', default=100)
        parser.add_argument('--mu', default=1)
        parser.add_argument('--Delta', default=0.2)
        parser.add_argument('--alpha_R', default=5)
        parser.add_argument('--muVar', default=0)
        parser.add_argument('--mulist', default=0)
        parser.add_argument('--NmuVar', default=0)
        parser.add_argument('--Vzmax', default=2.048)
        parser.add_argument('--Vbiasmax', default=0.3)
        parser.add_argument('--tot_vz', default=200)
        parser.add_argument('--tot_energy', default=200)
        args = parser.parse_args();

        print("dim = %s" % args.dim)
        print("mu = %s" % args.mu)
        print("Delta = %s" % args.Delta)
        print("alpha_R = %s" % args.alpha_R)
        print("muVar = %s" % args.muVar)
        print("mulist = %s" % args.mulist)
        print("NmuVar = %s" % args.NmuVar)
        print("Vzmax = %s" % args.Vzmax)
        print("Vbiasmax = %s" % args.Vbiasmax)

        parameters['dim']=int(args.dim)
        parameters['mu']=float(args.mu)
        parameters['Delta']=float(args.Delta)
        parameters['alpha_R']=float(args.alpha_R)
        parameters['muVar']=float(args.muVar)
        parameters['NmuVar']=float(args.NmuVar)
        parameters['Vzmax']=float(args.Vzmax)
        parameters['Vbiasmax']=float(args.Vbiasmax)
        parameters['tot_vz']=int(args.tot_vz)
        parameters['tot_energy']=int(args.tot_energy)
        if isinstance(args.mulist,str):
            muVarfn=args.mulist
            print('Use disorder file:',muVarfn)
            try:
                parameters['mulist']=np.loadtxt(muVarfn)
            except:
                print('Cannot find disorder file: ',muVarfn)
                parameters['error']=1
                sys.exit(1)
        elif parameters['muVar']!=0:
            mulist=np.random.normal(0,parameters['muVar'],int(parameters['NmuVar']))
            parameters['mulist']=[mulist.flatten()[int(parameters['NmuVar']/parameters['dim']*x)] for x in range(parameters['dim'])]
        else:
            parameters['mulist']=args.mulist


    parameters=comm.bcast(parameters,root=0)
    if parameters['error']!=0:   #for the slave to exit
        sys.exit(1)

    func=partial(LDOS_dis,a=1,mu=parameters['mu'],Delta=parameters['Delta'],alpha_R=parameters['alpha_R'],\
        mulist=parameters['mulist'],dim=parameters['dim'],delta=1e-3)
    tot_vz,tot_energy=parameters['tot_vz'],parameters['tot_energy']
    per=int(tot_vz/size)
    sendbuf_DOS=np.empty((per,tot_energy))
    sendbuf_LDOS_L=np.empty((per,tot_energy))
    sendbuf_LDOS_M=np.empty((per,tot_energy))
    sendbuf_LDOS_R=np.empty((per,tot_energy))
    energylist=np.linspace(-parameters['Vbiasmax'],parameters['Vbiasmax'],tot_energy)
    vzstep=parameters['Vzmax']/(tot_vz-1)
    for i in range(per):
        for j in range(tot_energy):
            Vz=(i+per*rank)*vzstep
            energy=energylist[j]
            # print('Vz: '+str(Vz)+' energy: '+str(energy)+' by rank: '+str(rank))
            # sys.stdout.flush()
            sendbuf_DOS[i,j],sendbuf_LDOS_L[i,j],sendbuf_LDOS_M[i,j],sendbuf_LDOS_R[i,j]=func((Vz,energy))


    if rank==0:
        recvbuf_DOS=np.empty((tot_vz,tot_energy))
        recvbuf_LDOS_L=np.empty((tot_vz,tot_energy))
        recvbuf_LDOS_M=np.empty((tot_vz,tot_energy))
        recvbuf_LDOS_R=np.empty((tot_vz,tot_energy))
    else:
        recvbuf_DOS=None
        recvbuf_LDOS_L=None
        recvbuf_LDOS_M=None
        recvbuf_LDOS_R=None
    comm.Gather(sendbuf_DOS,recvbuf_DOS,root=0)
    comm.Gather(sendbuf_LDOS_L,recvbuf_LDOS_L,root=0)
    comm.Gather(sendbuf_LDOS_M,recvbuf_LDOS_M,root=0)
    comm.Gather(sendbuf_LDOS_R,recvbuf_LDOS_R,root=0)

    if rank==0:
        fn='m'+str(parameters['mu'])+'D'+str(parameters['Delta'])+('muVar'+str(parameters['muVar']))*(parameters['muVar']!=0)+'L'+str(parameters['dim'])
        fig,ax=plt.subplots()
        ax.pcolormesh(np.linspace(0,parameters['Vzmax'],tot_vz),np.linspace(0,parameters['Vbiasmax'],tot_energy),recvbuf_DOS.T)
        fig.savefig(fn+'_DOS.png')

        fig,ax=plt.subplots()
        ax.pcolormesh(np.linspace(0,parameters['Vzmax'],tot_vz),np.linspace(0,parameters['Vbiasmax'],tot_energy),recvbuf_LDOS_L.T)
        fig.savefig(fn+'_LDOS_L.png')

        fig,ax=plt.subplots()
        ax.pcolormesh(np.linspace(0,parameters['Vzmax'],tot_vz),np.linspace(0,parameters['Vbiasmax'],tot_energy),recvbuf_LDOS_M.T)
        fig.savefig(fn+'_LDOS_M.png')

        fig,ax=plt.subplots()
        ax.pcolormesh(np.linspace(0,parameters['Vzmax'],tot_vz),np.linspace(0,parameters['Vbiasmax'],tot_energy),recvbuf_LDOS_R.T)
        fig.savefig(fn+'_LDOS_R.png')

        np.savetxt(fn+'_DOS.dat',recvbuf_DOS)
        np.savetxt(fn+'_LDOS_L.dat',recvbuf_DOS)
        np.savetxt(fn+'_LDOS_M.dat',recvbuf_DOS)
        np.savetxt(fn+'_LDOS_R.dat',recvbuf_DOS)


        if parameters['muVar']!=0:
            np.savetxt(fn+'_randlist.dat',parameters['mulist'])

if __name__=="__main__":
        main()
