import matplotlib
matplotlib.use('Agg')
import numpy as np
# import tinyarray
import matplotlib.pyplot as plt
from scipy.sparse import spdiags
from scipy.sparse import eye
from scipy.sparse import kron
from scipy.sparse.linalg import inv
from scipy.sparse import csc_matrix
from functools import partial
from scipy.interpolate import griddata
import sys
from mpi4py import MPI
import argparse


s0 = np.array([[1, 0], [0, 1]]);
sx = np.array([[0, 1], [1, 0]]);
sy = np.array([[0, -1j], [1j, 0]]);
sz = np.array([[1, 0], [0, -1]]);

comm = MPI.COMM_WORLD;
rank = comm.Get_rank();
size=comm.Get_size();

def hdis(param):
    a,mu,delta0,alpha_R,muVarList,dim,vz=param['a'],param['mu'],param['delta0'],param['alpha_R'],param['muVarList'],param['wireLength'],param['vz']
    t=25/a**2
    alpha=alpha_R/(2*a)
    band11sm=spdiags(np.vstack([np.ones(dim),np.ones(dim)]),np.array([-1,1]),dim,dim,format = 'csc')
    band1m1sm=spdiags(np.vstack([np.ones(dim),-np.ones(dim)]),np.array([-1,1]),dim,dim,format = 'csc')
    eyesm=eye(dim)
    mulist=mu*np.ones(dim)-muVarList
    diagmulist=spdiags(mulist,0,dim,dim)
    return kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1j*alpha*band1m1sm)))\
    +kron(eye(2),kron(sz,vz*eyesm))+kron(sx,kron(eye(2),delta0*eyesm))

def hsedis(param):
    a,mu,delta0,alpha_R,couplingSCSM,muVarList,dim,vz,vc,vBias=param['a'],param['mu'],param['delta0'],param['alpha_R'],param['couplingSCSM'],param['muVarList'],param['wireLength'],param['vz'],param['vc'],param['vBias']
    t=25/a**2
    alpha=alpha_R/(2*a)
    band11sm=spdiags(np.vstack([np.ones(dim),np.ones(dim)]),np.array([-1,1]),dim,dim,format = 'csc')
    band1m1sm=spdiags(np.vstack([np.ones(dim),-np.ones(dim)]),np.array([-1,1]),dim,dim,format = 'csc')
    eyesm=eye(dim)
    eye4sm=eye(4*dim)
    delta0=delta0*(np.sqrt(1-(vz/vc)**2))*(vz<vc)
    mulist=mu*np.ones(dim)-muVarList
    diagmulist=spdiags(mulist,0,dim,dim)
    return kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1j*alpha*band1m1sm)))\
    +kron(eye(2),kron(sz,vz*eyesm))-couplingSCSM*(vBias/np.sqrt(delta0**2-vBias**2-np.sign(vBias+1e-9)*1e-9j)*eye4sm\
    +kron(sx,kron(eye(2),delta0/np.sqrt(delta0**2-vBias**2-np.sign(vBias+1e-9)*1e-9j)*eyesm)))

def ldosall_dis(param,eta):
    if param['isSE']==0:
        ham=hdis(param)
    else:
        ham=hsedis(param)
    vBias,dim=param['vBias'],param['wireLength']
    hh=csc_matrix((vBias+1j*eta)*eye(4*dim)-ham)
    G=inv(hh)
    Gdiag=(G).diagonal()
    return -np.sum((np.reshape(Gdiag,(4,-1))),0).imag/np.pi

def LDOS_dis(param,eta):
    vz,vBias,dim=param['vz'],param['vBias'],param['wireLength']
    z=ldosall_dis(param,eta)
    return np.array([z.mean(),z[0],z[int(dim/2)],z[-1]])

def main():
    parameters={'error':0}
    if rank==0:
        parser = argparse.ArgumentParser()
        parser.add_argument('--a', default=1)
        parser.add_argument('--wireLength', default=100)
        parser.add_argument('--mu', default=1)
        parser.add_argument('--delta0', default=0.2)
        parser.add_argument('--alpha_R', default=5)
        parser.add_argument('--muVar', default=0)
        parser.add_argument('--muVarList', default=0)
        parser.add_argument('--muVarNum', default=0)
        parser.add_argument('--vzMax', default=2.048)
        parser.add_argument('--vBiasmax', default=0.3)
        parser.add_argument('--vzNum', default=200)
        parser.add_argument('--vBiasNum', default=200)
        parser.add_argument('--vz', default=0)
        parser.add_argument('--vc', default=float('inf'))
        parser.add_argument('--isSE', default=0)
        parser.add_argument('--couplingSCSM', default=0.2)

        args = parser.parse_args();

        parameters['a']=float(args.a)
        parameters['wireLength']=int(args.wireLength)
        parameters['mu']=float(args.mu)
        parameters['delta0']=float(args.delta0)
        parameters['alpha_R']=float(args.alpha_R)
        parameters['muVar']=float(args.muVar)
        parameters['muVarNum']=int(args.muVarNum)*(args.muVarNum>0)+int(args.wireLength)*(args.muVarNum==0)
        parameters['vzMax']=float(args.vzMax)
        parameters['vBiasmax']=float(args.vBiasmax)
        parameters['vzNum']=int(args.vzNum)
        parameters['vBiasNum']=int(args.vBiasNum)
        parameters['vz']=float(args.vz)
        parameters['vc']=float(args.vc)
        parameters['isSE']=int(args.isSE)
        parameters['couplingSCSM']=float(args.couplingSCSM)


        if isinstance(args.muVarList,str):
            muVarfn=args.muVarList
            print('Use disorder file:',muVarfn)
            try:
                parameters['muVarList']=np.loadtxt(muVarfn)
            except:
                print('Cannot find disorder file: ',muVarfn)
                parameters['error']=1
                sys.exit(1)
        elif parameters['muVar']!=0:
            muVarList=np.random.normal(0,parameters['muVar'],int(parameters['muVarNum']))
            parameters['muVarList']=[muVarList.flatten()[int(parameters['muVarNum']/parameters['wireLength']*x)] for x in range(parameters['wireLength'])]
        else:
            parameters['muVarList']=args.muVarList
        print(parameters)

    parameters=comm.bcast(parameters,root=0)
    if parameters['error']!=0:   #for the slave to exit
        sys.exit(1)

    vzNum,vBiasNum=parameters['vzNum'],parameters['vBiasNum']
    per=int(vzNum/size)
    sendbuf_DOS=np.empty((per,vBiasNum))
    sendbuf_LDOS_L=np.empty((per,vBiasNum))
    sendbuf_LDOS_M=np.empty((per,vBiasNum))
    sendbuf_LDOS_R=np.empty((per,vBiasNum))
    vBiasList=np.linspace(-parameters['vBiasmax'],parameters['vBiasmax'],vBiasNum)
    vzstep=parameters['vzMax']/(vzNum-1)
    for i in range(per):
        for j in range(vBiasNum):
            parameters['vz']=(i+per*rank)*vzstep
            parameters['vBias']=vBiasList[j]
            sys.stdout.flush()
            sendbuf_DOS[i,j],sendbuf_LDOS_L[i,j],sendbuf_LDOS_M[i,j],sendbuf_LDOS_R[i,j]=LDOS_dis(param=parameters,eta=1e-3)
    if rank==0:
        recvbuf_DOS=np.empty((vzNum,vBiasNum))
        recvbuf_LDOS_L=np.empty((vzNum,vBiasNum))
        recvbuf_LDOS_M=np.empty((vzNum,vBiasNum))
        recvbuf_LDOS_R=np.empty((vzNum,vBiasNum))
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
        plt.rcParams['pcolor.shading']='auto'
        fn_mu=('m'+str(parameters['mu']))
        fn_Delta='D'+str(parameters['delta0'])
        fn_alpha='a'+str(parameters['alpha_R'])
        fn_muVar=('mVar'+str(parameters['muVar']))*(parameters['muVar']!=0)
        fn_wl='L'+str(int(parameters['wireLength']))
        fn_couplingSCSM=('g'+str(parameters['couplingSCSM']))*(parameters['isSE']==1)
        fn_vc=('vc'+str(parameters['vc']))*(parameters['isSE']==1)*(parameters['vc']!=float('inf'))
        fn_range=('-'+str(parameters['vzMax'])+','+str(parameters['vBiasmax'])+'-')
        fn=fn_mu+fn_Delta+fn_muVar+fn_wl+fn_couplingSCSM+fn_vc+fn_range
        fig,ax=plt.subplots(2,2,tight_layout=True)
        ax[0,0].pcolormesh(np.linspace(0,parameters['vzMax'],vzNum),np.linspace(-parameters['vBiasmax'],parameters['vBiasmax'],vBiasNum),recvbuf_DOS.T)
        # ax[0,0].set_xlabel('Vz(meV)')
        ax[0,0].set_ylabel(r'$V_\mathrm{bias}$ (meV)')
        ax[0,0].set_title('DOS')
        # fig.savefig(fn+'_DOS.png')

        # fig,ax=plt.subplots()
        ax[0,1].pcolormesh(np.linspace(0,parameters['vzMax'],vzNum),np.linspace(-parameters['vBiasmax'],parameters['vBiasmax'],vBiasNum),recvbuf_LDOS_L.T)
        # ax[0,1].set_xlabel('Vz(meV)')
        # ax[0,1].set_ylabel(r'$V_\mathrm{bias}$ (meV)')
        ax[0,1].set_title('LDOS(Left)')
        # fig.savefig(fn+'_LDOS_L.png')

        # fig,ax=plt.subplots()
        ax[1,0].pcolormesh(np.linspace(0,parameters['vzMax'],vzNum),np.linspace(-parameters['vBiasmax'],parameters['vBiasmax'],vBiasNum),recvbuf_LDOS_M.T)
        ax[1,0].set_xlabel('Vz(meV)')
        ax[1,0].set_ylabel(r'$V_\mathrm{bias}$ (meV)')
        ax[1,0].set_title('LDOS(Middle)')
        # fig.savefig(fn+'_LDOS_M.png')

        # fig,ax=plt.subplots()
        ax[1,1].pcolormesh(np.linspace(0,parameters['vzMax'],vzNum),np.linspace(-parameters['vBiasmax'],parameters['vBiasmax'],vBiasNum),recvbuf_LDOS_R.T)
        ax[1,1].set_xlabel('Vz(meV)')
        # ax[1,1].set_ylabel(r'$V_\mathrm{bias}$ (meV)')
        ax[1,1].set_title('LDOS(Right)')
        # fig.savefig(fn+'_LDOS_R.png')

        fig.savefig(fn+'LDOS.png')

        np.savetxt(fn+'_DOS.dat',recvbuf_DOS)
        np.savetxt(fn+'_LDOS_L.dat',recvbuf_LDOS_L)
        np.savetxt(fn+'_LDOS_M.dat',recvbuf_LDOS_M)
        np.savetxt(fn+'_LDOS_R.dat',recvbuf_LDOS_R)


        if parameters['muVar']!=0:
            np.savetxt(fn+'_randlist.dat',parameters['muVarList'])

if __name__=="__main__":
        main()
