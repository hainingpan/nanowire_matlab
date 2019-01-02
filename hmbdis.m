function re=hmbdis(a,mu,delta11,vz,alpha_R,deltac,dim,vimp)
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alpha=alpha_R/(2*a);
band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim));
band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim));
eyesm=speye(dim);
epsilon=1;
mulist=mu-vimp;
diagmulist=spdiags(mulist,0,dim,dim);
H11=kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1i*alpha*band1m1sm)))+kron(eye(2),kron(sz,vz*eyesm))+kron(sx,kron(eye(2),delta11*eyesm));
H22=kron(sz,(kron(eye(2),-t*band11sm+(epsilon+2*t)*eyesm-diagmulist)+kron(sy,1i*alpha*band1m1sm)))+kron(eye(2),kron(sz,vz*eyesm))+kron(sx,kron(eye(2),delta11*eyesm));
H12=kron(sx,kron(eye(2),deltac*eyesm));
re=[H11,H12;H12,H22];
end

