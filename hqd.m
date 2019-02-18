%% Quantum dots (gaussian)
function re=hqd(a,mu,delta,vz,alpha_R,mumax,l0,dim)
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alpha=alpha_R/(2*a);
band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim));
band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim));
deltalist=[zeros(l0,1);delta*ones(dim-l0,1)];
diagdelta=spdiags(deltalist,0,dim,dim);
eyesm=speye(dim);
% mulist=mu*ones(dim,1);
x=(1:l0)';
mulist=[mu-mumax*exp(-x.*x/(2*l0*l0));mu*ones(dim-l0,1)];
diagmulist=spdiags(mulist,0,dim,dim);
re=kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1i*alpha*band1m1sm)))+kron(eye(2),kron(sz,vz*eyesm))+kron(sx,kron(eye(2),diagdelta));
end

