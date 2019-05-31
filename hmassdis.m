%% disorder in effective mass
function re=hmassdis(a,mu,delta,vz,alpha_R,dim,randlist)
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alpha=alpha_R/(2*a);
tlist=t./randlist;
% band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim));
band11sm=(spdiags([[tlist(1:end-1);0],[0;tlist(2:end)]],[-1,1],dim,dim));
band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim));
diagtlist=spdiags(2*tlist,0,dim,dim);

eyesm=speye(dim);
re=kron(sz,(kron(eye(2),-band11sm+(-mu)*eyesm+diagtlist)+kron(sy,1i*alpha*band1m1sm)))...
    +kron(eye(2),kron(sz,vz*eyesm))...
    +kron(sx,kron(eye(2),delta*eyesm));
end

