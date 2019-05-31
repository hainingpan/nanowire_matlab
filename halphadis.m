%% disorder in alpha
function re=halphadis(a,mu,delta,vz, dim,randlist)
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alphalist=randlist/(2*a);
band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim));
alphalist(1)=[];
% band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim));
band1m1sm=(spdiags([[alphalist;0],[0;-alphalist]],[-1,1],dim,dim));

eyesm=speye(dim);
re=kron(sz,(kron(eye(2),-t*band11sm+(2*t-mu)*eyesm)+kron(sy,1i*band1m1sm)))...
    +kron(eye(2),kron(sz,vz*eyesm))...
    +kron(sx,kron(eye(2),delta*eyesm));
end

