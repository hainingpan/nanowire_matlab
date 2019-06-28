%% random g factor
function re=hg(a,mu,delta,randlist,alpha_R,dim)
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alpha=alpha_R/(2*a);
band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim));
band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim));
eyesm=speye(dim);
% mulist=mu*ones(dim,1);
% mulist=mu-randlist;
diagvzlist=spdiags(randlist,0,dim,dim);
re=kron(sz,(kron(eye(2),-t*band11sm+(2*t-mu)*eyesm)+kron(sy,1i*alpha*band1m1sm)))+kron(speye(2),kron(sz,diagvzlist))+kron(sx,kron(eye(2),delta*eyesm));
end

