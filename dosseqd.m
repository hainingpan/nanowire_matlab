function re=dosseqd(mu,Delta,vz,alpha_R,gamma,vc,mumax,l0,dim,omega,delta)
a=1;
[~,ham]=hseqd(a,mu,Delta,vz,alpha_R,gamma,vc,omega,n,beta,mumax,l0,dim);
K=(omega+1i*delta)*speye(4*dim)-ham;
Kv=eigs(K,40,0,'Tolerance',1e-5,'MaxIterations',20000);
re=-sum(imag(1./Kv))/pi;
end