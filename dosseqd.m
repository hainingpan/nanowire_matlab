function re=dosseqd(a,mu,Delta,vz,alpha_R,gamma,vc,mumax,l0,dim,omega,delta,period)
% a=1;
[~,ham]=hseqd(a,mu,Delta,vz,alpha_R,gamma,vc,omega,1,0,mumax,l0,dim,period);
K=(omega+1i*delta)*speye(4*dim)-ham;
Kv=eigs(K,40,0,'Tolerance',1e-5,'MaxIterations',20000);
re=-sum(imag(1./Kv))/pi;
end