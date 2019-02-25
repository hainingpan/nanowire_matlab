function re=dossedis(a,mu,Delta,vz,alpha_R,gamma,vc,dim,vimp,omega,delta)
a=1;
[~,ham]=hsedis(a,mu,Delta,vz,alpha_R,gamma,vc,omega,1,0,dim,vimp);
K=(omega+1i*delta)*speye(4*dim)-ham;
Kv=eigs(K,40,0,'Tolerance',1e-5,'MaxIterations',20000);
re=-sum(imag(1./Kv))/pi;
end