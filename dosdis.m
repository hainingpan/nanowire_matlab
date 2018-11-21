%% dos vs x
function re=dosdis(mu,Delta,Vz,alpha_R,dim,vimp,omega,delta)
a=1;
ham=hdis(a,mu,Delta,Vz,alpha_R,dim,vimp);
% G=inv(full((omega+1i*delta)*speye(4*dim)-ham)); 
% re=-imag(trace(G))/pi;
K=(omega+1i*delta)*speye(4*dim)-ham;
Kv=eigs(K,40,0,'Tolerance',1e-5,'MaxIterations',20000);
re=-sum(imag(1./Kv))/pi;
end