%% ldos vs x
function re=ldosall_dis(a,mu,Delta,Vz,alpha_R,mulist,dim,omega,delta)
% ham=hmu(a,mu,Delta,Vz,alpha_R,dim,smoothpot,mumax,peakpos,sigma);
ham=hdis(a,mu,Delta,Vz,alpha_R,dim,mulist,0);
G=inv(full((omega+1i*delta)*speye(4*dim)-ham)); 
Gdiag=diag(G);
re=-imag(transpose(sum(reshape(Gdiag,[],4),2)))/pi;
end