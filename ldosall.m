function re=ldosall(mu,Delta,Vz,alpha_R,dim,smoothpot,mumax,peakpos,sigma,omega,delta)
a=1;
ham=hmu(a,mu,Delta,Vz,alpha_R,dim,smoothpot,mumax,peakpos,sigma);
G=inv(full((omega+1i*delta)*speye(4*dim)-ham)); %check when full should be changed
Gdiag=diag(G);
re=-imag(transpose(sum(reshape(Gdiag,[],4),2)))/pi;
end