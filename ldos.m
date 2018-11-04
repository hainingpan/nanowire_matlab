function re=ldos(mu,Delta,Vz,alpha_R,dim,smoothpot,mumax,peakpos,sigma,x,omega,delta)
a=1;
ham=hmu(a,mu,Delta,Vz,alpha_R,dim,smoothpot,mumax,peakpos,sigma);
G=inv(full((omega+1i*delta)*speye(4*dim)-ham)); %check when full should be changed
Gdiag=diag(G);
re=-imag(Gdiag(x)+Gdiag(x+dim)+Gdiag(x+2*dim)+Gdiag(x+3*dim))/pi;
end