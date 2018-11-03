function re=ldos(Vz,x,omega,delta)
a=1;
mu=0.2;
Delta=0.2;
dim=200;
alpha_R=5;
ham=hc(a,mu,Delta,Vz,alpha_R,dim);
G=inv(omega*speye(4*dim)-ham+1i*delta);
Gdiag=diag(G);
re=-imag(Gdiag(x)+Gdiag(x+dim)+Gdiag(x+2*dim)+Gdiag(x+3*dim))/pi;
end