function re=ldosdis(mu,Delta,Vz,alpha_R,dim,vimp,x,omega,delta)
a=1;
ham=hdis(a,mu,Delta,Vz,alpha_R,dim,vimp);
G=inv(full((omega+1i*delta)*speye(4*dim)-ham)); %check when full should be changed
Gdiag=diag(G);
re=-imag(Gdiag(x)+Gdiag(x+dim)+Gdiag(x+2*dim)+Gdiag(x+3*dim))/pi;
end