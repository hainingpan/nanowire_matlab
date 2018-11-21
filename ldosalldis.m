%% ldos vs x
function re=ldosalldis(mu,Delta,Vz,alpha_R,dim,pos,omega,delta)
a=1;
ham=hdis(a,mu,Delta,Vz,alpha_R,dim,pos);
G=inv(full((omega+1i*delta)*speye(4*dim)-ham)); 
Gdiag=diag(G);
re=-imag(transpose(sum(reshape(Gdiag,[],4),2)))/pi;
end