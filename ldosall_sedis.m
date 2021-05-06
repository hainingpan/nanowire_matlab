%% ldos vs x
function re=ldosall_sedis(a,mu,Delta,vz,alpha_R,gamma,vc,dim,mulist,omega,delta)
% a=1;
[~,ham]=hsedis(a,mu,Delta,vz,alpha_R,gamma,vc,omega,1,0,dim,mulist,0);
G=inv(full((omega+1i*delta)*speye(4*dim)-ham)); 
Gdiag=diag(G);
re=-imag(transpose(sum(reshape(Gdiag,[],4),2)))/pi;
end