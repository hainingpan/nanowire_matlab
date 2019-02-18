function re=iter_seqd(a,mu,delta,vz,alpha_R,gamma,vc,mumax,l0,n,dim)
omega=0;
omega_n=hseqd(a,mu,delta,vz,alpha_R,gamma,vc,omega,n,0,mumax,l0,dim);
% k=(omega_n-0)/(0-delta);
% beta=-k;
beta=4;
count=0;
while abs(omega_n-omega)>1e-5 && count<3000
    omega=omega_n;
    count=count+1;
    omega_n=hseqd(a,mu,delta,vz,alpha_R,gamma,vc,omega,n,beta,mumax,l0,dim);
%     fprintf("omega=%f\n",omega_n);
end
re=omega_n;
% fprintf("count=%d",count);
end
