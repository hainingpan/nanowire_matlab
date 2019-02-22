function re=iter_dis(a,mu,delta,vz,alpha_R,gamma,vc,n,dim,vimp,omega)
% omega=0;
omega_n=hsedis(a,mu,delta,vz,alpha_R,gamma,vc,omega,n,0,dim,vimp);
% k=(omega_n-0)/(0-delta);
% beta=-k;
beta=4;
count=0;
while abs(omega_n-omega)>1e-5 && count<3000
    omega=omega_n;
    count=count+1;
    omega_n=hsedis(a,mu,delta,vz,alpha_R,gamma,vc,omega,n,beta,dim,vimp);
%     fprintf("omega=%f\n",omega_n);
end
re=omega_n;
% fprintf("count=%d",count);
end
