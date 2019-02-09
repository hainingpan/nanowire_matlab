function re=iter(a,mu,delta,vz,alpha_R,gamma,n,dim)
omega=0;
omega_n=hse(a,mu,delta,vz,alpha_R,gamma,omega,n,0,dim);
k=(omega_n-0)/(0-delta);
beta=-k;
count=0;
while abs(omega_n-omega)>1e-6 && count<300
    omega=omega_n;
    count=count+1;
    omega_n=hse(a,mu,delta,vz,alpha_R,gamma,omega,n,beta,dim);
end
re=omega_n;
% fprintf("k=%f count=%d",k,count);
end
