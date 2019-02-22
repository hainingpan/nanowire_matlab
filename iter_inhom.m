function re=iter_inhom(a,mu,delta,vz,alpha_R,gamma,vc,n,dim,smoothpot,mumax,peakpos,sigma,omega)
% omega=0;
omega_n=hsemu(a,mu,delta,vz,alpha_R,gamma,vc,omega,n,0,dim,smoothpot,mumax,peakpos,sigma);
% k=(omega_n-0)/(0-delta);
% beta=-k;
beta=4;
count=0;
while abs(omega_n-omega)>1e-5 && count<3000
    omega=omega_n;
    count=count+1;
    omega_n=hsemu(a,mu,delta,vz,alpha_R,gamma,vc,omega,n,beta,dim,smoothpot,mumax,peakpos,sigma);
%     fprintf("omega=%f\n",omega_n);
end
re=omega_n;
% fprintf("count=%d",count);
end
