function re=WF(a,mu,delta,vz,alpha_R,dim,index)
nv=2*index;
ham=hc(a,mu,delta,vz,alpha_R,dim);
[V,D]=eigs(ham,nv,0,'Tolerance',1e-5,'MaxIterations',20000);
wfpos=V(:,1);
wfneg=V(:,2);
wfpos2=reshape(wfpos,dim,[]);
wfneg2=reshape(wfneg,dim,[]);
wfpos3=sum(abs(wfpos2).^2,2);
wfneg3=sum(abs(wfneg2).^2,2);
re=(wfpos3+wfneg3)/2;
end