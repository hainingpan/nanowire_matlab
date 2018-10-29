function re=hmu(a,mu,delta,vz,alpha_R,dim,smoothpot,mumax,peakpos)
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alpha=alpha_R/(2*a);
band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim));
band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim));
eyesm=speye(dim);
switch smoothpot
    case 'sin2'
        mulist=sin((0:dim-1)*2*pi/dim)*mumax+mu;
    case 'cos'
        mulist=cos((0:dim-1)*pi/dim)*mumax+mu;
    case 'sigmoid'
        mulist=1/(exp(-((0:dim-1)-0.5*dim))+1));
    case 'lorentz'
        mulist=mumax*1.0./((((0:dim-1)-peakpos*dim)).^2+.5)+mu;
    case 'lorentzsigmoid'
        mulist=(mumax*1.0./(((0:dim-1)-peakpos*dim).^2+.5)+(4-mu)/2./(exp(-((0:dim-1)-0.5*dim))+1))+mu;
end      
re=kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diag(mulist))+kron(sy,1i*alpha*band1m1sm)))+kron(eye(2),kron(sz,vz*eyesm))+kron(sx,kron(eye(2),delta*eyesm));
end

