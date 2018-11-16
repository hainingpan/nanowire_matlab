%%smooth potential confinement
function re=hmu(a,mu,delta,vz,alpha_R,dim,smoothpot,mumax,peakpos,sigma)
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alpha=alpha_R/(2*a);
band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim));
band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim));
eyesm=speye(dim);
site=(0:dim-1)';
switch smoothpot
    case 'const'
        mulist=ones(dim,1)*mu;
    case 'sin2'
        mulist=sin(site*2*pi/dim)*mumax+mu;
    case 'cos'
        mulist=cos(site*pi/dim)*mumax+mu;
    case 'sigmoid'
        mulist=mumax*1./(exp(-(site-0.5*dim)*a/sigma)+1)+mu;
    case 'lorentz'
        mulist=mumax*1.0./(((site-peakpos*dim)*a).^2+.5)+mu;
    case 'lorentzsigmoid'
        mulist=(mumax*1.0./(((site-peakpos*dim)*a).^2+.5)+(4-mu)/2./(exp(-((site-0.5*dim)*a))+1))+mu;
    case 'exp'
        mulist=(mumax*exp(-(site*a).^2/(2*sigma^2)))+mu;
    case 'disorder'
        disorderpos=randperm(dim,peakpos);
        mulist=mu*ones(dim,1);
        mulist(disorderpos)=mumax;
end      
diagmulist=spdiags(mulist,0,dim,dim);
re=kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1i*alpha*band1m1sm)))+kron(eye(2),kron(sz,vz*eyesm))+kron(sx,kron(eye(2),delta*eyesm));
end

