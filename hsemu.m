%%smooth potential confinement
function [re,ham]=hsemu(a,mu,delta,vz,alpha_R,gamma,vc,omega,n,beta,dim,smoothpot,mumax,peakpos,sigma)
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alpha=alpha_R/(2*a);
band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim));
band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim));
eyesm=speye(dim);
eye4sm=speye(4*dim);
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
%         mumaxR=-1.9;
%         lR=100;
%         sigmaR=10;
%         mulist=mulist+(mumaxR*exp(-(site*a-lR).^2/(2*sigmaR^2)));
    case 'disorder'
        disorderpos=randperm(dim,peakpos);
        mulist=mu*ones(dim,1);
        mulist(disorderpos)=mumax;
end      
diagmulist=spdiags(mulist,0,dim,dim);
delta=delta*(sqrt(1-(vz/vc)^2))*(vz<vc);

ham=kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1i*alpha*band1m1sm)))+kron(eye(2),kron(sz,vz*eyesm))...
    -gamma*(omega/sqrt(delta^2-omega^2-sign(omega+1e-9)*1e-9i)*eye4sm+kron(sx,kron(speye(2),delta/sqrt(delta^2-omega^2-sign(omega+1e-9)*1e-9i)*eyesm)));
%  try 
%        eigo=eigs(ham,n,'SM','Tolerance',1e-6,'MaxIterations',10000);     
%         if (prod(isnan(eigo))==1)
%             disp("Here goes 1");
%             error("1 is not enough");
%         end
%    catch            
%        try 
%            eigo=eigs(ham,20,'SM','Tolerance',1e-6,'MaxIterations',10000);
%             if (prod(isnan(eigo))==1)
%                  disp("Here goes 20");
%                  error("20 is not enough");
%             end
%        catch
%            eigo=eigs(ham,40,'SM','Tolerance',1e-6,'MaxIterations',10000);
%             if (prod(isnan(eigo))==1)
%                 disp("Here goes 40");
%                 error("40 is not enough");
%             end
%        end
%  end
%  re=(abs(eigo(n))+beta*omega)/(beta+1); 
re=0;
end

