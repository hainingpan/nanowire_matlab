%%ham for self energy and QD
function [re,ham]=hseqd(a,mu,delta,vz,alpha_R,gamma,vc,omega,n,beta,mumax,l0,dim,period)
% t=25/a^2;
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alpha=alpha_R/(2*a);
band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim))+period*sparse([1,dim],[dim,1],[1,1],dim,dim);
band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim))+period*sparse([1,dim],[dim,1],[1,-1],dim,dim);
eyesm=speye(dim);
x=(1:l0)'-1;
mulist=[mu-mumax*exp(-x.*x/(l0*l0));mu*ones(dim-l0,1)];
%uncomment for the second QD on the right
% xR=(1:15)'-1;
% lR=15;
% mulist(end:-1:end-lR+1)=[mu-2.3*exp(-xR.*xR/(lR*lR))];
diagmulist=spdiags(mulist,0,dim,dim);
delta=delta*(sqrt(1-(vz/vc)^2))*(vz<vc);
deltalist=[zeros(l0,1);delta*ones(dim-l0,1)];
% deltalist(end:-1:end-lR+1)=0;
% diagdelta=spdiags(deltalist,0,dim,dim);
% ham=kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1i*alpha*band1m1sm)))+kron(eye(2),kron(sz,vz*eyesm))...
%     -gamma*real(kron(speye(4),omega./sqrt(diagdelta.^2-omega^2-1e-9i))+kron(sx,kron(speye(2),diagdelta./sqrt(diagdelta.^2-omega^2-1e-9i))));
ham=kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1i*alpha*band1m1sm)))+kron(eye(2),kron(sz,vz*eyesm))...
-gamma*real(kron(speye(4),spdiags((deltalist~=0).*omega./sqrt(deltalist.^2-omega^2-sign(omega+1e-9)*1e-9i),0,dim,dim))+kron(sx,kron(speye(2),spdiags((deltalist~=0).*deltalist./sqrt(deltalist.^2-omega^2-sign(omega+1e-9)*1e-9i),0,dim,dim))));

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

