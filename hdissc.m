%% disdorder & self energy & collapsing superconducting gap (deltalist=delta_0*sqrt(1-(VZ/Vc)^2))
function re=hdissc(a,mu,delta,vz,alpha_R,dim,vimp,vc)
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alpha=alpha_R/(2*a);
band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim));
band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim));
eyesm=speye(dim);
% mulist=mu*ones(dim,1);
mulist=mu*ones(dim)-vimp;
diagmulist=spdiags(mulist,0,dim,dim);
delta=delta*sqrt(1-(vz/vc)^2)*(vz<vc);
re=kron(sz,(kron(eye(2),-t*band11sm+(2*t)*eyesm-diagmulist)+kron(sy,1i*alpha*band1m1sm)))+kron(eye(2),kron(sz,vz*eyesm))+kron(sx,kron(eye(2),delta*eyesm));
end

