function ham=hse_h(a,mu,delta,vz,alpha_R,gamma,omega,dim)
t=25/a^2;
sx=[0,1;1,0];
sy=[0,-1i;1i,0];
sz=[1,0;0,-1];
alpha=alpha_R/(2*a);
band11sm=(spdiags([ones(dim,1) ones(dim,1)],[-1,1],dim,dim));
band1m1sm=(spdiags([ones(dim,1) -ones(dim,1)],[-1,1],dim,dim));
eyesm=speye(dim);
eye4sm=speye(4*dim);
ham=kron(sz,(kron(eye(2),-t*band11sm+(-mu+2*t)*eyesm)+kron(sy,1i*alpha*band1m1sm)))+kron(eye(2),kron(sz,vz*eyesm))-gamma*(omega/sqrt(abs(delta^2-omega^2))*eye4sm+kron(sx,kron(eye(2),delta/sqrt(abs(delta^2-omega^2))*eyesm)));
ham=ham-omega*speye(4*dim);
end