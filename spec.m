function re=spec(mu,delta,alpha)
t=25;
vzm=1;
% vzm=((-alpha^2+2*t*mu)*sqrt(delta^2+mu^2)+sqrt(4*t*mu*alpha^2*(mu^2-delta^2)+(alpha^4+4*t^2*mu^2)*(delta^2+mu^2)))/(4*t*mu);
vzgrid=51;
vzset=linspace(0,vzm,vzgrid);
% vzstep=vzset(2)-vzset(1);
% vzset2=0:vzstep/100:vzm;
en=zeros(1,length(vzset));
parfor i=1:length(vzset)
     warning('off','all');
    vz=vzset(i);
    ham=hmb(1,mu,delta,vz,alpha,3000);    
    try 
       eigo=eigs(ham,10,'SM');     
        if (prod(isnan(eigo))==1)
            error("1 is not enough");
        end
   catch            
       try 
           eigo=eigs(ham,20,'SM');
            if (prod(isnan(eigo))==1)
                 error("20 is not enough");
            end
       catch
           eigo=eigs(ham,40,'SM');
            if (prod(isnan(eigo))==1)
                error("40 is not enough");
            end
       end
    end  
    en(i)=min(abs(eigo));
end
re=en;
end
