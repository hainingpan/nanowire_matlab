function [rev,re]=spec(mu,delta,alpha,deltac)
t=25;
vzm=1;
% vzm=((-alpha^2+2*t*mu)*sqrt(delta^2+mu^2)+sqrt(4*t*mu*alpha^2*(mu^2-delta^2)+(alpha^4+4*t^2*mu^2)*(delta^2+mu^2)))/(4*t*mu);
vzgrid=101;
vzset=linspace(0,vzm,vzgrid);
% vzstep=vzset(2)-vzset(1);
% vzset2=0:vzstep/100:vzm;
nv=30;
en=zeros(nv,length(vzset));
parfor i=1:length(vzset)
%      warning('off','all');
    vz=vzset(i);
    ham=hmb(1,mu,delta,vz,alpha,deltac,3000);    
%% For 1e-14 tolerance and single minimal value
%     try 
%        eigo=eigs(ham,1,'SM');     
%         if (prod(isnan(eigo))==1)
%             error("1 is not enough");
%         end
%    catch            
%        try 
%            eigo=eigs(ham,20,'SM');
%             if (prod(isnan(eigo))==1)
%                  error("20 is not enough");
%             end
%        catch
%            eigo=eigs(ham,40,'SM');
%             if (prod(isnan(eigo))==1)
%                 error("40 is not enough");
%             end
%        end
%     end  
%     en(i)=min(abs(eigo));

%% nv of smallest eigenvalue
    eigo=eigs(ham,200,0,'Tolerance',1e-5,'MaxIterations',10000);
    if (abs(eigo(1)) <1e-10 && abs(eigo(2))<1e-10 )
        if(eigo(1)*eigo(2)>0)
            eigo(2)=[];
            eigo(1)=abs(eigo(1));
        else
            if (eigo(1)<0)
                eigo(1)=[];
            else
                eigo(2)=[];
            end
        end
    end         
            
    eigo2=eigo(eigo>=0);
    en(:,i)=eigo2(1:60);
end
re=en;
rev=vzset;
end
