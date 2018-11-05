function [rev,re]=LDOS_Vz(mu,dim,smoothpot,mumax,peakpos,sigma,x)
Delta=0.2;
alpha=5;
delta=1e-3;
vzlist=linspace(0.,4,200);
omegalist=linspace(-.3,.3,200);
en=zeros(length(vzlist),length(omegalist));
parfor i=1:length(vzlist)    
    vz=vzlist(i);
    nn=zeros(1,length(omegalist));
    for j=1:length(omegalist)
        omega=omegalist(j);
        nn(j)=ldos(mu,Delta,vz,alpha,dim,smoothpot,mumax,peakpos,sigma,x,omega,delta);
    end
    en(i,:)=nn;
end
re=en;
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(Delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_smoothpot=num2str(smoothpot);
fn_mumax=strcat('mx',num2str(mumax));
fn_pos=strcat('x',num2str(x));
fn_sigma=strcat('sg',num2str(sigma));
fn_range=strcat('-',num2str(vzlist(end)),',',num2str(omegalist(end)),'-');
%  fn_delta=strcat('eta',num2str(delta));
if (strcmp(smoothpot,'lorentz')||strcmp(smoothpot,'lorentzsigmoid'))
    fn_peakpos=strcat('pk',num2str(peakpos));
else
    fn_peakpos='';
end
fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_smoothpot,fn_mumax,fn_peakpos,fn_pos,fn_sigma,fn_range);
save(strcat('LDOS',fn,'.dat'),'re','-ascii');
surf(vzlist,omegalist,en','edgecolor','none');colorbar;view(2);

saveas(gcf,strcat('LDOS',fn,'.png'))
end