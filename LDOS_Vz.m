function [rev,re]=LDOS_Vz(mu,dim,eta,smoothpot,mumax,peakpos)
delta=0.2;
alpha=5;
x=1;
% eta=1e-3;
vzlist=linspace(0.,0.4,100);
omegalist=linspace(-.3,.3,1001);
en=zeros(length(vzlist),length(omegalist));
for i=1:length(vzlist)    
    vz=vzlist(i);
    nn=zeros(1,length(omegalist));
    parfor j=1:length(omegalist)
        omega=omegalist(j);
        nn(j)=ldos(vz,x,omega,eta);
    end
    en(i,:)=nn;
end
re=en;
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_smoothpot=num2str(smoothpot);
fn_mumax=strcat('mx',num2str(mumax));
fn_pos=strcat('x',num2str(dim));
fn_eta=strcat('eta',num2str(eta));
if (strcmp(smoothpot,'lorentz')||strcmp(smoothpot,'lorentzsigmoid'))
    fn_peakpos=strcat('pk',num2str(peakpos));
else
    fn_peakpos='';
end
fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_smoothpot,fn_mumax,fn_peakpos,fn_pos,fn_eta);
save(strcat('LDOS',fn,'.dat'),'re','-ascii');
surf(vzlist,omegalist,en','edgecolor','none');colorbar;view(2);

saveas(gcf,strcat('LDOS',fn,'.png'))
end