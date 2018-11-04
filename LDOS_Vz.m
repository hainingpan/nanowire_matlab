function [rev,re]=LDOS_Vz(mu,dim,smoothpot,mumax,peakpos)
delta=0.2;
alpha=5;
vzlist=0:0.001:0.4;
omegalist=-0.3:0.001:0.3;
en=zeros(length(vzlist),length(omegalist));
for i=1:length(vzlist)    
    vz=vzlist(i);
    nn=zeros(1,length(omegalist));
    for j=1:length(omegalist)
        omega=omegalist(j);
        nn(j)=ldos(vz,1,omega,1e-5);
    end
    en(i,:)=nn;
end
re=en;
surf(vzlist,omegalist,en','edgecolor','none');colorbar;view(2)
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_smoothpot=num2str(smoothpot);
fn_mumax=strcat('mx',num2str(mumax));
fn_pos=strcat('x',num2str(x));
if (strcmp(smoothpot,'lorentz')||strcmp(smoothpot,'lorentzsigmoid'))
    fn_peakpos=strcat('pk',num2str(peakpos));
else
    fn_peakpos='';
end
fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_smoothpot,fn_mumax,fn_peakpos,fn_pos);
save(strcat(fn,'.dat'),'re','-ascii');
plot(vzlist,en)
xlabel('V_Z(meV)')
ylabel('V_{bias}(meV)')
axis([0,4,-.3,.3])
saveas(gcf,strcat(fn,'.png'))
end