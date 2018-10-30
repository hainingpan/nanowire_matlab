%%for single band 
function [rev,re]=spec0(mu,dim,smoothpot,mumax,peakpos)
a=.1;
delta=0.2;
alpha=5;
vzlist=0:0.01:4;
nv=30;
en=zeros(nv,length(vzlist));
for i=1:length(vzlist)
    vz=vzlist(i);
    ham=hmu(a,mu,delta,vz,alpha,dim,smoothpot,mumax,peakpos);
    eigo=eigs(ham,nv,0,'Tolerance',1e-5,'MaxIterations',20000);
    en(:,i)=sort(eigo(1:nv));
end
re=en;
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_smoothpot=num2str(smoothpot);
fn_mumax=strcat('mx',num2str(mumax));
if (strcmp(smoothpot,'lorentz')||strcmp(smoothpot,'lorentzsigmoid'))
    fn_peakpos=strcat('pk',num2str(peakpos));
else
    fn_peakpos='';
end
fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_smoothpot,fn_mumax,fn_peakpos);
save(strcat(fn,'.dat'),'re','-ascii');
plot(vzlist,en)
xlabel('V_Z(meV)')
ylabel('V_{bias}(meV)')
axis([0,4,-.3,.3])
saveas(gcf,strcat(fn,'.png'))
return 
