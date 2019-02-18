%%for self energy & inhomogenous potential
function [rev,re]=spec_seinhom(a,mu,delta,alpha,gamma,vc,dim,smoothpot,mumax,peakpos,sigma)
% a=1;
vzlist=linspace(0,2,101);
nv=4;
en=zeros(nv,length(vzlist));

parfor i=1:length(vzlist)
    vz=vzlist(i);
%     disp(i);
    for n=1:nv
        en(n,i)=iter_inhom(a,mu,delta,vz,alpha,gamma,vc,n,dim,smoothpot,mumax,peakpos,sigma);
    end
end
re=en;
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_gamma=strcat('g',num2str(gamma));
fn_vc=strcat('vc',num2str(vc))*(vc~=inf);
fn_smoothpot=num2str(smoothpot);
fn_mumax=strcat('mx',num2str(mumax));
fn_sigma=strcat('sg',num2str(sigma));
if (strcmp(smoothpot,'lorentz')||strcmp(smoothpot,'lorentzsigmoid'))
    fn_peakpos=strcat('pk',num2str(peakpos));
else
    fn_peakpos='';
end


fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_smoothpot,fn_mumax,fn_sigma,fn_peakpos,fn_gamma,fn_vc);
save(strcat(fn,'.dat'),'re','-ascii');
figure;
plot(vzlist,en)
hold on 
plot(vzlist,-en)
xlabel('V_Z(meV)')
ylabel('V_{bias}(meV)')
axis([0,vzlist(end),-.3,.3])
line([sqrt(mu^2+gamma^2),sqrt(mu^2+gamma^2)],[-0.3,0.3])
saveas(gcf,strcat(fn,'.png'))
end