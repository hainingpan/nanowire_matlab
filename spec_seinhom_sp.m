%%for self energy & inhomogenous potential
function [dosmap,rev]=spec_seinhom_sp(a,mu,delta,alpha,gamma,vc,dim,smoothpot,mumax,peakpos,sigma)
% a=1;
vzlist=linspace(0,2,201);
nv=20;
en=zeros(nv,length(vzlist));

parfor i=1:length(vzlist)
    vz=vzlist(i);
%     disp(i);
    enlist=linspace(-.3,.3,401);
    dos=arrayfun(@(w) dosseinhom(a,mu,delta,vz,alpha,gamma,vc,dim,smoothpot,mumax,peakpos,sigma,w,1e-3),enlist);
    [~,loc]=findpeaks(dos);
    init=enlist(loc);
    num_init=min(nv,length(init));
    tmp=init(1:num_init);
    if num_init<nv
        tmp=[tmp,zeros(1,nv-num_init)];
    end
    dosmap(:,i)=tmp(:);    
end
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
save(strcat(fn,'.dat'),'dosmap','-ascii');

dosmap(dosmap==0)=nan;
figure;
for i=1:nv
    scatter(vzlist,dosmap(i,:),'b','.');
    hold on
end
box on
hold off
xlabel('V_Z(meV)')
ylabel('V_{bias}(meV)')
axis([0,vzlist(end),-.3,.3])
line([sqrt(mu^2+gamma^2),sqrt(mu^2+gamma^2)],[-0.3,0.3])
saveas(gcf,strcat(fn,'.png'))
end