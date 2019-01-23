%%spectrum for quantum dot
function [rev,re]=spec_qd(a,mu,dim,mumax,l0)
% a=1;
delta=0.2;
alpha=5;
vzlist=linspace(0,2.048*1.5,100);
nv=80;
en=zeros(nv,length(vzlist));

for i=1:length(vzlist)
    vz=vzlist(i);
    ham=hqd(a,mu,delta,vz,alpha,mumax,l0,dim);
    eigo=eigs(ham,nv,0,'Tolerance',1e-5,'MaxIterations',20000);
    en(:,i)=sort(eigo(1:nv));
end
re=en;
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_mumax=strcat('mx',num2str(mumax));
fn_l0=strcat('l',num2str(l0));
% fn_num=strcat('N',num2str(numdis));
fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_mumax,fn_l0,fn_wl);
save(strcat(fn,'.dat'),'re','-ascii');
plot(vzlist,en)
xlabel('V_Z(meV)')
ylabel('V_{bias}(meV)')
axis([0,vzlist(end),-.3,.3])
line([sqrt(mu^2+delta^2),sqrt(mu^2+delta^2)],[-0.3,0.3])
saveas(gcf,strcat(fn,'.png'))
end