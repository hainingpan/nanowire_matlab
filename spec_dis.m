%%spectrum for disorder
function [rev,re]=spec_dis(a,mu,dim,numdis,v)
% a=1;
delta=0.2;
alpha=5;
vzlist=0:0.01:2;
nv=60;
en=zeros(nv,length(vzlist));
pos=randperm(dim,numdis);
for i=1:length(vzlist)
    vz=vzlist(i);
    ham=hdis(a,mu,delta,vz,alpha,dim,pos,v);
    eigo=eigs(ham,nv,0,'Tolerance',1e-5,'MaxIterations',20000);
    en(:,i)=sort(eigo(1:nv));
end
re=en;
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_num=strcat('N',num2str(numdis));
fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_num);
save(strcat(fn,'.dat'),'re','-ascii');
plot(vzlist,en)
xlabel('V_Z(meV)')
ylabel('V_{bias}(meV)')
axis([0,vzlist(end),-.3,.3])
saveas(gcf,strcat(fn,'.png'))
end