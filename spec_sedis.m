%%for self energy & disorder
function [rev,re,dosmap,vimp]=spec_sedis(a,mu,delta,alpha,gamma,vc,dim,v,vimp)
% a=1;
vzlist=linspace(0,2,101);
nv=4;
en=zeros(nv,length(vzlist));
if vimp==0
    vimp=v*randn(dim,1);
end
enlist=linspace(0,.2,201);
dosmap=zeros(nv,length(vzlist));
for i=1:length(vzlist)
    vz=vzlist(i);
    disp(i);
    dos=arrayfun(@(w) dossedis(1,.2,vz,5,.2,3,300,vimp,w,1e-3),enlist);
    [pk,loc]=findpeaks(dos);
    init=enlist(loc);
    dosmap(:,i)=init(1:4);
    for n=1:nv        
        if n<=length(init)
            en(n,i)=iter_dis(a,mu,delta,vz,alpha,gamma,vc,n,dim,vimp,init(n));
        end
    end
end
re=en;
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_gamma=strcat('g',num2str(gamma));
fn_v=strcat('v',num2str(v));
fn_vc=strcat('vc',num2str(vc));


fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_gamma,fn_v,fn_vc);
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