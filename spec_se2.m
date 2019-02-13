%%for self energy with delta_0 collapse at vzc
function [rev,re]=spec_se2(a,mu,delta,vzc,alpha,gamma,dim)
% a=1;
vzlist=linspace(0,vzc*0.999,200);
nv=10;
en=zeros(nv,length(vzlist));
deltalist=delta*(1-(vzlist/vzc).^2);
parfor i=1:length(vzlist)
    vz=vzlist(i);
%     disp(i);
    deltac=deltalist(i);
    for n=1:nv
        en(n,i)=iter(a,mu,deltac,vz,alpha,gamma,n,dim);
    end
end
re=en;
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_gamma=strcat('g',num2str(gamma));
fn_vzc=strcat('vzc',num2str(vzc));

fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_gamma,fn_vzc);
save(strcat(fn,'.dat'),'re','-ascii');
plot(vzlist,en)
xlabel('V_Z(meV)')
ylabel('V_{bias}(meV)')
axis([0,vzlist(end),-.3,.3])
line([sqrt(mu^2+gamma^2),sqrt(mu^2+gamma^2)],[-0.3,0.3])
saveas(gcf,strcat(fn,'.png'))
end