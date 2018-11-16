%%LDOS Vz cut
function re=LDOS_Vzcut(mu,dim,numdis,x,vz)
Delta=0.2;
alpha=5;
delta=1e-3;
omegalist=linspace(-.3,.3,200);
nn=zeros(1,length(omegalist));
pos=randperm(dim,numdis);
parfor j=1:length(omegalist)
    omega=omegalist(j);
    nn(j)=ldosdis(mu,Delta,vz,alpha,dim,pos,x,omega,delta);
end
re=nn;
rev=omegalist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(Delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_range=strcat('-',num2str(omegalist(end)),'-');
fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_range);
% save(strcat('LDOS',fn,'.dat'),'re','-ascii');
% plot(omegalist,nn);
% saveas(gcf,strcat('LDOS',fn,'.png'))
end