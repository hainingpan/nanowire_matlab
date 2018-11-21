%%DOS Vz cut
function re=DOS_Vzcut(mu,dim,v,vz)
Delta=0.2;
alpha=5;
delta=1e-3;
omegalist=linspace(-.3,.3,100);
nn=zeros(1,length(omegalist));
vimp=v*randn(dim,1);
parfor j=1:length(omegalist)
    omega=omegalist(j);
    nn(j)=dosdis(mu,Delta,vz,alpha,dim,vimp,omega,delta);
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
plot(omegalist,nn);
% saveas(gcf,strcat('LDOS',fn,'.png'))
end