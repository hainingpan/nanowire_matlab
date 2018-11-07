function LDOS_x(mu,dim,smoothpot,mumax,peakpos,sigma,vz)
Delta=0.2;
alpha=5;
delta=1e-3;
omegalist=linspace(-.3,.3,200);
en=zeros(length(omegalist),dim);
parfor i=1:length(omegalist)
    omega=omegalist(i);
    nn=ldosall(mu,Delta,vz,alpha,dim,smoothpot,mumax,peakpos,sigma,omega,delta);
    en(i,:)=nn;
end
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(Delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_smoothpot=num2str(smoothpot);
fn_mumax=strcat('mx',num2str(mumax));
fn_vz=strcat('v',num2str(vz));
fn_sigma=strcat('sg',num2str(sigma));
fn_range=strcat('-',num2str(dim),',',num2str(omegalist(end)),'-');
if (strcmp(smoothpot,'lorentz')||strcmp(smoothpot,'lorentzsigmoid'))
    fn_peakpos=strcat('pk',num2str(peakpos));
else
    fn_peakpos='';
end
fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_smoothpot,fn_mumax,fn_peakpos,fn_vz,fn_sigma,fn_range);
save(strcat('LDOS',fn,'.dat'),'en','-ascii');
surf((0:dim-1)/100,omegalist,en,'edgecolor','none');colorbar;view(2);
axis tight;
saveas(gcf,strcat('LDOS',fn,'.png'))
end



