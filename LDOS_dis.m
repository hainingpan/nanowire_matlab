function LDOS_dis(a,mu,Delta,muVar,N_muVar,mulist,dim)
% mu=1;
% Delta=.2;
alpha_R=5;
delta=1e-3;
Vzlist=linspace(0,2.048,201);
energylist=linspace(-.3,.3,201);
ldosmap=zeros(length(Vzlist),length(energylist),dim);
lenVz=length(Vzlist);
lenenergy=length(energylist);
% muVar=1;
writematrix(mulist,'m1.txt')
if length(mulist)==1
    mulist=muVar*randn(N_muVar,1);
%     mulist=muVar*(2*rand(n,1)-1);
end
writematrix(mulist,'m2.txt')

mulist=arrayfun(@(x) mulist(floor((x-1)*N_muVar/dim)+1),1:dim);

writematrix(mulist,'m3.txt')

parfor i=1:lenVz
%     tic;
    for j=1:lenenergy
    Vz=Vzlist(i);
    energy=energylist(j);    
    ldosmap(i,j,:)=ldosall_dis(a,mu,Delta,Vz,alpha_R,mulist,dim,energy,delta);       
    end
%     toc;
end

figure;
LDOS_L=(squeeze(ldosmap(:,:,1)))';
surf(Vzlist/Delta,energylist/Delta,LDOS_L,'edgecolor','none');
view(2);
axis tight;
colorbar;colormap hot;
xline(sqrt(mu^2+Delta^2)/Delta,'g');
xlabel('V_z/\Delta (meV)');
ylabel('E/\Delta (meV)');
title(strcat('LDOS on the left end \mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu)));
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(Delta));
fn_muVar=strcat('muVar',num2str(muVar));
fn_L=strcat('L',num2str(dim*a));
fn=strcat(fn_mu,fn_Delta,fn_muVar,fn_L,'_LDOS_L');
saveas(gcf,strcat(fn,'.png'));

figure;
LDOS_M=(squeeze(ldosmap(:,:,floor(dim/2))))';
surf(Vzlist/Delta,energylist/Delta,LDOS_M,'edgecolor','none');
view(2);
axis tight;
colorbar;colormap hot;
xline(sqrt(mu^2+Delta^2)/Delta,'g');
xlabel('V_z/\Delta');
ylabel('E/\Delta');
title(strcat('LDOS in the bulk \mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu)));
fn=strcat(fn_mu,fn_Delta,fn_muVar,fn_L,'_LDOS_M');
saveas(gcf,strcat(fn,'.png'));

figure;
LDOS_R=(squeeze(ldosmap(:,:,end)))';
surf(Vzlist/Delta,energylist/Delta,LDOS_R,'edgecolor','none');
view(2);
axis tight;
colorbar;colormap hot;
xline(sqrt(mu^2+Delta^2)/Delta,'g');
xlabel('V_z/\Delta');
ylabel('E/\Delta');
title(strcat('LDOS on the right end\mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu)));
fn=strcat(fn_mu,fn_Delta,fn_muVar,fn_L,'_LDOS_R');
saveas(gcf,strcat(fn,'.png'));

figure;
DOS=(mean(ldosmap,3))';
surf(Vzlist/Delta,energylist/Delta,DOS,'edgecolor','none');
view(2);
axis tight;
colorbar;colormap hot;
xline(sqrt(mu^2+Delta^2)/Delta,'g');
xlabel('V_z/\Delta');
ylabel('E/\Delta');
title(strcat('DOS \mu/\Delta=',num2str(mu/Delta),',\sigma_\mu/\mu=',num2str(muVar/mu)));
fn=strcat(fn_mu,fn_Delta,fn_muVar,fn_L,'_DOS');
saveas(gcf,strcat(fn,'.png'));

fn=strcat(fn_mu,fn_Delta,fn_muVar,fn_L);
save(strcat(fn,'.mat'),'LDOS_L','LDOS_M','LDOS_R','DOS','mulist','Vzlist','energylist');
