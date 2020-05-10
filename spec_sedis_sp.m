%%for self energy & disorder
function [dosmap,rev,vimp]=spec_sedis_sp(a,mu,delta,alpha,gamma,vc,dim,v,vimp)
% a=1;
vzlist=linspace(0,2.048,401);
% vzlist=0:0.001:1.05;
% nv=20;
if vimp==0
    vimp=v*randn(dim,1);
end
dosmap=cell(1,length(vzlist));
enlist=linspace(-.21,.21,1001);
parfor i=1:length(vzlist)
    vz=vzlist(i);
    disp(i);
    dos=arrayfun(@(w) dossedis(a,mu,delta,vz,alpha,gamma,vc,dim,vimp,w,1e-3),enlist);
    [~,loc]=findpeaks(dos);
    init=enlist(loc);
%     num_init=min(nv,length(init));
%     tmp=init(1:num_init);
%     if num_init<nv
%         tmp=[tmp,zeros(1,nv-num_init)];
%     end
    dc=delta*sqrt(1-(vz/vc)^2);
    dosmap{i}=init(abs(init)<dc); 
end
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_gamma=strcat('g',num2str(gamma));
fn_v=strcat('v',num2str(v));
fn_vc=strcat('vc',num2str(vc))*(vc~=inf);



fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_gamma,fn_v,fn_vc);
% save(strcat(fn,'.dat'),'dosmap','-ascii');

fid = fopen(strcat(fn,'.dat'),'w');
for i=1:length(vzlist)
    fprintf(fid,'%f ', dosmap{i});
    fprintf(fid,'\n');
end
fclose(fid);

% dosmap(dosmap==0)=nan;

figure;
for i=1:length(vzlist)
    scatter(ones(1,length(dosmap{i}))*vzlist(i),dosmap{i},'b','.');
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