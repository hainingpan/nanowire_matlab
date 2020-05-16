%%for self energy & quantum dots
function [dosmap,rev]=spec_seqddis_sp(a,mu,delta,alpha,gamma,vc,mumax,l0,dim,v,vimp,period)
% a=1;
vzlist=linspace(0,2.048,401);
%vzlist=0:0.0025:0.95;

% nv=20;
if vimp==0
    vimp=v*randn(dim,1);
end
dosmap=cell(1,length(vzlist));
enlist=linspace(-.21,.21,1001);
% dosmap2=zeros(length(vzlist),length(enlist));
parfor i=1:length(vzlist)
    vz=vzlist(i);
    disp(i);    
    dos=arrayfun(@(w) dosseqddis(a,mu,delta,vz,alpha,gamma,vc,mumax,l0,dim,vimp,w,1e-3,period),enlist);
    [~,loc]=findpeaks(dos);
    init=enlist(loc);
%     dosmap2(i,:)=dos;

% %     num_init=min(nv,length(init));
% %     tmp=init(1:num_init);
% %     if num_init<nv
% %         tmp=[tmp,zeros(1,nv-num_init)];
% %     end

%     func=@(w) log(dosseqd(a,mu,delta,vz,alpha,gamma,vc,mumax,l0,dim,w,1e-3)+100);
%     engrid=linspace(-.3,.3,11);
%     x=[];
%     y=[];
%     for j=1:10
%         opt.min=engrid(j);
%         opt.max=engrid(j+1);
%         opt.maxrecursion=3;
%         opt.points=12;
%         [xlist,ylist]=adaptive(func,opt);
%         x=[x,xlist];
%         y=[y,ylist];
%     end
    
%     [xorder,order]=sort(x);
%     yorder=y(order);
%     [pks1,loc1]=findpeaks(yorder);
%     init=xorder(loc1);
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
fn_mumax=strcat('mx',num2str(mumax));
fn_l0=strcat('l',num2str(l0));

fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_mumax,fn_l0,fn_gamma,fn_v,fn_vc);
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
axis([0,vzlist(end),enlist(1),enlist(end)])
line([sqrt(mu^2+gamma^2),sqrt(mu^2+gamma^2)],[enlist(1),enlist(end)])
saveas(gcf,strcat(fn,'.png'))
end