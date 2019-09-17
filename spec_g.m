%%g is randomly varying along the nanowire
function [rev,re,randlist]=spec_g(a,mu,dim,sigma,randlist)
% a=1;
delta=0.2;
alpha=5;
vzlist=linspace(0,2.048,100);
nv=60;
en=zeros(nv,length(vzlist));
if randlist==0    
    randlist=(sigma*randn(dim,1)+1);
    while (nnz(randlist<0)~=0)
        randlist=(sigma*randn(dim,1)+1);
    end
end

parfor i=1:length(vzlist)
    vz=vzlist(i);
    vzrandlist=vz*randlist;
    ham=hg(a,mu,delta,vzrandlist,alpha,dim);
    eigo=eigs(ham,nv,0,'Tolerance',1e-5,'MaxIterations',20000);
    en(:,i)=sort(eigo(1:nv));
end
re=en;
rev=vzlist;
fn_mu=strcat('m',num2str(mu));
fn_Delta=strcat('D',num2str(delta));
fn_alpha=strcat('a',num2str(alpha));
fn_wl=strcat('L',num2str(dim));
fn_sigma=strcat('sg',num2str(sigma));

fn=strcat(fn_mu,fn_Delta,fn_alpha,fn_wl,fn_sigma);
save(strcat(fn,'.dat'),'re','-ascii');
figure;
plot(vzlist,en)
xlabel('V_Z(meV)')
ylabel('V_{bias}(meV)')
axis([0,vzlist(end),-.3,.3])
line([sqrt(mu^2+delta^2),sqrt(mu^2+delta^2)],[-0.3,0.3])
saveas(gcf,strcat(fn,'.png'))
end