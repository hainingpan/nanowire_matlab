function [para,wf]=gen(murange,deltarange,vzrange)
a=1;
num=10000;
mulist=(murange(2)-murange(1))*rand(num,1)+murange(1);
deltalist=(deltarange(2)-deltarange(1))*rand(num,1)+deltarange(1);
vzlist=(vzrange(2)-vzrange(1))*rand(num,1)+vzrange(1);
para=[mulist,deltalist,vzlist];
alpha_R=5;
dim=500;
sp=zeros(num,dim);
sp2=zeros(num,dim);
parfor i=1:num
    vz=vzlist(i);
    mu=mulist(i);
    delta=deltalist(i);
    wfabs=WF(a,mu,delta,vz,alpha_R,dim,1);
    sp(i,:)=wfabs.';
    wfabs2=WF(a,mu,delta,vz,alpha_R,dim,2);
    sp2(i,:)=wfabs2.';
end
wf=[sp,sp2];
end
