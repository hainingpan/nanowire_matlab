mulist=0:.05:2;
b=zeros(length(mulist),50);
u=zeros(length(mulist),50);
delta=.2;
parfor i=1:length(mulist)
    mu=mulist(i);
    b(i,:)=specbottom(mu,delta,2);
    u(i,:)=specupper(mu,delta,2);
end
