mulist=0:.05:2;
b=zeros(length(mulist),51);
u=zeros(length(mulist),51);
delta=.2;
parfor i=1:length(mulist)
    mu=mulist(i);
    b(i,:)=specbottom(mu,delta,2);
    u(i,:)=specupper(mu,delta,2);
end
