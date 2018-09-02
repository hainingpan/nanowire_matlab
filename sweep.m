function sweep(deltac)
delta=.2;
mulist=0:.05:2;
store=zeros(length(mulist),101);
for i=1:length(mulist)
    mu=mulist(i);
    store(i,:)=spec(mu,delta,2,deltac);
end
% save(strcat('deltac',num2str(deltac),'.dat'),'store','-ascii');
end
