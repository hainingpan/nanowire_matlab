mulist=0:.1:2;
store=zeros(length(mulist),50);
parfor i=1:length(mulist)
    mu=mulist(i);
    store(i,:)=spec(mu,.2,2);
end
save('store.dat','store','-ascii');
