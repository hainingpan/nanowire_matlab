function sweep(delta)
mulist=0:.05:2;
store=zeros(length(mulist),50);
parfor i=1:length(mulist)
    mu=mulist(i);
    store(i,:)=spec(mu,delta,2);
end
save(strcat('deltac0.2',num2str(delta),'.dat'),'store','-ascii');
end
