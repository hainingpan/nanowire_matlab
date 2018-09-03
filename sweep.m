function sweep(deltac)
delta=.2;
mulist=0:.05:2;
% store=zeros(100,201);
for i=1:length(mulist)
    mu=mulist(i);
    disp(mu);
    [X,Y]=spec(mu,delta,2,deltac);
    save(strcat('mu',num2str(mu),'.dat','Y','-ascii'));
end
% save(strcat('deltac',num2str(deltac),'.dat'),'store','-ascii');
end
