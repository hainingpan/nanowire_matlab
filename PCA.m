norsp=normalize(sp);
sum=0;
parfor i=1:length(norsp)
    sum=sum+norsp(i,:).'*norsp(i,:);
end
sum=sum/length(norsp);
[v,d]=eigs(sum,1);
y=norsp*v;
figure;
% scatter(y(:,1),y(:,2),[],vzlist)
scatter(vzlist,y);
% scatter(y(:,1),y(:,2));
