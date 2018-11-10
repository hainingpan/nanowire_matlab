norsp=normalize(sp);
sum=0;
m=length(vzlist);
parfor i=1:m
%     disp(i);
    sum=sum+norsp(i,:).'*norsp(i,:);
end
sum=sum/m;
[v,d]=eigs(sum,2);
y=norsp*v;
% figure;
% scatter(y(:,1),y(:,2),[],vzlist)
% scatter(vzlist,y);
% scatter(y(:,1),y(:,2));
