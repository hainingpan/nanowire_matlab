function [y,d]=PCAn(data,n)
norsp=normalize(data);
sum=0;
m=size(data);
parfor i=1:m(1)
    sum=sum+norsp(i,:).'*norsp(i,:);
end
sum=sum/m(1);
[v,d]=eigs(sum,n);
y=norsp*v;
end
