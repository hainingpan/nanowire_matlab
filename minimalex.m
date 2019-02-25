xrange=linspace(-2,2,100);
Y=repmat(-xrange.^2,4,1)+repmat((-4:-1)',1,100);
Y(Y<-5)=0;
for i=1:100
    [~,~,v]=find(Y(:,i));
    Y(1:length(v),i)=v;
end
Y(Y==0)=nan;
%jump due to missing data
figure;
plot(xrange,Y);

figure;
%from bare eye, we see there are four lines
for i=1:4
    scatter(xrange,Y(i,:),'b');
    hold on
end