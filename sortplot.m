[xorder,order]=sort(x);
yorder=y(order);
plot(xorder,yorder,'-o');

xg=linspace(tt.min,tt.max,length(x));
yg=func(xg)
hold on;
plot(xg,yg,'-o');
hold off;