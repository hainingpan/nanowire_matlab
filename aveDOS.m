vzlist=[0.1,0.2,0.3,0.4];
for ii=1:4
    vz=vzlist(ii);
dosstore=zeros(100,100);
for i=1:100
    dosstore(i,:)=DOS_Vzcut(.2,300,0.1,vz);
end
plot(mean(dosstore));
saveas(gcf,strcat(num2str(vz),'V0.1.png'));

dosstore2=zeros(100,100);
for i=1:100
    dosstore2(i,:)=DOS_Vzcut(.2,300,0.5,vz);
end
plot(mean(dosstore2));
saveas(gcf,strcat(num2str(vz),'V0.5.png'));

dosstore3=zeros(100,100);
for i=1:100
    dosstore3(i,:)=DOS_Vzcut(.2,300,1,vz);
end
plot(mean(dosstore3));
saveas(gcf,strcat(num2str(vz),'V1.png'));
end