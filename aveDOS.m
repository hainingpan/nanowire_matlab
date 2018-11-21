vzlist=[0.1,0.2,0.3,0.4];
disorderlist=[0,0.1,0.5,1,2,3,5,10,20,30];
for ii=1:length(vzlist)
    vz=vzlist(ii);
    dosstore=zeros(100,100);
    for jj=1:length(disorderlist)
        disorder=disorderlist(jj);
        for i=1:100
            dosstore(i,:)=DOS_Vzcut(.2,300,disorder,vz);
        end
        plot(linspace(-0.3,0.3,100),mean(dosstore));
        xlabel('E(meV)');
        ylabel('DOS');
        saveas(gcf,strcat(num2str(vz),'V',num2str(disorder),'.png'));
    end
end