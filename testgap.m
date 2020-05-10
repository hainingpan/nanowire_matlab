[vznum,enum]=size(dosmat);
enlist=linspace(-.3,.3,enum);
vzlist=0:0.0025:2.048;
for i=1:vznum
    [~,loc]=findpeaks(dosmat(i,:),'MinPeakProminence',5);
    init=enlist(loc);
    dos{i}=init;
end

figure;
hold on

for i=1:length(vzlist)
    scatter(ones(1,length(dos{i}))*vzlist(i),dos{i},'b','.');
end

box on
xlim([0,2.048])

