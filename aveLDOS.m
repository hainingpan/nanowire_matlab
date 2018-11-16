dosstore=zeros(200,200);
for i=1:200
    dosstore(i,:)=LDOS_Vzcut(0.2,200,10,1,0.2);
end

dosstore2=zeros(200,200);
for i=1:200
    dosstore2(i,:)=LDOS_Vzcut(0.2,200,20,1,0.2);
end

dosstore3=zeros(200,200);
for i=1:200
    dosstore3(i,:)=LDOS_Vzcut(0.2,200,30,1,0.2);
end