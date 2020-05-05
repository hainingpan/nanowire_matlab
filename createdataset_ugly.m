dir='D:\CMTC\bothlead\Rp30\L500\';
datastore=[];
mulist=0.1:0.1:1;
for index=1:100
    fn=strcat(dir,num2str(index));
    for muindex=1:length(mulist)
        muVar=mulist(muindex);
        L=300;
        fnL=sprintf("m1.0D0.2a5L%dmVar%.1fLL-2.048,0.3-",L,muVar);
        fnR=sprintf("m1.0D0.2a5L%dmVar%.1fRR-2.048,0.3-",L,muVar);
        muVarfn=sprintf("muVar%.1flist.dat",muVar);        
        dataL=load(strcat(fn,'\',fnL,'.dat'));
        dataR=load(strcat(fn,'\',fnR,'.dat'));
        muVarlist=load(strcat(fn,'\',muVarfn));
        [nVz,ne]=size(dataL);
        vzmax=2.048;
        vzlist=linspace(0,vzmax,nVz);
        emax=0.3;
        elist=linspace(-emax,emax,ne);
        ispeakL=zeros(1,nVz);
        ispeakR=zeros(1,nVz);
        condL=zeros(1,nVz);
        condR=zeros(1,nVz);
        tol=0.5e-2;
        for j=1:nVz
            [ispeakL(j),condL(j)]=ispeak(elist(ceil(ne/2)-100:ceil(ne/2)+100),dataL(j,ceil(ne/2)-100:ceil(ne/2)+100),tol);
            [ispeakR(j),condR(j)]=ispeak(elist(ceil(ne/2)-100:ceil(ne/2)+100),dataR(j,ceil(ne/2)-100:ceil(ne/2)+100),tol);
        end
        peakindex=(logical(ispeakL)&logical(ispeakR));
        vzpeak=vzlist(peakindex);
        condL=condL(peakindex);
        condR=condR(peakindex);
        mu=1;
        Delta=0.2;
        for j=1:length(vzpeak)
            vz=vzpeak(j);
            ham=hdis(1,mu,Delta,vz,5,L,muVarlist);
            [wf1,wf2,x1,x2,overlap]=wf_majorana(ham,L);
            dp=[condL(j),condR(j),x1,x2,overlap,index,muVar,vz];
            datastore=[datastore;dp];
        end    
    end
end
%distance
figure;scatter3(datastore(:,1),datastore(:,2),abs(datastore(:,4)-datastore(:,3)),[],abs(datastore(:,4)-datastore(:,3)),'.')
%overlap
figure;scatter3(datastore(:,1),datastore(:,2),datastore(:,5),[],datastore(:,5),'.')
    