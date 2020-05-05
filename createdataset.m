dir='D:\CMTC\bothlead\Rp29\inhom\L3noSE\';
datastore=[];
for index=1:610
    fn=strcat(dir,num2str(index));
    lspng=ls(strcat(fn,'\*.png'));
    if contains(lspng(1,:),'L-')
        fnL=lspng(1,:);
        fnR=lspng(2,:);
    else
        fnL=lspng(2,:);
        fnR=lspng(1,:);
    end
    scan=sscanf(fnL,"m%fD%fa%fL%fexpmx%fsg%fL-%f,%f-.png");
    mu=scan(1);
    Delta=scan(2);
    a=scan(3);
    L=scan(4);
    mumax=scan(5);
    sigma=scan(6);
    vzmax=scan(7);
    emax=scan(8);    
    dataL=load(strcat(fn,'\',strrep(fnL,'png','dat')));
    dataR=load(strcat(fn,'\',strrep(fnR,'png','dat')));
    [nVz,ne]=size(dataL);
    vzlist=linspace(0,vzmax,nVz);
    elist=linspace(-emax,emax,ne);
    ispeakL=zeros(1,nVz);
    ispeakR=zeros(1,nVz);
    condL=zeros(1,nVz);
    condR=zeros(1,nVz);
    tol=.5e-2;
    for j=1:nVz
        [ispeakL(j),condL(j)]=ispeak(elist(ceil(ne/2)-100:ceil(ne/2)+100),dataL(j,ceil(ne/2)-100:ceil(ne/2)+100),tol);
        [ispeakR(j),condR(j)]=ispeak(elist(ceil(ne/2)-100:ceil(ne/2)+100),dataR(j,ceil(ne/2)-100:ceil(ne/2)+100),tol);
    end
    peakindex=(logical(ispeakL)&logical(ispeakR));
    vzpeak=vzlist(peakindex);
    condL=condL(peakindex);
    condR=condR(peakindex);
%     condL=dataL(peakindex,ceil(ne/2));
%     condR=dataR(peakindex,ceil(ne/2));
    
    for j=1:length(vzpeak)
        vz=vzpeak(j);
        ham=hmu(1,mu,Delta,vz,5,L,'exp',mumax,0,sigma);
        [~,~,x1,x2,overlap]=wf_majorana(ham,L);
        dp=[condL(j),condR(j),x1,x2,overlap,vz,index];
        datastore=[datastore;dp];
    end    
end
%distance
figure;scatter3(datastore(:,1),datastore(:,2),abs(datastore(:,4)-datastore(:,3)),[],abs(datastore(:,4)-datastore(:,3)),'.')
%overlap
figure;scatter3(datastore(:,1),datastore(:,2),datastore(:,5),[],datastore(:,5),'.')
    