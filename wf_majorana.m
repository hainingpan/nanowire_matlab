function [wf1,wf2,x1,x2,overlap]=wf_majorana(ham,L)
nv=10;
nsite=length(ham)/4;
[evec,eigo]=eigs(ham,nv,0,'Tolerance',1e-5,'MaxIterations',20000);
eigo=diag(eigo);
[eigoz,I]=mink(eigo,2,'ComparisonMethod','abs');
[~,larger]=max(eigoz);
[~,smaller]=min(eigoz);
vecpos=evec(:,I(larger));
vecnega=evec(:,I(smaller));
vec1=1/sqrt(2)*(vecpos+vecnega);
vec2=1i/sqrt(2)*(vecpos-vecnega);
wf1=sum(reshape(abs(vec1).^2,[],4),2);
wf2=sum(reshape(abs(vec2).^2,[],4),2);
x1=sum(linspace(0,L,nsite)'.*wf1);
x2=sum(linspace(0,L,nsite)'.*wf2);
overlap=sum(sqrt(wf1.*wf2));
end