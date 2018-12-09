murange=[0,0.4];
deltarange=[0.2,0.2];
vzrange=[0,0.4];
[pr,x]=gen(murange,deltarange,vzrange); %Generate data

y=PCAn(x,2); %PCA

label=kmeans(y,3); %kmeans

kk=3;%% number of classes
%SVM
SVMModels = cell(kk,1);
classes = unique(label);

for j = 1:numel(classes)
    indx=(label==classes(j));
    SVMModels{j} = fitcsvm( pr(:,[1,kk]),indx,'ClassNames',[false true],'KernelFunction','rbf','Standardize',true);
end

%predict

den= 1000;
vztest=linspace(0,vzrange(2),den)';
mutest=linspace(0,murange(2),den)';
[muGrid,vzGrid] = meshgrid(mutest,vztest);
prGrid = [muGrid(:),vzGrid(:)];
N = size(prGrid,1);
Scores = zeros(N,numel(classes));
for j = 1:numel(classes)
    [~,score] = predict(SVMModels{j},prGrid);
    Scores(:,j) = score(:,2); % Second column contains positive-class scores
end

%Visualization
[~,maxScore] = max(Scores,[],2);
figure;
h(1:kk) = gscatter(prGrid(:,1),prGrid(:,2),maxScore,[0 0 1;1 1 0 ; 0.0 0.5 1]);
hold on
h(kk+1:2*kk) = gscatter(pr(:,1),pr(:,kk),label);
axis tight
