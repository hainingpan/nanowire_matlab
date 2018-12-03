den= 1000;
vztest=linspace(0,0.6,den)';
mutest=linspace(0,0.6,den)';
[muGrid,vzGrid] = meshgrid(mutest,vztest);
prGrid = [muGrid(:),vzGrid(:)];
N = size(prGrid,1);
Scores = zeros(N,numel(classes));
for j = 1:numel(classes)
    [~,score] = predict(SVMModels{j},prGrid);
    Scores(:,j) = score(:,2); % Second column contains positive-class scores
end

[~,maxScore] = max(Scores,[],2);
figure;
h(1:3) = gscatter(prGrid(:,1),prGrid(:,2),maxScore,[0.1 0.5 0.5; 0.5 0.1 0.5; 0.5 0.5 0.1]);
hold on
h(4:6) = gscatter(pr(:,1),pr(:,3),label);
axis tight
