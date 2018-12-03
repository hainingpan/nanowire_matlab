SVMModels = cell(3,1);
classes = unique(label);

for j = 1:numel(classes)
    indx=(label==classes(j));
    SVMModels{j} = fitcsvm( pr(:,[1,3]),indx,'ClassNames',[false true],'KernelFunction','rbf','Standardize',true);
end
