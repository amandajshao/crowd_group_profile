function [classificationRes, accuracy, q] = svm_classifying(model, descriptorMatrix, classLabels, isLinear)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. the structure of descriptorMatrix is N*D, the number of descriptor is N,
% the length of one single descriptor is D
[Num_descriptor, ~] = size(descriptorMatrix);
if length(classLabels) ~= Num_descriptor
    error('number of labels and descriptors should be the same....');
end

if isLinear
    fprintf('Linear SVM is used .....\n');
%     addpath(genpath('liblinear-1.93'));
    
    if size(classLabels)==1
        classLabels = classLabels';
    end
    
    q = [];
    [classificationRes, accuracy] = predict(classLabels, sparse(descriptorMatrix), model);
else
    fprintf('Nonlinear SVM is used .....\n');
%     addpath(genpath('liblinear-1.93'));
    
    if size(classLabels)==1
        classLabels = classLabels';
    end
    
    [classificationRes, accuracy, q] = svmpredict(classLabels, sparse(descriptorMatrix), model, ['-b ' num2str(1)]);
end

end