function model = svm_learning(descriptorMatrix, training_labels, isLinear, alpha)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. the structure of descriptorMatrix is N*D, the number of descriptor is N,
% the length of one single descriptor is D
% 2. training_labels is a vector 

if isLinear
    fprintf('Linear SVM is used .....\n');
%     addpath(genpath('liblinear-1.93'));
    
%     alpha = 1000;
    options = ['-c ' num2str(alpha) ' -B ' num2str(1) ' -s ' num2str(4)];
    model = train(double(training_labels), sparse(descriptorMatrix), options);
    
    w = model.w(:, 1:end-1);
    b = model.w(:, end);
    
%     save('w_linear.mat', 'w');
%     save('b_linear.mat', 'b');
%     save('model_linear.mat', 'model');
else
    % seems not necessary
%     addpath(genpath('libsvm-3.17'));
%     alpha = 1e3;
    options = ['-c ' num2str(alpha) ' -s ' num2str(0) ' -t ' num2str(2) ' -b ' num2str(1)];
    model = svmtrain(double(training_labels), sparse(descriptorMatrix), options);
    
 
end

end