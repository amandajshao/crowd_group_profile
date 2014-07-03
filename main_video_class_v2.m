%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% main function for bag-of-word (BoW) codebook learning %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

clc;clear;

%% path setting
path = '.\';
path_descr = [path, 'result_groupDescr_new\'];

load([path, 'gt_labels.mat'], 'training_labels');
load([path, 'gt_videoName_bool.mat'], 'videoName_bool');


%% load data
load([path_descr, 'collectiveness.mat'], 'collectiveness');
load([path_descr, 'uniform.mat'], 'uniform');
load([path_descr, 'stability.mat'], 'stability');
load([path_descr, 'groupSize.mat'], 'group_size');
load([path_descr, 'CodeVector_kmeans_max_pooling.mat'], 'CodeVector');
conflict = CodeVector';

X1_temp = [[collectiveness.coll_v_mean]', [collectiveness.coll_v_var]', ...
    reshape([collectiveness.coll_count],[2,474])'];
X1 = zeros(size(collectiveness,1),3);
X1(:,1) = (X1_temp(:,1)-min(X1_temp(:,1))) / (max(X1_temp(:,1))-min(X1_temp(:,1)));
X1(:,2) = (X1_temp(:,2)-min(X1_temp(:,2))) / (max(X1_temp(:,2))-min(X1_temp(:,2)));
X1((X1_temp(:,3)./(X1_temp(:,4)+eps))>3, 3) = 1; 
X1((X1_temp(:,3)./(X1_temp(:,4)+eps))<0.45, 3) = -1; 

X2 = [[uniform.unif_mean]', [uniform.unif_var]'];
X2(:,1) = (X2(:,1)-min(X2(:,1))) / (max(X2(:,1))-min(X2(:,1)));
X2(:,2) = (X2(:,2)-min(X2(:,2))) / (max(X2(:,2))-min(X2(:,2)));

X3 = [stability.stab_invKnnNum]';
X3 = (X3-min(X3)) / (max(X3)-min(X3));
X4 = reshape([stability.stab_rankKnn],[10,474])';
X4 = (X4-min(min(X4))) / (max(max(X4))-min(min(X4)));
X5 = reshape([stability.stab_rwHist],[201,474])';
X5 = (X5-min(min(X5))) / (max(max(X5))-min(min(X5)));

X6 = conflict;
X6 = (X6-min(min(X6))) / (max(max(X6))-min(min(X6)));

X7 = reshape([group_size.gr_size_max_v],[6,474])';
X7 = (X7-min(min(X7))) / (max(max(X7))-min(min(X7)));



%%
% feature combination
X = [X1, X2, X3, X4, X5, X6, X7]; 
nGroup = size(X, 1);
label_unique = unique(training_labels);
test_num = 50;
perc = 1;
svm_n = 0;% 0-nonlinear(libsvm); 1-linear(liblinear)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TODO: leave-one-out svm: each time, select one scene to test its
%       independency. That is, first select a sample, and delete samples with
%       the same scene as this sample; and then for the other classes,
%       randomly select a sample from each class to build the final test set.
for test_n = 1 : test_num
    num = 0;
    param_alpha = 1e1;
    while param_alpha <= 1e8
        num = num + 1;
        param_alpha = param_alpha * 1e1;
        
        % generate data index for training
        for n = 1 : length(label_unique)
            class_num(n) = length(find(training_labels == label_unique(n)));
        end
                
        % random select one out to test
        rndIdx = []; rndIdx_1 = [];
        for n = 1 : length(label_unique)
            K = round(class_num(n) * perc);
            idx_rand = randperm(class_num(n));
            idx_curGr = find(training_labels == label_unique(n));
%             rndIdx_1 = [rndIdx_1; idx_curGr(idx_rand(2:K))]; % leave 1st in idx_rand
%             rndIdx = [rndIdx; idx_curGr(idx_rand(2:min(70,K)))];
            rndIdx = [rndIdx; idx_curGr(idx_rand(2:K))];
        end
        restIdx = setdiff(1:nGroup, rndIdx);
        
        % if there is other video clips from the same scene exists, remove them from training set
        for n = 1 : length(restIdx)
            video_bool = videoName_bool(restIdx(n),2);
            [class_idx,~] = find(videoName_bool(:,2) == video_bool);
            if length(class_idx) > 1
                [~,delete_idx] = intersect(rndIdx, class_idx);
                rndIdx(delete_idx) = [];
            end
        end
        
        % classification
        train_data = X(rndIdx, :);
        train_label = training_labels(rndIdx);
        test_data = X(restIdx, :);
        test_label = training_labels(restIdx);
        options = ['-c ' num2str(param_alpha) ' -s ' num2str(0) ' -t ' num2str(2) ' -b ' num2str(1)];
        
        model = svmtrain(train_label, train_data, options);
        [class_label, accuracy, q] = svmpredict(test_label, test_data, model);
        class_accuracy(num, :) = [svm_n, perc, param_alpha, accuracy(1)];
        
        if num == 1
            max_accuracy_para = class_accuracy(num,2:4); % [percentage, param_alpha, accruacy]
            label_test = [restIdx', test_label, class_label]; % [test_id, gt_label, test_label]
            label_train = [rndIdx, train_label];
        elseif class_accuracy(num,4) > max_accuracy_para(3)
            max_accuracy_para = class_accuracy(num,2:4);
            label_test = [restIdx', test_label, class_label]; % [test_id, gt_label, test_label]
            label_train = [rndIdx, train_label];
        end
    end
    record{test_n,1}.train = label_train;
    record{test_n,1}.test = label_test; % [test_id, gt_label, test_label]
    [C,order] = confusionmat(record{test_n,1}.test(:,2),record{test_n,1}.test(:,3));
    confusion_matrix = C./repmat(sum(C,2),1,size(C,1));
    accuracy_mean = mean(diag(confusion_matrix));
    record{test_n,1}.acc = [max_accuracy_para, accuracy_mean];
end

%% confusion matrix
C_sum = 0;
for i = 1 : length(record)
    [C,order] = confusionmat(record{i,1}.test(:,2), record{i,1}.test(:,3));
    C_sum = C_sum + C;
end
C_sum_percent = C_sum./test_num;
class_name = [{'1'}, {'2'}, {'3'}, {'4'}, {'5'}, {'6'}, {'7'}, {'8'}];
figure
draw_cm(C_sum_percent, class_name, length(class_name))
mean_acc = mean(diag(C_sum_percent));
fprintf('Average accuracy: %.4f\n', mean_acc);

%% 
% save([path,'v_class_record_old.mat'], 'record')




