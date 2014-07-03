%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% main function for bag-of-word (BoW) codebook learning %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

clc;clear;close all

%% Initialization & parameter setting
%==========================================================================
path = '.\result_groupDescr_new\';
param = struct;

% load training data
% conf_dict (62d) : [group_number(1d); conf_shape_context_of_each point(60d); file_number(1d)]
load([path, 'conf_dict.mat'], 'conf_dict_data');
conf_dict = double(conf_dict_data);
param.nsample = max(conf_dict(:,end));

% parameters for descriptor extraction
param.codebook = 'kmeans'; % methods for codebook construction
% param.codebook = 'sparse_coding';
param.num_basis = 1000; 
param.hist = 'max_pooling'; % histogram construction
% param.hist = 'sum_pooling';
if strcmp(param.codebook, 'sparse_coding')
    % define sparse coding parameters
    param.sc.func   = 'L1';
    param.sc.trials = 50;
    param.sc.beta   = 1e-1;
    param.sc.batch  = 5000;
    param.sc.gamma  = 0.15;
else
    param.knn = 10;
end


%==========================================================================
%% CodeBook computation
descripMtx = conf_dict(:,2:end-1)'; % descripMtx is a P \times N matrix
% descripMtx = descripMtx./max(max(descripMtx)); % normalize
fprintf('codebook computing begins......\n');
CodeBook = fun_codebook_computation(param, descripMtx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% you are recommended to save codebook before codevector computation,
% because it may take long time to compute codebook
savename_codebook = [path, 'CodeBook_', param.codebook, '_', param.hist, '.mat'];
save(savename_codebook, 'CodeBook');
% load(savename_codebook, 'CodeBook'); 


%==========================================================================
%% Code Vector (Coefficient) construction and histogram construction
fprintf('Code vector computing begins......\n');
fprintf('There are totally %d samples.\n', param.nsample);
for file_n = 1 : param.nsample
    fprintf('Calculating the code vector of sample %d......\n', file_n);
    
    descripMtx_cur = conf_dict(conf_dict(:,end)==file_n,2:end-1)'; % descripMtx is a P \times N matrix
    CodeVector(:,file_n) = fun_code_vector(descripMtx_cur, CodeBook, param);
end

conflict_bow = CodeVector; % BoW-based conflict representation

savename_codebook = [path, 'CodeVector_', param.codebook, '_', param.hist, '.mat'];
save(savename_codebook, 'CodeVector');



