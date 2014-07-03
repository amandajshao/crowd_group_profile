function [W_norm, S] = fun_rw(NNIndex)

% FUN_RW: Random walk
%         Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

%% 
% S -- Adjacency matrix
data_temp1 = [];
data_temp2 = [];
s = ones(1,(size(NNIndex,2)-1)* size(NNIndex,1));
for data_n = 1 : size(NNIndex,1)
    data_temp1 = [data_temp1, NNIndex(data_n,1) * ones(1,size(NNIndex,2)-1)];
    data_temp2 = [data_temp2, NNIndex(data_n,2:end)];
end
S_sparse = sparse(data_temp1, data_temp2, s, size(NNIndex,1), size(NNIndex,1));
S = full(S_sparse);

% C -- Connection matrix
c1 = tril(S); 
c2 = triu(S); 
C = c1 + c2'; 
C = C + C'; 

% Transition matrix
S_sum = sum(S,2);
S_sum = repmat(S_sum,1,size(S,2));
T = S./S_sum;

% random walk
alpha = 0.9; % can be tuned
W = pinv(eye(size(T)) - alpha*T) - eye(size(T)); % T+T^2+T^3+...
W_sum = sum(W,2);
W_sum = repmat(W_sum,1,size(W,2));
W_norm = W./W_sum;

% W_sum_all = sum(sum(W))/(size(W,1)^2);





