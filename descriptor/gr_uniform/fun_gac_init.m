function [K,distance_matrix,clusterNumber] = fun_gac_init(X, percentage, init_K)

% FUN_GAC_INIT: Summary of this function goes here
%               Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

% initial parameters
K = init_K;
distance_matrix = distmat(X);

[~, ~, ~, ~, rk] = graphAgglomerativeClustering(distance_matrix, ...
    1, 'gdl', K, 0, 1, false, X);
[~, topidx] = sort(rk, 'descend');

% NN indices
N = size(distance_matrix,1);
chosen_idx = topidx(1:round(N*percentage));
distance_matrix_sub = distance_matrix(chosen_idx, chosen_idx);
X_sub = X(chosen_idx, :);

[~, Q, ~, ~, ~] = graphAgglomerativeClustering(distance_matrix_sub, ...
    1, 'gdl', K, 1, 1, false, X_sub, true);

Q_diff = diff(Q);
if isempty(Q)
    clusterNumber = 1;
elseif length(Q_diff) == 1
    clusterNumber = 1;
else
    Q_evol = diff(diff(Q));
    [~, clusterNumPred] = min(Q_evol);
    clusterNumber = length(Q)-clusterNumPred+1;
end

end
