function CodeBook = fun_codebook_computation(param, descripMtx)

% FUN_CODEBOOK_COMPUTATION: Compute codebook
%                           Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

if strcmp(param.codebook, 'sparse_coding')
    func       = param.sc.func;
    num_basis  = param.num_basis;
    num_trials = param.sc.trials;
    beta       = param.sc.beta;
    batch_size = param.sc.batch;
    [CodeBook U stat] = sparse_coding(descripMtx, num_basis, beta, func, [], num_trials, batch_size);
elseif strcmp(param.codebook, 'kmeans')
    num_basis  = param.num_basis;
    [CodeBook, ~] = litekmeans(descripMtx, num_basis);
    param.num_basis = size(CodeBook, 2); 
end

end