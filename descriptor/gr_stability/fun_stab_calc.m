function [stat_trk_ind, stat_dist, curSubdataTrkWAll, group_size] = ...
    fun_stab_calc(trks, t_0, t_1, nTrks, trkTime, cur_clusterV, K, trkClusterTimeLine)
% FUN_STAB_CALC: Summary of this function goes here
%                Detailed explanation goes here
% --------------------------------------------------------------------- %
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

stat_trk_ind = []; % trk index set
stat_dist = [];    % distance set
curSubdataTrkWAll = [];

group_size = cell(1,t_1);
for curTime = t_0 : t_1
    t_oo = curTime - t_0 + 1;
    
    % prepare data
    cur_trk_ind = find(trkClusterTimeLine(:,curTime)~=0);
    cur_gr_ind = trkClusterTimeLine(cur_trk_ind,curTime);
    data = fun_curX(trks, nTrks, trkTime, curTime, cur_trk_ind);
    
    % preprocess data
    [cur_trk_ind, cur_gr_ind, data] = fun_curX_preprocess(data, cur_gr_ind, cur_trk_ind);
    ind = find(cur_gr_ind==cur_clusterV);
    trk_mem = cur_trk_ind(ind);
    if isempty(trk_mem)
        cur_trk_ind = [];
        cur_gr_ind = [];
        data = [];
        break;
    end
    
    subdata = data(ind,:);
    group_size{1,curTime} = size(subdata,1);
    
    if size(subdata,1) > K % if delete this condition, when size(subdata,1)<=K, gacMink.m will cause error
        % compute sorted distances 
        dist_data = distmat(subdata); % distance matrix
        [sortedDist,NNIndex] = gacMink(dist_data, K+1, 2); % sorted distance and corresponding KNN index
                
        % random walk over this group, W_norm is the random walk matrix
        [W_norm, adjacency] = fun_rw(NNIndex); 
        curSubdataTrkWAll{t_oo, 1} = [[0;trk_mem(NNIndex(:,1))], [trk_mem(NNIndex(:,1))';W_norm]];
        curSubdataTrkWAll{t_oo, 2} = adjacency;
        curSubdataTrkWAll{t_oo, 3} = trk_mem(NNIndex);
        curSubdataTrkWAll{t_oo, 4} = subdata;
                        
        stat_trk_ind = [stat_trk_ind; [ones(size(trk_mem))*t_oo, reshape(trk_mem(NNIndex), size(NNIndex))]];
        stat_dist = [stat_dist; [ones(size(trk_mem))*t_oo sortedDist]];
    end
end

end

