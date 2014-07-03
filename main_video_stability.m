%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% main function for stability %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.
clc;clear;close all

%% Descriptor -- stability (3 cues)
path = '.\';
path_gr = [path,'result_groupDet_new\'];
path_xls = [path, 'video_info_t0.xls'];
[~,~,xls] = xlsread(path_xls);
path_img_dir = xls(2:end,1);

% initialization and parameter setting
group_size_th = 25; 
K = 10;
invar_th = 0.8;
hist_bin = [0:0.005:1];

for file_n = 1 : length(path_img_dir)
    file_name = path_img_dir{file_n};
    fprintf('Group descriptor "Stability" for [%d:%s].\n', file_n, file_name);
        
    %% load collective result from group detection
    load([path_gr, '\trkClusterTimeLine_1_', file_name, '.mat'], 'trkClusterTimeLine');
    load([path_gr, '\trks_', file_name, '.mat'], 'trks');
    load([path_gr, '\A_1_', file_name, '.mat'], 'A');
    load([path_gr, '\color_1_', file_name, '.mat'], 'color_ind')
    
    %% 
    trkClusterNumTime = max(trkClusterTimeLine);
    [trkTime, ~, nTrks, ~] = fun_trkInfo(trks);
    t_seq = find(trkClusterNumTime ~= 0);
    
    %% Do not need too long time (can be tuned)
    loca = cellfun(@findstr, xls(:,1), repmat({file_name}, size(xls(:,1))), 'UniformOutput', false);
    [t_loc, ~, ~] = find(~cellfun(@isempty, loca) == 1);
    t_start = fun_cell2num(xls(t_loc,5));
    t_end = min(t_seq(end),fun_cell2num(xls(t_loc,6)));

    %% initial
    cur_trk_ind = find(trkClusterTimeLine(:,t_start)~=0);
    cur_gr_ind = trkClusterTimeLine(cur_trk_ind,t_start);
    clusterV = unique(cur_gr_ind);
    
    %% stability factors computation
    %----------------------------------------------------------------------------------%
    % stat_trk_ind -- [1st_col:Time; 2nd_col:Trk_index; 3rd~17th_col:NN_Trk_index]
    %                  (rank by "Time")
    % sort_stat_trk_ind -- [1st_col:Time; 2nd_col:Trk_index; 3rd~17th_col:NN_Trk_index]
    %                  (rank by "Trk_index")
    %----------------------------------------------------------------------------------%
    
    gr_num = 1;
    group_size = {}; stab_invKnnNum_gr = []; stab_rankKnn_gr = []; stab_rwHist_gr = [];
    for grSele = 1 : length(clusterV)
        fprintf('Group %d.\n', grSele);
        stat_nn_inv_num = []; hist_sim = []; stab_rw = [];
        
        % compute random_walk:curSubdataTrkWAll; each group member's KNN at each time
        [stat_trk_ind, stat_dist, curSubdataTrkWAll, group_size_cur] = ...
            fun_stab_calc(trks, t_start, t_end, nTrks, trkTime, clusterV(grSele), K, trkClusterTimeLine);
        if isempty(stat_trk_ind) || length(unique(stat_trk_ind(:,1)))==1
            continue;
        end
        
        % initial set
        [sort_stat_trk_ind, ~] = sortrows(stat_trk_ind, 2);
        trk_ind_all = unique(sort_stat_trk_ind(:,2));
        
        num = 1;
        for i = 1 : length(trk_ind_all)
            stat_nn_num = [];
            nn_index = find(sort_stat_trk_ind(:,2)==trk_ind_all(i));
            if length(nn_index) <= 2
                num = num - 1;
            else
                nn_all = sort_stat_trk_ind(nn_index, 3:end); % all KNN indice of i_temp th trk
                nn_ind_all = unique(nn_all); % unique KNN indice of i_temp th trk
                
                %% stat1: invariant NN index
                for j = 1 : length(nn_ind_all)
                    [r_temp, c_temp] = find(nn_all == nn_ind_all(j));
                    stat_nn_num(1,j) = length(r_temp); % invariant neighbors
                end
                stat_nn_inv_seq = (stat_nn_num >= size(nn_all,1)*invar_th);
                stat_nn_inv_num(num) = length(find(stat_nn_inv_seq ~= 0));
                
                %% stat2: ranked Knn index similarity
                nn_sim = zeros(1, size(nn_all,1)-1);
                for j = 2 : size(nn_all,1)
                    nn_sim(1,j-1) = edit_distance_levenshtein(nn_all(1,:)', nn_all(j,:)');
                end
                hist_sim(num,:) = hist(nn_sim, 1:K);
            end
            num = num + 1;
        end
        
        %% stat3: group member relative stability: graph-based (random walk) over invariant KNN
        [trkWAll, stab_rw] = fun_stab_rw(sort_stat_trk_ind, curSubdataTrkWAll);
        
        if num ~= 1
            stab_invKnnNum_gr(gr_num, :) = sum(stat_nn_inv_num)/length(stat_nn_inv_num); % feature of invariant Knn number
            stab_rankKnn_gr(gr_num, :) = sum(hist_sim, 1)/size(hist_sim,1); % feature of rankKnn variance
            stab_rwHist_gr(gr_num,:) = hist(stab_rw,hist_bin); % kl-divergence of random walk
            group_size = [group_size; group_size_cur];
            gr_num = gr_num + 1;
        end
    end
    
    %% record: only groupSize larger than threshold
    group_size_mean = (sum(cellfun(@sum, group_size),2))./(sum(~cellfun(@isempty,group_size),2));
    stab_invKnnNum_all = []; stab_rankKnn_all = []; stab_rwHist_all = [];
    group_size_all = [];
    for gr_n = 1 :length(group_size_mean)
        if group_size_mean(gr_n) > group_size_th
            stab_invKnnNum_all =  [stab_invKnnNum_all; stab_invKnnNum_gr(gr_n,:)];
            stab_rankKnn_all = [stab_rankKnn_all; stab_rankKnn_gr(gr_n,:)];
            stab_rwHist_all = [stab_rwHist_all; stab_rwHist_gr(gr_n,:)];
            group_size_all = [group_size_all; group_size_mean(gr_n)];
        end
    end
    stability(file_n).stab_invKnnNum = mean(stab_invKnnNum_all);
    stability(file_n).stab_rankKnn = sum(stab_rankKnn_all,1)./size(stab_rankKnn_all,1);
    stability(file_n).stab_rwHist = sum(stab_rwHist_all,1)./size(stab_rwHist_all,1);
    % other records
    stability(file_n).stab_all.stab_invKnnNum_all = stab_invKnnNum_all;
    stability(file_n).stab_all.stab_rankKnn_all = stab_rankKnn_all;
    stability(file_n).stab_all.stab_rwHist_all = stab_rwHist_all;
    stability(file_n).group_size = group_size_all;
    stability(file_n).file_name = file_name;
    
    fprintf('Done!\n');
end

path_result = [path, 'result_groupDescr_new\'];
mkdir(path_result);
save([path_result,'stability.mat'], 'stability')







