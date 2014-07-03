%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% main function for collectiveness %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

% clc;clear;close all

%% Descriptor -- Collectiveness (A-fit error within each group)
path = '.\';
file_name = '1_8_groupSplit-festivalwalk_1_2-1';
path_img = [path, file_name, '\'];
path_xls = [path, 'video_info_t0.xls'];
[~,~,xls] = xlsread(path_xls);
fprintf('Group descriptor "Collectiveness" for [%s].\n', file_name);

%% load collective result from group detection
load(['.\', file_name, '\trkClusterTimeLine_1_', file_name, '.mat'], 'trkClusterTimeLine');
load(['.\', file_name, '\trks_', file_name, '.mat'], 'trks');
load(['.\', file_name, '\A_1_', file_name, '.mat'], 'A');
load(['.\', file_name, '\color_1_', file_name, '.mat'], 'color_ind')

%% initialization and parameter setting
fit_time_len = 10; % the same as the coherent group detection's parameter

%%
trkClusterNumTime = max(trkClusterTimeLine);
[trkTime, lenTime, nTrks, trkTimeLine] = fun_trkInfo(trks);
t_seq = find(trkClusterNumTime ~= 0);

%% Do not need too long time (can be tuned)
loca = cellfun(@findstr, xls(:,1), repmat({file_name}, size(xls(:,1))), 'UniformOutput', false);
[t_loc, ~, ~] = find(~cellfun(@isempty, loca) == 1);
t_start = fun_cell2num(xls(t_loc,5));
t_end = min(t_seq(end),fun_cell2num(xls(t_loc,6)));

%% 
cur_trk_ind = find(trkClusterTimeLine(:,t_start)~=0);
cur_gr_ind = trkClusterTimeLine(cur_trk_ind,t_start);
data = fun_curX(trks, nTrks, trkTime, t_start, cur_trk_ind);
[cur_trk_ind, cur_gr_ind, data] = fun_curX_preprocess(data, cur_gr_ind, cur_trk_ind);
clusterValue = unique(cur_gr_ind);

%% collectiveness computation
rmse_trk_mem = []; 
for grSele = 1 : length(clusterValue)
    clusterV = clusterValue(grSele);
    for curTime = t_start : t_end
        t_oo = find(t_seq==curTime);
        if isempty(t_oo)
            continue;
        end
        A_cur = A(:,t_oo);
        cur_color_ind = color_ind{t_oo};
        cur_color_ind_cluster = find(cur_color_ind==clusterV,1);
        
        % prepare data
        cur_trk_ind = find(trkClusterTimeLine(:,curTime)~=0);
        cur_gr_ind = trkClusterTimeLine(cur_trk_ind,curTime);
        data = fun_curX(trks, nTrks, trkTime, curTime, cur_trk_ind);
        
        % preprocess data
        [cur_trk_ind, cur_gr_ind, data] = fun_curX_preprocess(data, cur_gr_ind, cur_trk_ind);
        ind = find(cur_gr_ind==clusterValue(grSele));
        trk_mem = cur_trk_ind(ind);
        if isempty(trk_mem)
            cur_trk_ind = [];
            cur_gr_ind = [];
            data = [];
            break;
        end
        subdata = data(ind,:);
        
        %%
        trk_mem_cur = [];
        if ~isempty(cur_color_ind_cluster)
            A_cur_cluster = A_cur{cur_color_ind_cluster, 1};
            rmse_trk_mem_cur = zeros(size(ind,1),1);
            for i = 1 : size(ind,1)
                trk_mem = trks(cur_trk_ind(ind(i)));
                trk_mem_t = find(trk_mem.t >= curTime & trk_mem.t <= curTime+fit_time_len-1);
                trk_mem_cur = [trk_mem.x(trk_mem_t)'; trk_mem.y(trk_mem_t)'; ones(1, length(trk_mem_t))];
                
                % use A_cur_cluster to generate trk_mem_cur_gen
                trk_mem_cur_gen = ones(3, size(trk_mem_cur,2));
                trk_mem_cur_gen(:,1) = trk_mem_cur(:,1);
                for fit_t = 2 : size(trk_mem_cur,2)
                    trk_mem_cur_gen(:,fit_t) = A_cur_cluster * trk_mem_cur_gen(:,fit_t-1);
                end
                
                % compute fit error -- RMSE
                rmse_trk_mem_cur(i,1) = sqrt(sum(sum((trk_mem_cur_gen-trk_mem_cur).^2)) / size(trk_mem_cur,2));
            end
            
            if isempty(ind)
                rmse_trk_mem{clusterV,curTime} = [];
            else
                rmse_trk_mem{clusterV,curTime} = sum(rmse_trk_mem_cur)/length(rmse_trk_mem_cur);
            end
        end
    end
end
rmse_trk_mem_mean = (sum(cellfun(@sum, rmse_trk_mem),2))./(sum(~cellfun(@isempty,rmse_trk_mem),2));

%% record
% group descriptor
collectiveness.coll_gr = rmse_trk_mem_mean;
% video descriptor
collectiveness.coll_v_mean = mean(rmse_trk_mem_mean);
collectiveness.coll_v_var = std(rmse_trk_mem_mean);

fprintf('Done!\n');







