%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% main function for uniformity %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

% clc;clear;close all

%% Descriptor -- Uniformity (the number of graph-cuts)
path = '.\';
file_name = '1_8_groupSplit-festivalwalk_1_2-1';
path_img = [path, file_name, '\'];
path_xls = [path, 'video_info_t0.xls'];
[~,~,xls] = xlsread(path_xls);
fprintf('Group descriptor "Uniformity" for [%s].\n', file_name);

%% load collective result from group detection
load(['.\', file_name, '\trkClusterTimeLine_1_', file_name, '.mat'], 'trkClusterTimeLine');
load(['.\', file_name, '\trks_', file_name, '.mat'], 'trks');
load(['.\', file_name, '\A_1_', file_name, '.mat'], 'A');
load(['.\', file_name, '\color_1_', file_name, '.mat'], 'color_ind');

%%
trkClusterNumTime = max(trkClusterTimeLine);
[trkTime, ~, nTrks, ~] = fun_trkInfo(trks);
t_seq = find(trkClusterNumTime ~= 0);

%% Do not need too long time (can be tuned)
loca = cellfun(@findstr, xls(:,1), repmat({file_name}, size(xls(:,1))), 'UniformOutput', false);
[t_loc, ~, ~] = find(~cellfun(@isempty, loca) == 1);
t_start = fun_cell2num(xls(t_loc,5));
t_end = min(t_seq(end),fun_cell2num(xls(t_loc,6)));

%%
for curTime = t_start : t_end    
    % prepare data
    cur_trk_ind = find(trkClusterTimeLine(:,curTime)~=0);
    cur_gr_ind = trkClusterTimeLine(cur_trk_ind,curTime);
    clusterValue = unique(cur_gr_ind);
    data = fun_curX(trks, nTrks, trkTime, curTime, cur_trk_ind);
    
    % preprocess data
    [cur_trk_ind, cur_gr_ind, data] = fun_curX_preprocess(data, cur_gr_ind, cur_trk_ind);
    
    % clusterNum calculation
    for grSele = 1 : length(clusterValue)
        clusterV = clusterValue(grSele);
        ind = find(cur_gr_ind==clusterV);
        trk_mem = cur_trk_ind(ind);
        
        subdata = data(ind,:);
        group_size{clusterV,curTime} = size(subdata,1);
        if isempty(trk_mem)
            cur_trk_ind = [];
            cur_gr_ind = [];
            data = [];
            break;
        end
        
        if size(subdata,1) <= 5
            data_sub = [subdata [1:size(subdata,1)]'];
            clusteredLabels = ones(size(subdata,1),1);
            percentage = 1;
            clusterNum = 1;
        else
            K = max(2,floor(length(subdata(:,1))/10));
            distance_matrix = distmat(subdata);
            [~, distance_matrix, clusterNum] = fun_gac_init(subdata, 1, K);
            if isempty(clusterNum)
                data_sub = [subdata [1:size(subdata,1)]'];
                clusteredLabels = ones(size(subdata,1),1);
                percentage = 1;
                clusterNum = 1;
            end
        end
        clusterNum_record{clusterV,curTime} = clusterNum;
    end
end

%% record
uniform.unif_gr_mean = (sum(cellfun(@sum, clusterNum_record),2))./(sum(~cellfun(@isempty,clusterNum_record),2));
temp = cellfun(@sum, clusterNum_record);
for i = 1 : size(clusterNum_record,1)    
    uniform.unif_gr_var(i,1) = std(temp(i,temp(i,:)~=0));
end
uniform.unif_v_mean = mean(uniform.unif_gr_mean);
uniform.unif_v_var = mean(uniform.unif_gr_var);

fprintf('Done!\n');



