%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% main function for group size %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

% clc;clear;close all

%% Descriptor -- Group size 
path = '.\';
file_name = '1_8_groupSplit-festivalwalk_1_2-1';
path_xls = [path, 'video_info_t0.xls'];
[~,~,xls] = xlsread(path_xls);
fprintf('Group descriptor "GroupSize" for [%s].\n', file_name);

%% load collective result from group detection
load(['.\', file_name, '\trkClusterTimeLine_1_', file_name, '.mat'], 'trkClusterTimeLine');
load(['.\', file_name, '\trks_', file_name, '.mat'], 'trks');
load(['.\', file_name, '\A_1_', file_name, '.mat'], 'A');
load(['.\', file_name, '\color_1_', file_name, '.mat'], 'color_ind')

%% initialization and parameter setting
hist_bin = [0:0.2:1];

%%
trkClusterNumTime = max(trkClusterTimeLine);
[trkTime, lenTime, nTrks, trkTimeLine] = fun_trkInfo(trks);
t_seq = find(trkClusterNumTime ~= 0);

%% Do not need too long time (can be tuned)
loca = cellfun(@findstr, xls(:,1), repmat({file_name}, size(xls(:,1))), 'UniformOutput', false);
[t_loc, ~, ~] = find(~cellfun(@isempty, loca) == 1);
t_start = fun_cell2num(xls(t_loc,5));
t_end = min(t_seq(end),fun_cell2num(xls(t_loc,6)));

%% group size computation
gr_size = [];
for curTime = t_start : t_end
    t_oo = find(t_seq==curTime);
    if isempty(t_oo)
        continue;
    end
    cur_color_ind = color_ind{t_oo};
    
    % prepare data
    cur_trk_ind = find(trkClusterTimeLine(:,curTime)~=0);
    cur_gr_ind = trkClusterTimeLine(cur_trk_ind,curTime);
    data = fun_curX(trks, nTrks, trkTime, curTime, cur_trk_ind);
    
    % preprocess data
    [cur_trk_ind, cur_gr_ind, data] = fun_curX_preprocess(data, cur_gr_ind, cur_trk_ind);
    clusterValue = unique(cur_gr_ind);
    
    for grSele = 1 : length(cur_color_ind)
        clusterV = cur_color_ind(grSele);
        ind = find(cur_gr_ind == clusterV);
        subdata = data(ind,:);
        gr_size{clusterV,curTime} = size(subdata,1);
    end
end

%% record
group_size_mean = (sum(cellfun(@sum, gr_size),2))./(sum(~cellfun(@isempty,gr_size),2));
group_size_max = []; group_size_sum = []; group_size_mean_all = [];
for gr_n = 1 :length(group_size_mean)
    group_size_max = [group_size_max; group_size_mean(gr_n)./max(group_size_mean)];
    group_size_sum = [group_size_sum; group_size_mean(gr_n)./sum(group_size_mean)];
end
% group descriptor
group_size.gr_size_max_gr = group_size_max;
group_size.gr_size_sum_gr = group_size_sum;
% video descriptor
group_size.gr_size_max_v = hist(group_size_max, hist_bin);
group_size.gr_size_sum_v = hist(group_size_sum, hist_bin);

fprintf('Done!\n');

