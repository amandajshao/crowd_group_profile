%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% main function for uniformity %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.
clc;clear;close all

%% Descriptor -- Uniformity (the number of graph-cuts)
path = '.\';
path_gr = [path,'result_groupDet_new\'];
path_xls = [path, 'video_info_t0.xls'];
[~,~,xls] = xlsread(path_xls);
path_img_dir = xls(2:end,1);

% initialization and parameter setting
group_size_th = 25; 

for file_n = 1 : length(path_img_dir)
    file_name = path_img_dir{file_n};
    fprintf('Group descriptor "Uniformity" for [%d:%s].\n', file_n, file_name);
    
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
    
    %%
    group_size = []; clusterNum_record = [];
    for curTime = t_start : t_end
        fprintf('Frame %d.\n', curTime);
        
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
    
    %% record: only groupSize larger than threshold
    group_size_mean = (sum(cellfun(@sum, group_size),2))./(sum(~cellfun(@isempty,group_size),2));
    unif_all = []; unif_mean = []; unif_var = [];
    group_size_all = [];
    for gr_n = 1 :length(group_size_mean)
        if group_size_mean(gr_n) > group_size_th
            cluster_temp = cellfun(@sum,clusterNum_record(gr_n,:));
            unif_all =  [unif_all; cluster_temp];
            unif_mean = [unif_mean; sum(cluster_temp)./length(find(cluster_temp~=0))];
            unif_var = [unif_var; std(cluster_temp(cluster_temp~=0))];
            group_size_all = [group_size_all; group_size_mean(gr_n)];
        end
    end
    uniform(file_n).unif_mean = mean(unif_mean);
    uniform(file_n).unif_var = mean(unif_var);
    % other records
    uniform(file_n).unif_all = unif_all;
    uniform(file_n).group_size = group_size_all;
    uniform(file_n).file_name = file_name;   
end

path_result = [path, 'result_groupDescr_new\'];
mkdir(path_result);
save([path_result,'v_uniform_fill.mat'], 'uniform')

