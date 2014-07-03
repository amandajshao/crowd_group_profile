%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% main function for conflict %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.
clc;clear;close all

%% Descriptor -- shape context based conftions distribution
% TIPs: after computing the conflict descriptors for all videos, use
% main_BoW.m under folder of '.\descriptor\gr_conflict\' to construct
% codebook and codevector, finally get the BoW-based conflict descriptor
path = '.\';
path_gr = [path,'result_groupDet_new\'];
path_xls = [path, 'video_info_t0.xls'];
[~,~,xls] = xlsread(path_xls);
path_img_dir = xls(2:end,1);

% initialization and parameter setting
K_percent = 1/3;
conf_norm = 0; % 1--normalized

fit_time_len = 10; %the same as the coherent group detection's parameter & collectiveness
gen_mem_time_len = 2;

cali_para.smooth_size = 10;
cali_para.smooth_sigma = 10;
cali_para.pattern_scale = 100;
cali_para.contour_th = 0.05;
sc_para.nbins_theta = 12;
sc_para.nbins_r = 5;
sc_para.r_inner = 0.125;
sc_para.r_outer = 4;
sc_para.mean_dist_global = [];

group_size_th = 25; % 10
gr_num = 0;
gv_flag = 0; % 1 for group conflict calculation; 0 for video conflict calculation

conf_dict_data = []; 
for file_n = 1 : length(path_img_dir)
    file_name = path_img_dir{file_n};
    fprintf('Group descriptor "Conflict" for [%d:%s].\n', file_n, file_name);
    
    %% load collective result from group detection
    load([path_gr, '\trkClusterTimeLine_1_', file_name, '.mat'], 'trkClusterTimeLine');
    load([path_gr, '\trks_', file_name, '.mat'], 'trks');
    load([path_gr, '\A_1_', file_name, '.mat'], 'A');
    load([path_gr, '\color_1_', file_name, '.mat'], 'color_ind')
    
    %%
    trkClusterNumTime = max(trkClusterTimeLine);
    t_seq = find(trkClusterNumTime ~= 0);
    [trkTime, ~, nTrks, ~] = fun_trkInfo(trks);

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
    
    %%  conftion computation
    groupPointAll = []; group_size = [];
    for grSele = 1 : length(clusterValue)
        fprintf('Group %d.\n', grSele);
        
        %% compute conflict value
        % groupPoint: [curTime x y v_x v_y confValue confNum]
        [groupPoint, group_size_cur] = fun_conflict_val(trkClusterTimeLine, A, trks, color_ind, ...
            t_start, t_end, t_seq, nTrks, trkTime, K_percent, clusterValue, ...
            grSele, gen_mem_time_len, fit_time_len, conf_norm);
        group_size_mean = (sum(cellfun(@sum, group_size_cur),2))./(sum(~cellfun(@isempty,group_size_cur),2));
        if group_size_mean > group_size_th
            group_size = [group_size; group_size_mean];
            groupPointAll = [groupPointAll; groupPoint];
            gr_num = gr_num + 1;
        end      
    end
        
    %% compute conflict distribution with calibration and shape context computation
    conf = fun_conflict_sc(groupPointAll, cali_para, sc_para, gv_flag);

    %% record
    conf_dict_data = [conf_dict_data; [gr_num*ones(size(conf.sc_conf_contour,1),1), conf.sc_conf_contour], repmat(file_n-2,size(conf.sc_conf_contour,1),1)];
    conflict(file_n).conf = conf;
    conflict(file_n).group_size = group_size;
end

path_result = [path, 'result_groupDescr_new\'];
mkdir(path_result);
save([path_result,'conf_dict.mat'], 'conf_dict_data')
save([path_result,'conflict.mat'], 'conflict')






