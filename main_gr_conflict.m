%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% main function for conflict %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

% clc;clear;close all

%% Descriptor -- shape context based conftions distribution
% TIPs: after computing the conflict descriptors for all videos, use
% main_BoW.m under folder of '.\descriptor\gr_conflict\' to construct 
% codebook and codevector, finally get the BoW-based conflict descriptor
path = '.\';
file_name = '1_8_groupSplit-festivalwalk_1_2-1';
path_img = [path, file_name, '\'];
path_xls = [path, 'video_info_t0.xls'];
[~,~,xls] = xlsread(path_xls);
fprintf('Group descriptor "Conflict" for [%s].\n', file_name);

%% load collective result from group detection
load(['.\', file_name, '\trkClusterTimeLine_1_', file_name, '.mat'], 'trkClusterTimeLine');
load(['.\', file_name, '\trks_', file_name, '.mat'], 'trks');
load(['.\', file_name, '\A_1_', file_name, '.mat'], 'A');
load(['.\', file_name, '\color_1_', file_name, '.mat'], 'color_ind');

%% initialization and parameter setting
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

gv_flag = 1; % 1 for group conflict calculation; 0 for video conflict calculation

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
cur_trk_ind = find(trkClusterTimeLine(:,t_start)~=0);
cur_gr_ind = trkClusterTimeLine(cur_trk_ind,t_start);
data = fun_curX(trks, nTrks, trkTime, t_start, cur_trk_ind);
[cur_trk_ind, cur_gr_ind, data] = fun_curX_preprocess(data, cur_gr_ind, cur_trk_ind);
clusterValue = unique(cur_gr_ind);

%% conftion computation
conf_dict_data = [];
for grSele = 1 : length(clusterValue)
    fprintf('Group %d.\n', grSele);
    
    %% compute conflict value
    % groupPoint: [curTime x y v_x v_y confValue confNum]
    [groupPoint, group_size_cur] = fun_conflict_val(trkClusterTimeLine, A, trks, color_ind, ...
        t_start, t_end, t_seq, nTrks, trkTime, K_percent, clusterValue, ...
        grSele, gen_mem_time_len, fit_time_len, conf_norm);
    
    %% compute conflict distribution with calibration and shape context computation    
    conf = fun_conflict_sc(groupPoint, cali_para, sc_para, gv_flag);
    
    conf_dict_data = [conf_dict_data; [grSele*ones(size(conf.sc_conf_contour,1),1), conf.sc_conf_contour]];
    conf_all{1,grSele} = conf;
end

%% record
conflict.conf = conf_all;
conflict.conf_dict = conf_dict_data;

fprintf('Done!\n');









