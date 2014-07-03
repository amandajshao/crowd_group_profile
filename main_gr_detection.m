%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% main function of group detection %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% May 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.
clc;clear;close all

%% group detection parameter setting
path = '.\';
file_name = '1_8_groupSplit-festivalwalk_1_2-1'; 
load([path, file_name, '\trks_1_smooth.mat'], 'trks'); % load trks
path_img = [path, file_name, '\'];
frames = dir([path_img, '0*.jpg']);

%% parameters
paraSet.start_fr = 1;
paraSet.len_fr = 0; % predefined length of frames for group detection, if 0, no predefinition
% parameter of CF (coherent filtering)
paraSet.d = 3;
paraSet.K = 10;
paraSet.lamda = 0.6;
% parameter of Kalman Filter
paraSet.max_iter = 10;
paraSet.t_intv_kf = 10;
% parameter of group detection
paraSet.err_th = 5; % threshold for group refinement
paraSet.t_intv = 10; % collect t_intv frame info for group detection
paraSet.dv_th_seed = 0.985; % threshold for seeding group

%% detect dynamic group
fprintf('Group detection for [%s].\n', file_name);
[CF_init, group, group_id, color_ind, A, t_seq, trkClusterTimeLine, trks] = cohGr_det(trks, path_img, frames, paraSet);
% visualization
fun_grVisual(group, path_img, frames, t_seq{1,1}, t_seq{1,end}, color_ind, t_seq)
% save group detection result for group descriptor implement
save([path,file_name,'\group_1_',file_name,'.mat'], 'group');
save([path,file_name,'\group_id_1_',file_name,'.mat'], 'group_id');
save([path,file_name,'\color_1_',file_name,'.mat'], 'color_ind');
save([path,file_name,'\A_1_',file_name,'.mat'], 'A');
save([path,file_name,'\t_seq_1_',file_name,'.mat'], 't_seq');
save([path,file_name,'\trkClusterTimeLine_1_',file_name,'.mat'], 'trkClusterTimeLine');
save([path,file_name,'\trks_',file_name,'.mat'], 'trks');











