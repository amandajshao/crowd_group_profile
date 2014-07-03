%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% main function of group detection %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% June 15, 2014, Jing Shao
% If you use this code, please cite the paper:
% J. Shao, C. C. Loy, X. Wang, "Scene-Independent Group Profiling in Crowd", CVPR, 2014.

clc;clear;close all

%% Please make sure add three folders and their subfolders to the path before running code
% They are: CoheGroupDetect, descriptor, util
% CoheGroupDetect contains group detection code, except main function
% descriptor contains descriptor code, except main function
% util contains some utilized code

%% group detection
main_gr_detection; 

%% group descriptor
main_gr_collectiveness;
main_gr_stability;
main_gr_uniform;
main_gr_conflict;
main_gr_size;




